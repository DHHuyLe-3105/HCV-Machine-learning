"""
Firth penalized logistic regression (Firth 1993; Heinze & Schemper 2002).

Penalised log-likelihood:  l*(b) = l(b) + 0.5 * log|I(b)|
Modified score:            U*(b_j) = sum_i [ y_i - p_i + h_i (0.5 - p_i) ] x_ij
where h_i = diag of the hat matrix H = W^.5 X (X'WX)^-1 X' W^.5,  W = diag(p(1-p)).

Inference is by penalised likelihood ratio test and profile penalised
likelihood confidence intervals (NOT Wald), which is the point of using
Firth in the first place: with complete separation the Wald SE explodes.
"""
import numpy as np
from scipy import optimize, stats


def _pen_loglik(beta, X, y):
    eta = X @ beta
    # log-likelihood, computed stably
    ll = np.sum(y * eta - np.logaddexp(0.0, eta))
    p = 1.0 / (1.0 + np.exp(-eta))
    w = p * (1.0 - p)
    XtWX = X.T @ (X * w[:, None])
    sign, logdet = np.linalg.slogdet(XtWX)
    if sign <= 0:
        return -np.inf
    return ll + 0.5 * logdet


def _fit(X, y, fixed_idx=None, fixed_val=None, max_iter=500, tol=1e-10):
    """Newton-Raphson on the modified score. If fixed_idx is given, that
    coefficient is held at fixed_val (used for profiling)."""
    n, k = X.shape
    beta = np.zeros(k)
    if fixed_idx is not None:
        beta[fixed_idx] = fixed_val
    free = [j for j in range(k) if j != fixed_idx]
    for _ in range(max_iter):
        eta = X @ beta
        p = 1.0 / (1.0 + np.exp(-eta))
        w = p * (1.0 - p)
        w = np.clip(w, 1e-12, None)
        XtWX = X.T @ (X * w[:, None])
        XtWX_inv = np.linalg.pinv(XtWX)
        # hat diagonal
        Xw = X * np.sqrt(w)[:, None]
        h = np.einsum('ij,jk,ik->i', Xw, XtWX_inv, Xw)
        u = X.T @ (y - p + h * (0.5 - p))
        if fixed_idx is None:
            step = XtWX_inv @ u
        else:
            sub = np.ix_(free, free)
            step_free = np.linalg.pinv(XtWX[sub]) @ u[free]
            step = np.zeros(k)
            step[free] = step_free
        # step halving on the penalised log-likelihood
        base = _pen_loglik(beta, X, y)
        factor, ok = 1.0, False
        for _h in range(30):
            cand = beta + factor * step
            if _pen_loglik(cand, X, y) >= base - 1e-9:
                ok = True
                break
            factor /= 2.0
        if not ok:
            break
        beta = cand
        if np.max(np.abs(factor * step)) < tol:
            break
    return beta


def _profile_ci(X, y, j, beta_hat, ll_hat, alpha=0.05):
    """Profile penalised-likelihood CI for coefficient j."""
    crit = stats.chi2.ppf(1 - alpha, 1) / 2.0

    def g(val):
        b = _fit(X, y, fixed_idx=j, fixed_val=val)
        return ll_hat - _pen_loglik(b, X, y) - crit

    out = []
    for direction in (-1, 1):
        lo, hi = beta_hat[j], beta_hat[j]
        step = 0.5
        for _ in range(200):
            hi = hi + direction * step
            if g(hi) > 0:
                break
            lo = hi
            step *= 1.5
            if abs(hi - beta_hat[j]) > 60:
                out.append(np.nan)
                break
        else:
            out.append(np.nan)
            continue
        if len(out) and np.isnan(out[-1]):
            continue
        a, b = (min(lo, hi), max(lo, hi))
        try:
            root = optimize.brentq(g, a, b, xtol=1e-6)
        except ValueError:
            root = np.nan
        out.append(root)
    return out[0], out[1]


def firth_logit(X, y, names=None, alpha=0.05, ci=True):
    """Fit Firth logistic regression. X must already include an intercept
    column. Returns a dict of arrays."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    k = X.shape[1]
    names = list(names) if names is not None else [f'x{j}' for j in range(k)]
    beta = _fit(X, y)
    ll_hat = _pen_loglik(beta, X, y)

    lo = np.full(k, np.nan)
    hi = np.full(k, np.nan)
    pval = np.full(k, np.nan)
    for j in range(k):
        # penalised likelihood ratio test of H0: beta_j = 0
        b0 = _fit(X, y, fixed_idx=j, fixed_val=0.0)
        lr = 2.0 * (ll_hat - _pen_loglik(b0, X, y))
        pval[j] = stats.chi2.sf(max(lr, 0.0), 1)
        if ci:
            lo[j], hi[j] = _profile_ci(X, y, j, beta, ll_hat, alpha)
    return {'names': names, 'beta': beta, 'or': np.exp(beta),
            'ci_low': np.exp(lo), 'ci_high': np.exp(hi), 'p': pval,
            'loglik': ll_hat, 'n': X.shape[0], 'events': int(y.sum())}


if __name__ == '__main__':
    # Validation. For a saturated model on a 2x2 table (intercept + one binary
    # predictor) the Firth penalty is exactly equivalent to adding 0.5 to every
    # cell, so the Firth OR must equal ((a+.5)(d+.5))/((b+.5)(c+.5)).
    rng = np.random.default_rng(0)
    for (a, b, c, d) in [(16, 5, 6, 37), (21, 0, 0, 38), (3, 7, 2, 40), (10, 0, 5, 5)]:
        # x=1: a events, b non-events; x=0: c events, d non-events
        yv = np.r_[np.ones(a), np.zeros(b), np.ones(c), np.zeros(d)]
        xv = np.r_[np.ones(a + b), np.zeros(c + d)]
        Xm = np.column_stack([np.ones_like(xv), xv])
        res = firth_logit(Xm, yv, ['intercept', 'x'], ci=False)
        expected = ((a + .5) * (d + .5)) / ((b + .5) * (c + .5))
        print(f'cells {a},{b},{c},{d}  firth OR={res["or"][1]:.6f}  '
              f'closed form={expected:.6f}  diff={abs(res["or"][1]-expected):.2e}')

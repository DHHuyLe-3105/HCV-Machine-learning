from __future__ import annotations

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class FiniteClipper(BaseEstimator, TransformerMixin):
    def __init__(
        self,
        *,
        lower_quantile: float = 0.001,
        upper_quantile: float = 0.999,
        fill_strategy: str = "median",
        fill_value: float = 0.0,
    ) -> None:
        self.lower_quantile = lower_quantile
        self.upper_quantile = upper_quantile
        self.fill_strategy = fill_strategy
        self.fill_value = fill_value

    def fit(self, X, y=None):
        X_arr = np.asarray(X, dtype=float)
        X_arr = np.where(np.isfinite(X_arr), X_arr, np.nan)

        if X_arr.ndim != 2:
            raise ValueError(f"Expected 2D array, got shape {X_arr.shape}.")

        if not (0.0 <= self.lower_quantile <= 1.0 and 0.0 <= self.upper_quantile <= 1.0):
            raise ValueError("Quantiles must be between 0 and 1.")
        if self.lower_quantile > self.upper_quantile:
            raise ValueError("lower_quantile must be <= upper_quantile.")

        lower = np.nanquantile(X_arr, self.lower_quantile, axis=0)
        upper = np.nanquantile(X_arr, self.upper_quantile, axis=0)

        invalid = ~np.isfinite(lower) | ~np.isfinite(upper)
        if invalid.any():
            lower = np.where(invalid, 0.0, lower)
            upper = np.where(invalid, 0.0, upper)

        if self.fill_strategy == "median":
            fill = np.nanmedian(X_arr, axis=0)
            fill = np.where(np.isfinite(fill), fill, self.fill_value)
        elif self.fill_strategy == "constant":
            fill = np.full(X_arr.shape[1], float(self.fill_value))
        else:
            raise ValueError("fill_strategy must be 'median' or 'constant'.")

        self.lower_ = lower
        self.upper_ = upper
        self.fill_ = fill
        return self

    def transform(self, X):
        X_arr = np.asarray(X, dtype=float)
        X_arr = np.where(np.isfinite(X_arr), X_arr, np.nan)

        if X_arr.ndim != 2:
            raise ValueError(f"Expected 2D array, got shape {X_arr.shape}.")

        if not hasattr(self, "lower_"):
            raise RuntimeError("FiniteClipper is not fitted.")

        X_arr = np.where(np.isnan(X_arr), self.fill_, X_arr)
        X_arr = np.clip(X_arr, self.lower_, self.upper_)
        return X_arr

"""
Run the ORIGINAL sequence-level pipeline across the same five seeds, in the same
environment, so that SA5 compares like with like.

Why this is needed: the grouped multi-seed numbers were produced on this
machine, but the sequence-level numbers currently in the manuscript
(0.888 +/- 0.030 for voting, 0.887 +/- 0.027 for stacking) came from a different
machine with different library versions. A sensitivity analysis has to be
compared against a baseline run under identical conditions.

This script changes nothing about the pipeline. It only runs the notebook's own
multi-seed function with the original splitting, then merges the result with the
grouped summary already produced and prints the comparison.

    python rerun_ungrouped_multiseed.py

Put it in the repository root next to genbank_metadata.csv. It writes into
outputs/Ungrouped_<timestamp>/ and never touches the existing folders.
"""
import json
import os
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) or '.')
import glob
from datetime import datetime

import numpy as np
import pandas as pd

NOTEBOOK = 'Predictive_model_ras.ipynb'
RAS_PATH = 'data/data_RAS.xlsx'
SEEDS = (42, 7, 13, 21, 31)          # same five seeds as the grouped run
OUTDIR = os.path.join('outputs', 'Ungrouped_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
os.makedirs(OUTDIR, exist_ok=True)

# ------------------------------------------------- load the notebook functions
with open(NOTEBOOK, encoding='utf-8') as fh:
    nb = json.load(fh)
cells = [''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code']
ns = {'__name__': '__main__'}
for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11]:
    exec(compile(cells[i], f'<cell {i}>', 'exec'), ns)

# -------------------- same rank-safe FDA as the grouped run, for comparability
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as _QDA
from sklearn.decomposition import PCA as _PCA


class SafeQDA(_QDA):
    def fit(self, X, y):
        self.pca_ = None
        try:
            return super().fit(X, y)
        except Exception as exc:
            X = np.asarray(X, dtype=float)
            counts = np.bincount(np.asarray(y).astype(int))
            for k in (20, 10, 5, 3):
                k = int(min(k, X.shape[1] - 1, counts.min() - 1))
                if k < 2:
                    continue
                for reg in (max(self.reg_param, 0.1), 0.3, 0.6):
                    try:
                        self.pca_ = _PCA(n_components=k, random_state=0).fit(X)
                        self.reg_param = reg
                        out = super().fit(self.pca_.transform(X), y)
                        print(f'    FDA: plain fit failed ({type(exc).__name__}); '
                              f'refitted with {k} components, reg_param={reg}')
                        return out
                    except Exception:
                        continue
            self.pca_ = None
            raise

    def _maybe(self, X):
        return self.pca_.transform(np.asarray(X, dtype=float)) if self.pca_ is not None else X

    def predict(self, X):
        return super().predict(self._maybe(X))

    def predict_proba(self, X):
        return super().predict_proba(self._maybe(X))

    def decision_function(self, X):
        return super().decision_function(self._maybe(X))


ns['QuadraticDiscriminantAnalysis'] = SafeQDA

import sys
import sklearn
try:
    import xgboost
    import imblearn
    print(f'versions: python {sys.version.split()[0]}  sklearn {sklearn.__version__}  '
          f'xgboost {xgboost.__version__}  imblearn {imblearn.__version__}')
except ImportError:
    print(f'versions: python {sys.version.split()[0]}  sklearn {sklearn.__version__}')

# ------------------------------------------------------------------------ run
X, y, ras_df, _, feature_names = ns['build_feature_matrix'](RAS_PATH, ras_sheet='Sheet1')
print(f'feature matrix {X.shape}; ORIGINAL sequence-level splitting, seeds {SEEDS}')

ns['multi_seed_sensitivity'](X, y, ras_df, feature_names, seeds=SEEDS,
                             label='RAS_ungrouped_multiseed', save_path=OUTDIR)

# ------------------------------------------------------------------- compare
def load_summary(path):
    d = pd.read_excel(path).set_index('Model')
    return d[['AUC_Score_mean', 'AUC_Score_std', 'AUC_Score_min', 'AUC_Score_max']]


ung = load_summary(os.path.join(OUTDIR, 'multiseed_summary_RAS_ungrouped_multiseed.xlsx'))

grouped_files = sorted(glob.glob(os.path.join(
    'outputs', 'Grouped_*', 'multiseed_summary_RAS_grouped_multiseed.xlsx')))
if not grouped_files:
    print('\nNo grouped summary found; run rerun_grouped.py first.')
    raise SystemExit

grp = load_summary(grouped_files[-1])
print(f'\ncomparing against {grouped_files[-1]}')

cmp = pd.DataFrame({
    'seq_mean': ung['AUC_Score_mean'], 'seq_sd': ung['AUC_Score_std'],
    'pat_mean': grp['AUC_Score_mean'], 'pat_sd': grp['AUC_Score_std'],
})
cmp['delta_mean'] = (cmp['pat_mean'] - cmp['seq_mean']).round(3)
cmp['sd_ratio'] = (cmp['pat_sd'] / cmp['seq_sd']).round(2)
cmp = cmp.round(3).sort_values('seq_mean', ascending=False)
cmp.to_csv(os.path.join(OUTDIR, 'SA5_like_for_like_comparison.csv'))

print('\n=== multi-seed AUC, five seeds, same machine ===')
print('  seq_* = partition drawn over sequences (original)')
print('  pat_* = partition drawn over patients (SA5)\n')
print(cmp.to_string())
print(f'\nwritten to {OUTDIR}')

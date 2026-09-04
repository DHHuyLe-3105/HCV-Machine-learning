"""
Re-run the RAS machine-learning pipeline with PATIENT-GROUPED splitting.

Why: the 162 sequences come from 146 patients. Sixteen patients contributed two
samples each (isolates named TF##A/TF##B and HCV_##A/HCV_##B, drawn 169-282 days
apart). Under the original random split, the two samples of one patient could
land on opposite sides of the train/test boundary, so a model could be tested on
a virus it had already seen. This script keeps every sample from a patient on
the same side and reports what that costs.

Nothing in the original notebook is edited. The notebook's own function
definitions are executed as written, and only the splitting function they call
is replaced.

    python rerun_grouped.py

Requires the repository layout (data/data_RAS.xlsx, model_utils.py,
Predictive_model_ras.ipynb) plus genbank_metadata.csv in the same folder.
"""
import json
import os
import re
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split as _sk_split

NOTEBOOK = 'Predictive_model_ras.ipynb'
RAS_PATH = 'data/data_RAS.xlsx'
META = 'genbank_metadata.csv'
OUTDIR = os.path.join('outputs', 'Grouped_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
SEEDS = (42, 7, 13, 21, 31)

os.makedirs(OUTDIR, exist_ok=True)

# ----------------------------------------------------------------- 1. patients
meta = pd.read_csv(META)
meta['accession'] = meta['accession'].astype(str).str.strip()


def patient_of(accession, isolate):
    iso = str(isolate)
    m = re.match(r'^(TF\d+)[AB]?$', iso)
    if m:
        return 'CAM_' + m.group(1)
    m = re.match(r'^(HCV_\d+)[AB]_S\d+$', iso)
    if m:
        return 'CAM_' + m.group(1)
    return accession


meta['patient'] = [patient_of(a, i) for a, i in zip(meta['accession'], meta['isolate'])]
acc2pat = dict(zip(meta['accession'], meta['patient']))

n_dup = (meta['patient'].value_counts() > 1).sum()
print(f'{len(meta)} sequences from {meta["patient"].nunique()} patients '
      f'({n_dup} patients contributed more than one sequence)')

# ------------------------------------------------- 2. run the notebook's code
with open(NOTEBOOK, encoding='utf-8') as fh:
    nb = json.load(fh)
cells = [''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code']

# cells that only define functions; the driver cells are deliberately skipped
DEF_CELLS = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11]
ns = {'__name__': '__main__'}
for i in DEF_CELLS:
    exec(compile(cells[i], f'<cell {i}>', 'exec'), ns)
print(f'loaded {len(DEF_CELLS)} definition cells from {NOTEBOOK}')

build_feature_matrix = ns['build_feature_matrix']
train_and_evaluate_models = ns['train_and_evaluate_models']
bootstrap_sensitivity_analysis = ns['bootstrap_sensitivity_analysis']
multi_seed_sensitivity = ns['multi_seed_sensitivity']

X, y, ras_df, ras_features, feature_names = build_feature_matrix(RAS_PATH, ras_sheet='Sheet1')
accessions = ras_df['Accession number'].astype(str).str.strip().tolist()
groups = np.array([acc2pat.get(a, a) for a in accessions])
missing = [a for a in accessions if a not in acc2pat]
if missing:
    print(f'WARNING: {len(missing)} accessions absent from {META}; each is treated '
          f'as its own patient: {missing[:5]}')
print(f'feature matrix {X.shape}, {len(set(groups))} patient groups')

# ------------------------------------------- 3. group-aware replacement split
def grouped_train_test_split(*arrays, test_size=None, random_state=None,
                             stratify=None, **kwargs):
    """Drop-in replacement for train_test_split that never separates two samples
    of the same patient. The notebook always calls it as
    (X, y, indices, test_size=..., random_state=..., stratify=...), so the row
    indices are read from the third array and mapped to patients."""
    if len(arrays) != 3:
        return _sk_split(*arrays, test_size=test_size,
                         random_state=random_state, stratify=stratify, **kwargs)

    Xa, ya, idx = arrays
    idx = np.asarray(idx)
    g = groups[idx]

    # one row per patient, labelled by that patient's outcome
    gdf = pd.DataFrame({'g': g, 'y': np.asarray(ya)}).drop_duplicates('g')
    strat = gdf['y'].values if stratify is not None else None
    g_tr, g_te = _sk_split(gdf['g'].values, test_size=test_size,
                           random_state=random_state, stratify=strat)
    in_test = np.isin(g, g_te)

    return (Xa[~in_test], Xa[in_test],
            np.asarray(ya)[~in_test], np.asarray(ya)[in_test],
            idx[~in_test], idx[in_test])


ns['train_test_split'] = grouped_train_test_split
print('splitting function replaced with the patient-grouped version')

# --------------------------------------------- 3b. keep FDA fittable
# QuadraticDiscriminantAnalysis needs more samples per class than features. The
# grouped split can leave a training class smaller than the surviving feature
# count, and sklearn then refuses to fit, which aborts the whole grid search.
# This subclass tries the ordinary fit first, so results are unchanged wherever
# the original would have worked, and only falls back to fitting the same model
# in a reduced principal-component space when the plain fit is impossible.
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as _QDA
from sklearn.decomposition import PCA as _PCA


class SafeQDA(_QDA):
    """QDA that degrades gracefully instead of aborting the run."""

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
print('FDA made rank-safe (falls back to a PCA space only if the plain fit fails)')

import sklearn, sys
try:
    import xgboost, imblearn
    print(f'versions: python {sys.version.split()[0]}  sklearn {sklearn.__version__}  '
          f'xgboost {xgboost.__version__}  imblearn {imblearn.__version__}')
except ImportError:
    print(f'versions: python {sys.version.split()[0]}  sklearn {sklearn.__version__}')

# ------------------------------------------------------------------ 4. run it
print('\n=== primary run, seed 42, patient-grouped 60/20/20 ===')
results = train_and_evaluate_models(
    X, y, seq_df=ras_df, feature_names=feature_names,
    label='RAS_grouped', save_path=OUTDIR, random_state=42)

print('\n=== bootstrap sensitivity ===')
bootstrap_sensitivity_analysis(results, n_bootstraps=1000,
                               label='RAS_grouped', save_path=OUTDIR)

print('\n=== multi-seed sensitivity ===')
multi_seed_sensitivity(X, y, ras_df, feature_names, seeds=SEEDS,
                       label='RAS_grouped_multiseed', save_path=OUTDIR)

# ------------------------------------------------------- 5. leakage diagnostic
idx_tr, idx_te = results['idx_train'], results['idx_test']
shared = set(groups[idx_tr]) & set(groups[idx_te])
print(f'\npatients appearing on BOTH sides of the split: {len(shared)} (should be 0)')

# how bad was it before? repeat the original split for the same seed
_, _, _, _, i_tr0, i_te0 = _sk_split(X, y, np.arange(len(y)), test_size=0.2,
                                     random_state=42, stratify=y)
shared0 = set(groups[i_tr0]) & set(groups[i_te0])
print(f'under the ORIGINAL random split, seed 42: {len(shared0)} patients were '
      f'split across train and test')
if shared0:
    print('  ', sorted(shared0))

print(f'\nall outputs written to {OUTDIR}')

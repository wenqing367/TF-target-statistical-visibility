from __future__ import annotations
import re
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from _paths import PROJECT_ROOT as ROOT, project_relative

RAW = ROOT / 'data' / 'zju_cellatlas' / 'raw'
DGE = RAW / 'Human_dge.h5ad'
OUT_DATA = ROOT / 'data' / 'processed' / 'zju_cellatlas_adult_tissue_subsets'
OUT_TABLES = ROOT / 'results' / 'tables_adult_tissue'
OUT_INV = ROOT / 'results' / 'extension_forward_visibility' / 'data_inventory'
MIN_CELLS = 1000

ADULT_TISSUES = {
    'AdultLung': ['AdultLung'],
    'AdultKidney': ['AdultKidney'],
    'AdultLiver': ['AdultLiver'],
    'AdultColon': ['AdultTransverseColon', 'AdultSigmoidColon', 'AdultAscendingColon'],
}

def slug(x: str) -> str:
    return re.sub(r'[^A-Za-z0-9]+', '_', str(x)).strip('_').lower()

def mean_detect(mat):
    if sparse.issparse(mat):
        mean = np.asarray(mat.mean(axis=0)).ravel()
        det = np.asarray((mat > 0).sum(axis=0)).ravel()
    else:
        arr = np.asarray(mat)
        mean = arr.mean(axis=0)
        det = (arr > 0).sum(axis=0)
    return mean, det / mat.shape[0], det

def make_summary(sub, path):
    counts = sub.layers['counts'] if 'counts' in sub.layers else sub.X
    mean, rate, det = mean_detect(counts)
    pd.DataFrame({
        'gene': sub.var_names.astype(str).str.upper(),
        'mean_counts': mean,
        'detection_rate': rate,
        'n_detected_cells': det.astype(int),
        'n_cells': sub.n_obs,
    }).drop_duplicates('gene').to_csv(path, index=False)

def to_csr(x):
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)

def main():
    if not DGE.exists():
        raise SystemExit(f'Missing HCL h5ad: {DGE}')
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    OUT_TABLES.mkdir(parents=True, exist_ok=True)
    OUT_INV.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(DGE, backed='r')
    if 'tissue' not in adata.obs.columns:
        raise SystemExit('HCL h5ad does not contain obs["tissue"]')
    tissue = adata.obs['tissue'].astype(str)
    rows=[]
    print('[raw_shape]', adata.shape)
    print('[adult_tissue_counts]')
    for context, labels in ADULT_TISSUES.items():
        mask = tissue.isin(labels)
        cells = list(adata.obs_names[mask])
        n = len(cells)
        status = 'prepared' if n >= MIN_CELLS else 'skipped_low_cells'
        prefix = f'zju_adult_tissue_{slug(context)}'
        out_h5ad = OUT_DATA / f'{prefix}_normalized.h5ad'
        out_summary = OUT_TABLES / f'{prefix}_gene_summary.csv'
        rows.append({
            'source':'ZJU_HCL',
            'analysis_level':'adult_tissue',
            'context':context,
            'source_tissue_labels': ';'.join(labels),
            'n_cells':n,
            'status':status,
            'h5ad':project_relative(out_h5ad),
            'gene_summary':project_relative(out_summary),
        })
        print(context, n, labels, status, flush=True)
        if n < MIN_CELLS:
            continue
        sub = adata[cells, :].to_memory()
        sub.var_names = sub.var_names.astype(str).str.upper()
        sub.var_names_make_unique()
        sub.X = to_csr(sub.X)
        sub.layers['counts'] = sub.X.copy()
        sub.write_h5ad(out_h5ad, compression='gzip')
        make_summary(sub, out_summary)
        print('[write]', out_h5ad, flush=True)
        print('[write]', out_summary, flush=True)
    inv = OUT_INV / 'zju_cellatlas_adult_tissue_contexts.csv'
    pd.DataFrame(rows).to_csv(inv, index=False)
    print('[write]', inv)
    adata.file.close()

if __name__ == '__main__':
    main()

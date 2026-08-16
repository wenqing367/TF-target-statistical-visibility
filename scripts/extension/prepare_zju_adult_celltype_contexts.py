from __future__ import annotations
import re
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from _paths import PROJECT_ROOT as ROOT, project_relative

RAW = ROOT / 'data' / 'zju_cellatlas' / 'raw'
DGE = RAW / 'Human_dge.h5ad'
INFO = RAW / 'Human_cell_info.xlsx'
OUT_DATA = ROOT / 'data' / 'processed' / 'zju_cellatlas_adult_celltype_subsets'
OUT_TABLES = ROOT / 'results' / 'tables_adult_celltype'
OUT_INV = ROOT / 'results' / 'extension_forward_visibility' / 'data_inventory'
MIN_CELLS = 1000

ADULT_CELLTYPES = [
    'Epithelial cell',
    'Endothelial cell',
    'Fibroblast',
    'Smooth muscle cell',
    'Enterocyte',
    'AT2 cell',
    'Loop of Henle',
]

def slug(x: str) -> str:
    return re.sub(r'[^A-Za-z0-9]+', '_', str(x)).strip('_').lower()

def normalize_cell_id(x: str) -> str:
    return re.sub(r'-\d+$', '', str(x))

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
    if not DGE.exists() or not INFO.exists():
        raise SystemExit(f'Missing HCL files under {RAW}')
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    OUT_TABLES.mkdir(parents=True, exist_ok=True)
    OUT_INV.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(DGE, backed='r')
    obs_names = pd.Index(adata.obs_names.astype(str))
    obs_base_to_name = pd.Series(obs_names.values, index=[normalize_cell_id(x) for x in obs_names])
    obs_base_to_name = obs_base_to_name[~obs_base_to_name.index.duplicated(keep='first')]

    info = pd.read_excel(INFO)
    info['cellnames'] = info['cellnames'].astype(str)
    info['_cell_base'] = info['cellnames'].map(normalize_cell_id)
    info['_obs_name'] = info['_cell_base'].map(obs_base_to_name)
    info = info.dropna(subset=['_obs_name']).copy()
    if info.empty:
        raise SystemExit('No overlap between cell info and h5ad after cell id normalization')
    if 'celltype' not in info.columns or 'stage' not in info.columns:
        raise SystemExit(f'Missing celltype/stage in HCL info: {list(info.columns)}')

    print('[raw_shape]', adata.shape)
    rows=[]
    for ct in ADULT_CELLTYPES:
        labels = info[(info['stage'].astype(str).eq('Adult')) & (info['celltype'].astype(str).str.lower().eq(ct.lower()))].copy()
        cells = list(dict.fromkeys(labels['_obs_name'].tolist()))
        n = len(cells)
        tissues = labels['sample'].astype(str).value_counts()
        status = 'prepared' if n >= MIN_CELLS else 'skipped_low_cells'
        prefix = f'zju_adult_celltype_{slug(ct)}'
        out_h5ad = OUT_DATA / f'{prefix}_normalized.h5ad'
        out_summary = OUT_TABLES / f'{prefix}_gene_summary.csv'
        rows.append({
            'source':'ZJU_HCL',
            'analysis_level':'adult_celltype',
            'context':ct,
            'n_cells':n,
            'n_source_tissues':len(tissues),
            'top_source_tissues': ';'.join([f'{k}:{int(v)}' for k,v in tissues.head(10).items()]),
            'status':status,
            'h5ad':project_relative(out_h5ad),
            'gene_summary':project_relative(out_summary),
        })
        print(f'[subset] {ct}: {n} adult cells, {len(tissues)} tissues, status={status}', flush=True)
        print('[top_tissues]', '; '.join([f'{k}:{int(v)}' for k,v in tissues.head(8).items()]), flush=True)
        if n < MIN_CELLS:
            continue
        sub = adata[cells, :].to_memory()
        sub.obs['adult_celltype_context'] = ct
        sub.var_names = sub.var_names.astype(str).str.upper()
        sub.var_names_make_unique()
        sub.X = to_csr(sub.X)
        sub.layers['counts'] = sub.X.copy()
        sub.write_h5ad(out_h5ad, compression='gzip')
        make_summary(sub, out_summary)
        print('[write]', out_h5ad, flush=True)
        print('[write]', out_summary, flush=True)
    inv = OUT_INV / 'zju_cellatlas_adult_celltype_contexts.csv'
    pd.DataFrame(rows).to_csv(inv, index=False)
    print('[write]', inv)
    adata.file.close()

if __name__ == '__main__':
    main()

from __future__ import annotations
import argparse
import importlib.util
import re
from pathlib import Path
import pandas as pd
from _paths import MATCHING_ROOT, PROJECT_ROOT as ROOT, resolve_project_path


def load_runner(script: Path):
    if not script.exists():
        raise SystemExit(f'Missing upstream matching script: {script}')
    spec = importlib.util.spec_from_file_location('refined_runner', script)
    runner = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(runner)
    return runner

def slug(x: str) -> str:
    return re.sub(r'[^A-Za-z0-9]+', '_', str(x)).strip('_').lower()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--context', required=True)
    ap.add_argument('--matching-root', default=str(MATCHING_ROOT))
    ap.add_argument('--out-root', default=str(ROOT / 'results' / 'zju_cellatlas_adult_tissue_refined_matching_by_context'))
    args = ap.parse_args()
    matching_root = Path(args.matching_root).expanduser().resolve()
    runner = load_runner(matching_root / 'scripts' / '43_known_edge_visibility_refined_matching.py')
    edge_sets = {
        'trrust': matching_root / 'data' / 'ground_truth' / 'trrust_edges_detection_ge_05pct.csv',
        'dorothea_ab': matching_root / 'data' / 'ground_truth' / 'dorothea_ab_edges_detection_ge_05pct.csv',
        'trrust_dorothea_intersection': matching_root / 'data' / 'ground_truth' / 'trrust_dorothea_ab_intersection_detection_ge_05pct.csv',
        'trrust_dorothea_union': matching_root / 'data' / 'ground_truth' / 'trrust_dorothea_ab_union_detection_ge_05pct.csv',
    }
    inv = pd.read_csv(ROOT / 'results' / 'extension_forward_visibility' / 'data_inventory' / 'zju_cellatlas_adult_tissue_contexts.csv')
    row = inv[(inv['status'].eq('prepared')) & (inv['analysis_level'].eq('adult_tissue')) & (inv['context'].eq(args.context))]
    if row.empty:
        raise SystemExit(f'No prepared adult tissue context found: {args.context}')
    r = row.iloc[0]
    config = {
        'dataset': 'ZJU_HCL_adult_tissue',
        'condition': str(r['context']),
        'h5ad': resolve_project_path(r['h5ad']),
        'gene_summary': resolve_project_path(r['gene_summary']),
        'edge_sets': edge_sets,
    }
    out_dir = Path(args.out_root) / f"adult_tissue__{slug(args.context)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print('[context]', config['dataset'], config['condition'])
    print('[h5ad]', config['h5ad'])
    print('[out]', out_dir)
    runner.helpers.DATASETS = [config]
    runner.run_analysis(out_dir, use_existing_components=False)

if __name__ == '__main__':
    main()

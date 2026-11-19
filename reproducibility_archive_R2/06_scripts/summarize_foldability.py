#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""summarize_foldability.py —— 本地运行，汇总四候选单体预测置信度。
自动识别 result_nomsa/（single_sequence 无外网版）或 result/。
输出 foldability_summary.csv / .md。
字段：mean/median pLDDT、frac>50/70、pTM、max PAE、独立模型间 CA RMSD。
"""
import os
import csv
import json
import glob
import math
import statistics

BASE = os.path.dirname(os.path.abspath(__file__))
_nomsa = os.path.join(BASE, 'result_nomsa')
RESULT = _nomsa if os.path.isdir(_nomsa) and any(
    os.path.isdir(os.path.join(_nomsa, c)) for c in ('66', '477', '499', '399')
) else os.path.join(BASE, 'result')
print('RESULT dir:', RESULT)
CAND = ['66', '477', '499', '399']


def plddt_from_pdb(path):
    """从 PDB 的 B-factor 列取 CA 的 pLDDT。"""
    vals = []
    for line in open(path, encoding='utf-8', errors='ignore'):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                pass
    return vals


def metrics_from_json(path):
    """从 scores json 取 pTM 与 max_pae。"""
    try:
        d = json.load(open(path, encoding='utf-8'))
    except Exception:
        return None, None
    ptm = d.get('ptm') if isinstance(d.get('ptm'), (int, float)) else None
    mpae = d.get('max_pae') if isinstance(d.get('max_pae'), (int, float)) else None
    if mpae is None and isinstance(d.get('pae'), list):
        try:
            mpae = max(max(r) for r in d['pae'])
        except Exception:
            pass
    return ptm, mpae


def ca_coords(path):
    d = {}
    for line in open(path, encoding='utf-8', errors='ignore'):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            try:
                r = int(line[22:26])
                d[r] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                pass
    return d


def rmsd(a, b):
    common = sorted(set(a) & set(b))
    if not common:
        return None
    ss = 0.0
    for k in common:
        for i in range(3):
            ss += (a[k][i] - b[k][i]) ** 2
    return math.sqrt(ss / len(common))


rows = []
for cid in CAND:
    d = os.path.join(RESULT, cid)
    pdbs = sorted(glob.glob(os.path.join(d, '**', '*model*.pdb'), recursive=True)) or \
           sorted(glob.glob(os.path.join(d, '**', '*.pdb'), recursive=True))
    jsons = sorted(glob.glob(os.path.join(d, '**', '*scores*.json'), recursive=True))

    pref = [p for p in pdbs if 'rank_001' in os.path.basename(p)]
    main_pdb = pref[0] if pref else (pdbs[0] if pdbs else None)
    pref_j = [j for j in jsons if 'rank_001' in os.path.basename(j)]
    ptm, mpae = metrics_from_json(pref_j[0] if pref_j else (jsons[0] if jsons else None))

    plddt = plddt_from_pdb(main_pdb) if main_pdb else []
    # 独立模型间 RMSD：取排名前二的不同 model
    inter_rmsd = None
    if len(pdbs) >= 2:
        inter_rmsd = rmsd(ca_coords(pdbs[0]), ca_coords(pdbs[1]))

    row = {
        'candidate': 'Ligand %s' % cid,
        'n_models': len(pdbs),
        'mean_pLDDT': round(statistics.mean(plddt), 2) if plddt else '',
        'median_pLDDT': round(statistics.median(plddt), 2) if plddt else '',
        'frac_pLDDT_gt50': round(sum(1 for v in plddt if v > 50) / len(plddt), 3) if plddt else '',
        'frac_pLDDT_gt70': round(sum(1 for v in plddt if v > 70) / len(plddt), 3) if plddt else '',
        'pTM': round(ptm, 3) if ptm is not None else '',
        'max_PAE_A': round(mpae, 2) if mpae is not None else '',
        'intermodel_CA_RMSD_A': round(inter_rmsd, 2) if inter_rmsd is not None else '',
    }
    rows.append(row)
    print(row)

keys = list(rows[0].keys())
with open(os.path.join(BASE, 'foldability_summary.csv'), 'w', newline='', encoding='utf-8-sig') as fh:
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    for r in rows:
        w.writerow(r)

md = ['# 问题3：四候选单体折叠置信度（ColabFold single_sequence 无外网版）', '',
      '| 候选 | 模型数 | mean pLDDT | median pLDDT | frac>50 | frac>70 | pTM | max PAE (Å) | 模型间CA RMSD (Å) |',
      '|---|---|---|---|---|---|---|---|---|']
for r in rows:
    md.append('| %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
        r['candidate'], r['n_models'], r['mean_pLDDT'], r['median_pLDDT'],
        r['frac_pLDDT_gt50'], r['frac_pLDDT_gt70'], r['pTM'],
        r['max_PAE_A'], r['intermodel_CA_RMSD_A']))
md += ['', '注：pLDDT 取 rank_001 模型 PDB 的 B-factor 列；pTM 与 max PAE 取对应 scores json；'
      '模型间 RMSD 为前两个独立模型的 CA 对齐偏差。']
open(os.path.join(BASE, 'foldability_summary.md'), 'w', encoding='utf-8').write('\n'.join(md))
print('\n输出: foldability_summary.csv / foldability_summary.md')

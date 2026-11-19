#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_membrane.py  ——  CSBJ R2 膜可及性重分析（路线A）
=========================================================================
回应审稿意见 #1：原稿只统计了【受体界面残基】相对膜平面的位置，没有量化
【配体本身】有多少原子插进膜、占据多大体积。本脚本从受体跨膜螺旋自动拟合
膜平面，对未定向的复合物做膜取向校正，再逐候选量化配体侧膜插入。

四个候选：66 / 477 / 499 / 399。链约定：A=配体 mini-protein，B=受体 Claudin-18。

输出（写入脚本所在目录）：
  aligned/<id>_aligned.pdb   受体TM轴旋正到z、膜心置z=0的复合物
  per_residue/<id>_per_residue.csv   逐残基入膜明细
  pymol/<id>_membrane_view.pml       PyMOL 全复合物膜取向图脚本
  summary.csv / summary.json         候选级入膜汇总
  membrane_plane_fit.txt             每个候选的膜平面拟合参数（方法学记录）

膜平面拟合方法：
  1) Kyte-Doolittle 疏水性(窗口19)在受体链上识别跨膜螺旋段；
  2) 取每段跨膜螺旋 CA 坐标，用最小二乘拟合螺旋轴向量；
  3) 统一轴向后平均得到膜法向 n；
  4) 膜心 = 各跨膜螺旋中点沿 n 的中位数；
  5) 膜半厚默认 16.2 Å（与原稿 PPM 值一致）；
  6) 把整个复合物旋转/平移，使 n→+z、膜心→z=0，再量化。
"""

import argparse
import csv
import json
import math
import os
from collections import defaultdict

KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5,
      'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
      'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
      'Y': -1.3, 'V': 4.2}
AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
       'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}
VDW = {'H': 1.10, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80,
       'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98, 'SE': 1.90}


def element_of(name, field):
    if field and field.strip():
        return field.strip().upper()
    letters = ''.join(c for c in name.upper() if c.isalpha())
    if len(letters) >= 2 and letters[:2] in ('CL', 'BR', 'CA', 'MG', 'ZN', 'FE', 'SE'):
        return letters[:2]
    return letters[:1] if letters else 'C'


def parse_pdb(path):
    """返回原子 dict 列表，保留原始行以便写回坐标。"""
    atoms = []
    for line in open(path, encoding='utf-8', errors='ignore'):
        if not (line.startswith('ATOM  ') or line.startswith('HETATM')):
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        atoms.append({'line': line.rstrip('\n'), 'name': line[12:16].strip(),
                      'resname': line[17:20].strip(), 'chain': line[21:22].strip(),
                      'resseq': line[22:27].strip(), 'x': x, 'y': y, 'z': z,
                      'element': element_of(line[12:16], line[76:78]),
                      'is_ca': line[12:16].strip() == 'CA'})
    return atoms


def detect_tm_segments(atoms, rec_chains, win=19, thr=1.6, minlen=8):
    """KD 疏水性识别受体跨膜螺旋段，返回 [(start,end),...]。"""
    seq = {}
    for a in atoms:
        if a['chain'] in rec_chains and a['is_ca']:
            try:
                seq[int(a['resseq'])] = AA3.get(a['resname'], 'A')
            except ValueError:
                pass
    nums = sorted(seq)
    if not nums:
        raise ValueError('受体链未找到 CA 原子')
    s = ''.join(seq[n] for n in nums)
    scores = []
    for i in range(len(s) - win + 1):
        scores.append((sum(KD.get(c, 0) for c in s[i:i + win]) / win,
                       nums[i], nums[i + win - 1]))
    segs = []
    for sc, a, b in [x for x in scores if x[0] > thr]:
        if segs and a <= segs[-1][1] + 1:
            segs[-1] = (segs[-1][0], b)
        else:
            segs.append((a, b))
    return [(a, b) for a, b in segs if b - a + 1 >= minlen]


def unit(v):
    n = math.sqrt(sum(t * t for t in v))
    return tuple(t / n for t in v) if n else v


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def axis_from_points(pts):
    """最小二乘直线方向（协方差矩阵最大特征向量，用幂迭代近似）。"""
    n = len(pts)
    cx = sum(p[0] for p in pts) / n
    cy = sum(p[1] for p in pts) / n
    cz = sum(p[2] for p in pts) / n
    xx = sum((p[0] - cx) ** 2 for p in pts) / n
    yy = sum((p[1] - cy) ** 2 for p in pts) / n
    zz = sum((p[2] - cz) ** 2 for p in pts) / n
    xy = sum((p[0] - cx) * (p[1] - cy) for p in pts) / n
    xz = sum((p[0] - cx) * (p[2] - cz) for p in pts) / n
    yz = sum((p[1] - cy) * (p[2] - cz) for p in pts) / n
    v = (1.0, 0.3, 0.2)
    for _ in range(100):
        nv = (xx * v[0] + xy * v[1] + xz * v[2],
              xy * v[0] + yy * v[1] + yz * v[2],
              xz * v[0] + yz * v[1] + zz * v[2])
        nv = unit(nv)
        if dot(nv, v) > 0.999999:
            v = nv
            break
        v = nv
    return v, (cx, cy, cz)


def median(vals):
    s = sorted(vals)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0


def rotation_align(n):
    """构造把单位向量 n 旋到 +z 的旋转矩阵 Rodrigues。返回 3x3。"""
    z = (0.0, 0.0, 1.0)
    d = dot(n, z)
    if d > 0.999999:
        return ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    if d < -0.999999:  # 翻转180°绕x轴
        return ((1, 0, 0), (0, -1, 0), (0, 0, -1))
    k = unit(cross(n, z))
    kx, ky, kz = k
    c, s, t = d, math.sqrt(1 - d * d), 1 - d
    R = (
        (c + t * kx * kx, t * kx * ky - s * kz, t * kx * kz + s * ky),
        (t * ky * kx + s * kz, c + t * ky * ky, t * ky * kz - s * kx),
        (t * kz * kx - s * ky, t * kz * ky + s * kx, c + t * kz * kz),
    )
    return R


def apply_R(p, R):
    return (dot(R[0], p), dot(R[1], p), dot(R[2], p))


def fit_membrane(atoms, rec_chains):
    """拟合膜法向、膜心、跨膜段。"""
    segs = detect_tm_segments(atoms, rec_chains)
    if len(segs) < 2:
        raise ValueError('跨膜螺旋段不足2条，无法拟合膜平面: %s' % segs)
    ca = defaultdict(dict)
    for a in atoms:
        if a['chain'] in rec_chains and a['is_ca']:
            try:
                ca[a['chain']][int(a['resseq'])] = (a['x'], a['y'], a['z'])
            except ValueError:
                pass
    normal0 = None
    axes, midpoints, all_tm = [], [], []
    for a, b in segs:
        pts = []
        for ch in ca:
            for r in range(a, b + 1):
                if r in ca[ch]:
                    pts.append(ca[ch][r])
        if len(pts) < 4:
            continue
        v, mid = axis_from_points(pts)
        if normal0 is not None and dot(v, normal0) < 0:
            v = tuple(-t for t in v)
        else:
            normal0 = v
        axes.append(v)
        midpoints.append(mid)
        all_tm.extend(pts)
    n = unit(tuple(sum(v[i] for v in axes) for i in range(3)))
    center = tuple(median([dot(p, n) for p in all_tm]) * n[i] for i in range(3))
    return n, center, segs


def transform(atoms, center, normal):
    """把复合物旋正：n→+z，膜心→z=0。返回新坐标列表。(不写PDB)"""
    R = rotation_align(normal)
    coords = []
    for a in atoms:
        p = (a['x'] - center[0], a['y'] - center[1], a['z'] - center[2])
        coords.append(apply_R(p, R))
    return coords


def write_aligned_pdb(atoms, coords, path):
    with open(path, 'w', encoding='utf-8') as fh:
        for a, (x, y, z) in zip(atoms, coords):
            line = a['line']
            newline = line[:30] + '%8.3f%8.3f%8.3f' % (x, y, z) + line[54:]
            fh.write(newline + '\n')
        fh.write('END\n')


def quantify(atoms, coords, lig_chains, half):
    """在已旋正坐标上量化配体入膜。"""
    lig = [(a, c) for a, c in zip(atoms, coords)
           if a['chain'] in lig_chains and a['element'] != 'H']
    if not lig:
        raise ValueError('配体链无重原子')
    embedded = [(a, c) for a, c in lig if abs(c[2]) <= half]
    n_heavy = len(lig)
    n_emb = len(embedded)
    vol = sum((4.0 / 3.0) * math.pi * VDW.get(a['element'], 1.70) ** 3
              for a, c in embedded)
    max_depth = max(abs(c[2]) for a, c in lig)
    centroid_z = sum(c[2] for a, c in lig) / n_heavy

    res = defaultdict(list)
    for a, c in lig:
        res[(a['chain'], a['resname'], a['resseq'])].append(abs(c[2]) <= half)
    per_res = []
    for key, flags in sorted(res.items(), key=lambda kv: int(kv[0][2]) if kv[0][2].lstrip('-').isdigit() else 0):
        ne = sum(flags)
        per_res.append({'chain': key[0], 'resname': key[1], 'resseq': key[2],
                        'heavy': len(flags), 'embedded': ne,
                        'frac': round(ne / len(flags), 3),
                        'buried': ne / len(flags) >= 0.5})
    return {
        'n_ligand_heavy': n_heavy, 'n_embedded_heavy': n_emb,
        'embedded_fraction': round(n_emb / n_heavy, 4),
        'embedded_volume_A3': round(vol, 1),
        'max_abs_z_A': round(max_depth, 2), 'centroid_z_A': round(centroid_z, 2),
        'n_buried_residues': sum(r['buried'] for r in per_res),
        'n_ligand_residues': len(per_res), 'per_residue': per_res,
    }


def write_pymol(path, aligned_pdb, lig_chains, half, atoms_aligned):
    L = []
    L.append('# PyMOL 膜取向图：受体TM轴已旋正到z轴，膜心z=0')
    L.append('reinitialize')
    L.append('load "%s", complex' % os.path.abspath(aligned_pdb).replace('\\', '/'))
    L.append('hide everything')
    L.append('show cartoon, chain B')
    L.append('show sticks, chain %s' % '+'.join(sorted(lig_chains)))
    L.append('color grey, chain B')
    L.append('color lightblue, chain %s' % '+'.join(sorted(lig_chains)))
    emb = []
    for a in atoms_aligned:
        if a['chain'] in lig_chains and a['element'] != 'H':
            if abs(a['z']) <= half:
                emb.append('(chain %s and resi %s and name "%s")' %
                           (a['chain'], a['resseq'].strip(), a['name']))
    if emb:
        L.append('select ligand_in_membrane, ' + ' or '.join(emb))
        L.append('color red, ligand_in_membrane')
    span = 40.0
    L.append('pseudoatom upper, pos=[0,0,%g]' % half)
    L.append('pseudoatom lower, pos=[0,0,%g]' % (-half))
    L.append('show spheres, upper or lower')
    L.append('set sphere_scale, %.2f, upper or lower' % (span / 1.7))
    L.append('color yellow, upper or lower')
    L.append('set transparency, 0.85, upper or lower')
    L.append('bg white')
    L.append('set orthoscopic, on')
    L.append('viewport 1100,800')
    L.append("label upper, 'upper leaflet  z=+%.1f A'" % half)
    L.append("label lower, 'lower leaflet  z=%.1f A'" % (-half))
    open(path, 'w', encoding='utf-8').write('\n'.join(L) + '\n')


def main():
    base = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description='CSBJ R2 膜可及性重分析')
    ap.add_argument('--rec', default='B', help='受体链，逗号分隔 (默认 B)')
    ap.add_argument('--lig', default='A', help='配体链，逗号分隔 (默认 A)')
    ap.add_argument('--half', type=float, default=16.2, help='膜半厚 Å (默认16.2)')
    ap.add_argument('pdbs', nargs='*', help='复合物PDB；留空则用内置4个候选路径')
    args = ap.parse_args()

    rec_chains = set(args.rec.split(','))
    lig_chains = set(args.lig.split(','))
    pdbs = args.pdbs or [
        'D:/calculation/批量操作代码/筛选/11.17/66/co-folding/66.pdb',
        'D:/calculation/批量操作代码/筛选/11.17/399/399.pdb',
        'D:/calculation/批量操作代码/筛选/11.17/477/477.pdb',
        'D:/calculation/批量操作代码/筛选/11.17/499/499.pdb',
    ]

    for d in ('aligned', 'per_residue', 'pymol'):
        os.makedirs(os.path.join(base, d), exist_ok=True)

    rows, fitlog = [], []
    for pdb in pdbs:
        cid = os.path.splitext(os.path.basename(pdb))[0]
        atoms = parse_pdb(pdb)
        n, center, segs = fit_membrane(atoms, rec_chains)
        coords = transform(atoms, center, n)
        aligned = os.path.join(base, 'aligned', cid + '_aligned.pdb')
        write_aligned_pdb(atoms, coords, aligned)
        aligned_atoms = parse_pdb(aligned)

        q = quantify(aligned_atoms,
                     [(a['x'], a['y'], a['z']) for a in aligned_atoms],
                     lig_chains, args.half)

        pr_path = os.path.join(base, 'per_residue', cid + '_per_residue.csv')
        with open(pr_path, 'w', newline='', encoding='utf-8-sig') as fh:
            w = csv.writer(fh)
            w.writerow(['chain', 'resname', 'resseq', 'heavy_atoms',
                        'embedded_atoms', 'embedded_frac', 'buried_>=50%'])
            for r in q['per_residue']:
                w.writerow([r['chain'], r['resname'], r['resseq'], r['heavy'],
                            r['embedded'], r['frac'], r['buried']])

        write_pymol(os.path.join(base, 'pymol', cid + '_membrane_view.pml'),
                    aligned, lig_chains, args.half, aligned_atoms)

        print('\n=== %s ===' % cid)
        print('  跨膜段(KD):', segs)
        print('  膜法向: (%.3f, %.3f, %.3f)' % n)
        print('  配体重原子 %d | 入膜 %d (%.1f%%)' %
              (q['n_ligand_heavy'], q['n_embedded_heavy'],
               q['embedded_fraction'] * 100))
        print('  入膜体积(加和近似) %.1f Å^3 | 埋膜残基 %d/%d' %
              (q['embedded_volume_A3'], q['n_buried_residues'], q['n_ligand_residues']))

        row = {'candidate': cid,
               'TM_segments': ';'.join('%d-%d' % s for s in segs),
               'normal_x': round(n[0], 4), 'normal_y': round(n[1], 4),
               'normal_z': round(n[2], 4),
               'center_x': round(center[0], 3), 'center_y': round(center[1], 3),
               'center_z': round(center[2], 3), 'half_thickness_A': args.half}
        for k in ('n_ligand_heavy', 'n_embedded_heavy', 'embedded_fraction',
                  'embedded_volume_A3', 'max_abs_z_A', 'centroid_z_A',
                  'n_buried_residues', 'n_ligand_residues'):
            row[k] = q[k]
        rows.append(row)
        fitlog.append('%s: TM=%s  n=(%.3f,%.3f,%.3f)  center=(%.2f,%.2f,%.2f)' %
                      (cid, segs, n[0], n[1], n[2], center[0], center[1], center[2]))

    keys = list(rows[0].keys()) if rows else []
    with open(os.path.join(base, 'summary.csv'), 'w', newline='', encoding='utf-8-sig') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    json.dump(rows, open(os.path.join(base, 'summary.json'), 'w', encoding='utf-8'),
              ensure_ascii=False, indent=2)
    open(os.path.join(base, 'membrane_plane_fit.txt'), 'w', encoding='utf-8').write(
        '膜平面拟合记录（Kyte-Doolittle TM识别 + CA最小二乘轴向 + Rodrigues旋正）\n' +
        '膜半厚 half-thickness = %.1f Å（与原稿PPM一致）\n\n' % args.half +
        '\n'.join(fitlog) + '\n')
    print('\n输出目录:', base)


if __name__ == '__main__':
    main()

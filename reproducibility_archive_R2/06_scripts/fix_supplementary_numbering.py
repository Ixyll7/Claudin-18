# -*- coding: utf-8 -*-
"""
解决导师批注34：Supplementary Figure S1 呢？要按照顺序来排。
方案：Figure 与 Table 分开独立连续编号；正文引用同步对齐；补回丢失的两处引用。
所有改动文字加黄色高亮。生成新文件，不覆盖原件。
"""
import zipfile, xml.etree.ElementTree as ET, os, shutil

W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
Wp = W + 'p'
Wr = W + 'r'
Wt = W + 't'
WrPr = W + 'rPr'
Whl = W + 'highlight'

SI_SRC = 'E:/蛋白对接/写文章用/我的/CSBJ/Supplementary Information_R3.docx'
SI_DST = 'E:/蛋白对接/写文章用/我的/CSBJ/Supplementary Information_R4.docx'
MS_SRC = 'E:/蛋白对接/写文章用/我的/CSBJ/返修/R2/membrane_reanalysis_R2/manuscript(r5)_foldability.docx'
MS_DST = 'E:/蛋白对接/写文章用/我的/CSBJ/返修/R2/membrane_reanalysis_R2/manuscript(r6)_suppOrder.docx'

# 编号映射：旧 -> 新（Figure/Table 分开独立连续编号）
# Figure: S2 -> S1
# Table: S1->S1, S3->S2, S4->S3, S5->S4, S6->S5
FIG_MAP = {'2': '1'}
TBL_MAP = {'1': '1', '3': '2', '4': '3', '5': '4', '6': '5'}

def set_highlight(run):
    rPr = run.find(WrPr)
    if rPr is None:
        rPr = ET.Element(WrPr)
        run.insert(0, rPr)
    hl = rPr.find(Whl)
    if hl is None:
        hl = ET.SubElement(rPr, Whl)
    hl.set(W + 'val', 'yellow')

def para_text(p):
    return ''.join((t.text or '') for t in p.iter(Wt))

def rebuild_para(p, text, keep_bold=False):
    """把段落正文 run 重建为单个高亮 run，保留 pPr，保留首个 run 的加粗属性"""
    pPr = p.find(W + 'pPr')
    # 记录原首 run 是否加粗
    bold = False
    for r in p.findall(Wr):
        rPr = r.find(WrPr)
        if rPr is not None and rPr.find(W + 'b') is not None:
            bold = True
            break
    # 删除所有 run
    for r in list(p.findall(Wr)):
        p.remove(r)
    newr = ET.Element(Wr)
    rPr = ET.SubElement(newr, WrPr)
    if bold or keep_bold:
        ET.SubElement(rPr, W + 'b')
    hl = ET.SubElement(rPr, Whl)
    hl.set(W + 'val', 'yellow')
    t = ET.SubElement(newr, Wt)
    t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    t.text = text
    p.append(newr)

def replace_in_runs(p, replacer, mark=True):
    """在每个 w:t 上做字符串替换；mark=True 则给改动的 run 加高亮"""
    changed = False
    for r in p.iter(Wr):
        for t in r.findall(Wt):
            if t.text:
                new = replacer(t.text)
                if new != t.text:
                    t.text = new
                    changed = True
                    if mark:
                        set_highlight(r)
    return changed

log = []

# ========== 1. 处理 SI ==========
shutil.copy(SI_SRC, SI_DST)
z = zipfile.ZipFile(SI_SRC)
root = ET.fromstring(z.read('word/document.xml'))
body = root.find(W + 'body')
children = list(body)

import re
# 匹配 "Supplementary Figure S2" 或 "Supplementary Table S3" 等（编号数字单独 run，故需在整段层面处理）
pat_full = re.compile(r'(Supplementary\s+(Figure|Table)\s+S)(\d+)', re.I)

for idx in [0, 6, 7, 10, 13, 16]:
    p = children[idx]
    full = para_text(p)
    def repl(m):
        prefix, kind, num = m.group(1), m.group(2), m.group(3)
        newnum = FIG_MAP.get(num, num) if kind.lower() == 'figure' else TBL_MAP.get(num, num)
        return prefix + newnum
    newtext = pat_full.sub(repl, full)
    if newtext != full:
        rebuild_para(p, newtext)
        log.append('SI body#%d: %s -> %s' % (idx, full[:45], newtext[:45]))
    else:
        log.append('SI body#%d: NO CHANGE' % idx)

# 写回 SI
zin = zipfile.ZipFile(SI_SRC, 'r')
tmp = SI_DST + '.tmp'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/document.xml':
            data = ET.tostring(root, encoding='utf-8', xml_declaration=True)
        zout.writestr(item, data)
zin.close()
os.replace(tmp, SI_DST)
log.append('SI written: %s' % SI_DST)

# ========== 2. 处理 manuscript ==========
shutil.copy(MS_SRC, MS_DST)
z2 = zipfile.ZipFile(MS_SRC)
root2 = ET.fromstring(z2.read('word/document.xml'))
body2 = root2.find(W + 'body')
children2 = list(body2)

# 正文引用替换（注意顺序：先长后短，用一次性正则避免链式替换）
# 正文出现: S5(membership), Figure S2, Table S3, S5, S5, S4
# 映射: Figure S2->Figure S1; Table S3->Table S2; Table S4->Table S3; Table S5->Table S4
ms_pat = re.compile(r'(Supplementary\s+(?:Figure|Table)\s+S)(\d+)', re.I)
# 需要根据 Figure/Table 类型区分的映射
def ms_repl(m):
    prefix, num = m.group(1), m.group(2)
    if 'figure' in prefix.lower():
        newnum = FIG_MAP.get(num, num)
    else:
        newnum = TBL_MAP.get(num, num)
    return prefix + newnum

for i, p in enumerate(children2):
    if p.tag != Wp:
        continue
    full = para_text(p)
    newtext = ms_pat.sub(ms_repl, full)
    if newtext != full:
        rebuild_para(p, newtext)
        log.append('MS body#%d: numbering updated' % i)

# 补回丢失的引用：
# (a) 膜分析段 p#（含 "fraction and volume of ligand atoms occupying the bilayer"）-> 补 Table S1
# (b) 折叠段（含 "five independent single-sequence"）-> 补 Table S5（新编号，原S6）
adds = [
    ('fraction and volume of ligand atoms occupying the bilayer',
     ' The residue-level interface map of the Claudin-18 / Ligand 66 complex is provided in Supplementary Table S1.'),
    ('five independent single-sequence',
     ' The full set of monomer folding-confidence metrics is provided in Supplementary Table S5.'),
]
for keyword, addition in adds:
    done = False
    for i, p in enumerate(children2):
        if p.tag != Wp:
            continue
        full = para_text(p)
        if keyword in full and addition.strip() not in full:
            newtext = full.rstrip() + addition
            rebuild_para(p, newtext)
            log.append('MS body#%d: added citation -> %s' % (i, addition[:55]))
            done = True
            break
    if not done:
        log.append('WARN: add-citation keyword not found or already present: %s' % keyword[:40])

zin2 = zipfile.ZipFile(MS_SRC, 'r')
tmp2 = MS_DST + '.tmp'
with zipfile.ZipFile(tmp2, 'w', zipfile.ZIP_DEFLATED) as zout:
    for item in zin2.infolist():
        data = zin2.read(item.filename)
        if item.filename == 'word/document.xml':
            data = ET.tostring(root2, encoding='utf-8', xml_declaration=True)
        zout.writestr(item, data)
zin2.close()
os.replace(tmp2, MS_DST)
log.append('MS written: %s' % MS_DST)

print('\n'.join(log))

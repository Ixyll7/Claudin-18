# 问题 7 操作指南：把数据仓库做到审稿人满意

审稿人 #7 的核心要求：出版前提供一个**固定的、公开可访问的归档**，包含完整打分数据、
界面分析、精确配置文件、预测输出和复现候选筛选所需的脚本。

我实地核查了你的仓库 `https://github.com/Ixyll7/Claudin-18` 的当前状态：

| 审稿人说的问题 | 实际情况 | 是否已解决 |
|---|---|---|
| "no releases available"（没建 release） | Releases 页确实显示 "There aren't any releases here" | ❌ 需要你建 |
| README 提到的打分表看不到 | `results/γ.xlsx` 其实**已在 10-23 上传到 `results/` 目录**（不在根目录，审稿人没展开） | ⚠️ 需让路径更醒目 |
| 补充表没在 101 页 PDF 里 | 补充表是**单独提交的 SI 文档**，本体没问题 | ✅ 提交时确保 SI 一起交 |
| 要 fixed archive | GitHub main 分支是可变的，不算 fixed | ❌ 建议加 Zenodo DOI |

---

## 你要做的 4 步（按顺序）

### 第 1 步：把本归档目录上传到 GitHub

把本文件夹 `reproducibility_archive_R2/` 里的 6 个子目录（01–06）整个上传到仓库，
建议在仓库根目录新建一个文件夹，例如 `R2_reproducibility/`，保持子目录结构。

用网页上传：GitHub 仓库页 → `Add file` → `Upload files` → 拖拽 01–06 六个文件夹 →
`Commit changes`。

### 第 2 步：创建 v1.0 Release（这是审稿人明确点名的）

1. 打开 https://github.com/Ixyll7/Claudin-18/releases
2. 点 **Draft a new release**
3. Tag 填 `v1.0`（这就是审稿人找的东西）
4. Title 填 `v1.0 — R2 reproducibility archive`
5. Description 可以写：
   `Fixed archive of scoring data, interface analyses, prediction outputs and scripts supporting the R2 revision. Includes the complete candidate-selection dataset and post-selection analyses.`
6. 勾选 **Target: main**
7. 点 **Publish release**

做完这一步，审稿人再打开 releases 页就能看到 v1.0，"no releases" 的问题就没了。

### 第 3 步：更新 README，让打分表路径醒目

README 里已经写了 `results/γ.xlsx`，但建议再加一段显眼的说明，让审稿人一眼看到。
在 README 开头或 "Candidate scoring results" 段加：

```
## Reproducibility archive (v1.0)

A fixed, versioned archive is available as the `v1.0` release and in directory
`R2_reproducibility/`. It contains:
- Complete per-candidate scoring tables (MEGADOCK, docking γ, co-folding γ) for all 500 candidates
- Candidate-level membership table for the consensus selection
- Library properties, membrane-accessibility, foldability and FoldX mutation-scan data
- Scripts to reproduce each analysis

The complete scoring spreadsheet is `results/γ.xlsx` (in the `results/` subdirectory).
```

### 第 4 步：生成永久 DOI（强烈建议，让 "fixed" 名副其实）

GitHub release 虽然带 tag，但严格说还不是"fixed citable archive"。最稳妥的做法是用 Zenodo
给仓库绑一个 DOI：

1. 打开 https://zenodo.org ，用 GitHub 账号登录
2. Settings → GitHub → 勾选 `Ixyll7/Claudin-18` 仓库，点 **Connect**
3. 回到仓库，去 releases 页对 v1.0 点编辑并保存一次（会触发 Zenodo 归档）；
   或者直接在 Zenodo 上点 **New deposit** → 授权选该 release → 填标题/作者 → **Publish**
4. 发布后会得到一个 DOI，形如 `10.5281/zenodo.xxxxxxx`
5. 把这个 DOI 发我，我写进稿件的 Data availability 段落

---

## 我这边已经帮你做好的

- 建立了归档目录 `reproducibility_archive_R2/`，含 6 个子目录、14 个文件：
  - `01_screening_scores/` 三张全库打分表 + 候选 membership 表
  - `02_library_properties/` pTM + 序列性质表
  - `03_membrane_analysis/` 四候选入膜汇总
  - `04_foldability/` 折叠置信度汇总
  - `05_mutation_scan/` Cys109 + Gln65 PositionScan 原始输出
  - `06_scripts/` 4 个可复现脚本
- 写了归档 README，逐文件说明内容与软件版本

## 稿件这边我会做的（等你确认 DOI 后）

Data availability 段落会从现在的一句光秃秃链接，改成规范表述，大致是：

> The complete scoring data, interface analyses, prediction outputs, exact software configurations
> and analysis scripts supporting this study are available in a fixed, versioned archive at
> [GitHub v1.0 release URL] and [Zenodo DOI]. Supplementary Tables S1–S6 and Supplementary Figure S1
> are provided in the Supplementary Information.

你把 GitHub release 建好、（如果做了）Zenodo DOI 拿到后发我，我就把这段写进稿件并高亮。

---

## 优先级说明

如果时间紧，**第 1、2 步是必须的**（审稿人点名了 release 和数据可得性），第 4 步 Zenodo DOI 是
锦上添花但很值得做（让 "fixed" 站得住）。第 3 步一分钟就能改，顺手做了。

注意：真正"复现候选筛选"还需要 500 个 mini-protein 的三维结构——这些已经在仓库的
`chai--500_mini_proteins-0--250.zip` 和 `chai--500_mini_proteins-251--499.zip` 里（审稿人 ignored
的那批），所以数据链是完整的，不用重新生成。

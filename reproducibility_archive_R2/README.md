# Reproducibility archive — R2 revision

Fixed public archive of the data and scripts required to reproduce the candidate-selection and
post-selection analyses reported in the manuscript *AI-Driven Virtual Screening of Mini-Protein
Candidates against Claudin-18: Membrane Accessibility as a Critical Filter in Computational
Prioritization* (CSBJ-D-26-00765R2).

This archive is intended to be deposited as a **versioned, citable release** (see UPLOAD_GUIDE.md).
Files are grouped by the analysis they support. All tabular data are plain CSV/TXT and can be read
with standard tools.

## 01_screening_scores/ — consensus candidate selection

| File | Content |
|---|---|
| `megadock_scores-1.csv` | MEGADOCK docking score for all 500 candidates (higher = better complementarity). Field `number` = candidate ID; `megadock` = score. |
| `docking_rosetta_gamma.csv` | Docking-derived Rosetta γ (= dG_separated/dSASA×100) for all 500 complexes. Field `folder_03` = candidate ID; `dG_separated/dSASAx100` = γ (more negative = better). |
| `cofolding_rosetta_gamma.csv` | Co-folding-derived Rosetta γ for all 500 complexes (same fields as above). |
| `screening_membership.csv` | Candidate-level membership table: each of the four final candidates with its three scores, ranks, top-set membership and three-way intersection. |

**Selection algorithm reproduced:** three parallel filters on the 500-candidate library —
MEGADOCK score top 20% (n = 100), docking-derived γ top 10% (n = 50), co-folding-derived γ
top 20% (n = 100); the four final candidates are the three-way intersection. Entries that did
not yield a valid score in a given filter are excluded from that filter's ranked list.
`build_selection_table.py` regenerates the membership table from the three full-library files.

## 02_library_properties/ — random library characterization

| File | Content |
|---|---|
| `ptm_scores.csv` | Chai-1 pTM for all 500 random mini-proteins (IDs `03_0`–`03_499`). |
| `merged_properties.csv` | Per-sequence length, GRAVY, Shannon entropy, normalized entropy and pTM. |

## 03_membrane_analysis/ — membrane accessibility

| File | Content |
|---|---|
| `membrane_occupancy_summary.csv` | For each of the four candidates: fraction and volume of ligand atoms occupying the bilayer, and number of buried residues, after orienting each complex to a fitted membrane plane (half-thickness 16.2 Å). |

`analyze_membrane.py` performs the orientation and occupancy calculation from the complex PDBs.

## 04_foldability/ — monomer folding confidence

| File | Content |
|---|---|
| `foldability_summary.csv` | Per candidate: mean/median pLDDT, fraction of residues with pLDDT > 50 / > 70, pTM, maximum PAE and inter-model Cα RMSD, from five independent single-sequence ColabFold predictions (alphafold2_ptm, no external MSA). |

`summarize_foldability.py` extracts these values from raw ColabFold outputs (PDB B-factor column
and scores JSON).

## 05_mutation_scan/ — FoldX PositionScan

| File | Content |
|---|---|
| `Cys109_positionscan_raw.txt` | Raw FoldX PositionScan output at Ligand 66 residue 109 (all substitutions; ΔE = E_mut − E_WT). |
| `Gln65_positionscan_raw.txt` | Raw FoldX PositionScan output at residue 65, included for completeness. |

Note: ΔE is a single-structure mutation-scan energy, **not** a binding free energy (ΔΔG_bind);
no bound/unbound thermodynamic cycle was computed.

## 06_scripts/ — analysis scripts

Python scripts used to generate the tables above from the raw data files. See the header comment
in each script for its inputs and outputs. `fix_supplementary_numbering.py` is an editorial helper
that renumbers the supplementary items and is not part of the scientific analysis.

## Software versions

Chai-1 v0.6.1 (structure prediction); MEGADOCK v4.1.1 (rigid docking); ColabFold 1.5.x with
AlphaFold2-multimer default parameters (co-folding and single-sequence prediction); Rosetta v3.14
(fastrelax and InterfaceAnalyzer); FoldX 5.0 (mutation scanning).

## Corresponding supplementary items

Supplementary Tables S1–S6 and Supplementary Figure S1 are provided in the separately submitted
Supplementary Information document. Complete residue-level interface mapping and energy results
are in Tables S1 and S3, respectively.

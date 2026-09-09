
# The series of genetic codes that preceded our own

Data and code for Masel, Douglas, Wehbi & McShea, *The series of genetic codes
that preceded our own*.

The order of amino acid recruitment in Table 1 is a revision of Wehbi et al.
(2024) *Proc. Natl. Acad. Sci. USA* **121**:e2410311121. The **base method** —
Pfam age classification, ancestral sequence reconstruction, and clan-level amino
acid usage — is unchanged from that paper, and its code lives at
[sawsanwehbi/Pfam-age-classification](https://github.com/sawsanwehbi/Pfam-age-classification);
Pfam sequences, alignments, trees and protein-ID mappings are on figshare
(<<DOI>>). **The modification made here** is to the post-LUCA control
denominator only, and is implemented in `<<table1_pipeline.R>>` — see
"Modification to Wehbi et al. (2024)" below.

## Contents

| Path | What it is |
|---|---|
| `<<table1_pipeline.R>>` | Applies the revised post-LUCA denominator and produces the Table 1 values |
| `<<clan_classification.csv>>` | Every clan, its LUCA/post-LUCA assignment, and which exclusion criterion (if any) removed it from the denominator |
| `AA_properties_Fig1_SupFig1.csv` | Per-amino-acid data for Figures 1, 6 and Supplementary Figure 1. See the data dictionary below |
| `Scatter_Plots.R` | Figures 1a, 1b, 6 and Supplementary Figure 1 |
| `phylogenetics/` | BEAST 2 analyses of the Class I and Class II aaRS catalytic domains (Figure 5, Supplementary Files 3–4). See `phylogenetics/README.md` |
| `sessionInfo.txt` | R and package versions used |

## Modification to Wehbi et al. (2024)

Wehbi et al. (2024) compared reconstructed amino acid usage in LUCA-era protein
clans against clans presumed to be ancient but post-LUCA, which served as a
methodological control. That paper was careful to establish the antiquity of the
LUCA-era clans, but did not take equivalent precautions to *exclude* antiquity
among the control clans.

Here the numerator is unchanged and the denominator is made more stringent. Two
classes of clan are excluded from the post-LUCA control:

1. clans that had already diversified into multiple Pfams by the time of the
   Last Archaeal Common Ancestor (LACA) or the Last Bacterial Common Ancestor
   (LBCA); and
2. clans with at least one Pfam, annotated as a LACA or LBCA candidate, that had
   duplicated prior to LACA or LBCA respectively.

This removes <<n>> of <<m>> post-LUCA control clans, leaving <<k>>. LUCA usage is
then the reconstructed amino acid frequency at confidently inferred sites for
clans inferred to be present in one copy in LUCA, divided by similarly
reconstructed frequencies in the retained control clans. Nothing else about the
method differs from Wehbi et al. (2024).

## Reproducing the figures

```r
install.packages(c("ggplot2", "broom", "common"))
source("Scatter_Plots.R")   # reads AA_properties_Fig1_SupFig1.csv from this directory
```

## Data dictionary: `AA_properties_Fig1_SupFig1.csv`

| Column | Description |
|---|---|
| `AA` | One-letter amino acid code |
| `Nb_nonH_sidechain` | Number of non-hydrogen atoms in the side chain (G = 0, A = 1, … W = 10) |
| `Trifonov_order` | Rank in the consensus order of recruitment of Trifonov (2000), with ties given the same rank: G/A, V/D, P, S, E/L, T, R, N, K, Q, I, C, H, F, M, Y, W |
| `LUCA_usage` | Revised LUCA clan usage, as printed in Table 1: reconstructed amino acid frequency in LUCA clans divided by that in the revised post-LUCA control clans |
| `LUCA_usageSE` | Standard error of `LUCA_usage`, as printed in Table 1 |
| `aaRS_Class` | Structural class of the aminoacyl-tRNA synthetase that activates the amino acid: `1`, `2`, or `1_2` for lysine, which has both LysRS-I and LysRS-II |
| `<<Recruitment_order>>` | Rank order of recruitment implied by `LUCA_usage`, i.e. the Table 1 order, with ties as represented in Figure 2 (1 = G/V/T … 13 = W). **Not** the order published in Wehbi et al. (2024) |
| `mean_protozymedistance` | Number of insertion modules acquired in the aaRS catalytic domain since the protozyme (Douglas et al. 2024, plus one additional insertion in LysRS-I as discussed in Figure 4), averaged across the forms of that aaRS. Supplementary File 5 |
| `min_protozymedistance` | As above, but the minimum across forms rather than the mean |


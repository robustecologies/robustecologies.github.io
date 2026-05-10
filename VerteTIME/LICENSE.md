# VerteTIME licensing

VerteTIME ships two licence categories:

  

## Code (R source, namespace, build infrastructure)

The R package code (everything in `R/`, `src/`, `tests/`, `NAMESPACE`,
the `License` field of `DESCRIPTION`) is licensed under the GNU General
Public License v3.0 or later. Full text:
<https://www.gnu.org/licenses/gpl-3.0.html>.

  

## Data, vignettes, manuscript, figures

The harmonised data compilation (`data/*.rda` shipped with the package,
every file under `web-export/vertetime-v1.0/`, and every figure produced
by
[`vt_publish()`](https://robustecologies.github.io/VerteTIME/reference/vt_publish.md)
or by knitting the manuscript), the package vignettes and the long-form
journal manuscript are licensed under Creative Commons Attribution 4.0
International (CC-BY 4.0). Full text:
<https://creativecommons.org/licenses/by/4.0/legalcode>.

  

## Provenance and attribution

VerteTIME’s community time series are independent re-curations from
peer-reviewed primary literature and grey literature. The compilation is
the intellectual property of Sergio Picó and Pablo Almaraz. Secondary
biodiversity time-series compilations (LPI, BioTIME, GPDD, RivFishTIME,
PREDICTS, BBS, PECBMS, CBC, eBird Trends, TEAM/Wildlife Insights and
others) are referenced only as positional comparators in the manuscript
timeline figure; they are not the source of any VerteTIME series.
Per-dataset primary references and DOIs are recorded in
`data_provenance.csv` and reproduced in the appendix of the manuscript.

When citing VerteTIME, please cite both the package (R
`citation("VerteTIME")`) and the long-form manuscript once it is
published. Per-dataset reuse must additionally cite the primary
reference recorded in `data_provenance.csv` for the dataset in question.

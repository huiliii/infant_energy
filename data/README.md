# Data

Datasets are not redistributed here. Download or generate each item
below, then point `code/config.m` at the resulting folders.

Expected layout (defaults assumed by `config.m`):

```
data/
├── dhcp_fetal.mat                       # fetal cohort (see below)
├── dhcp_neo.mat                         # neonatal cohort
├── dhcp_fet_ce_alff.mat                 # produced by s01_compute_ce_alff.m
├── dhcp_neo_ce_alff.mat                 # produced by s01_compute_ce_alff.m
├── group_corr.mat                       # produced by s01_compute_ce_alff.m
├── residuals.mat                        # ALFF residuals + net_index
├── infant_aal_spin_permutations_50k.mat # spin-test permutations
├── brainspan/
│   ├── expression_matrix.csv
│   ├── columns_metadata.csv
│   ├── rows_metadata.csv
│   ├── suppletable13.csv                # Kang et al. 2011 gene-set table
│   └── trsp_brain_regions.mat           # aal_rns1, aal_rns2
├── ahba/
│   └── AHBA_adult_expression.csv
└── cmrglc/
    └── adult_CMRglc.mat                 # adult_CMRglc_aal (90 x 1)
```

## Sources

- **dHCP (fetal + neonatal MRI).** Apply for access at
  https://www.developingconnectome.org/data-release/. The `mat` field in
  `dhcp_fetal.mat` / `dhcp_neo.mat` is a `90 x 90 x N` structural
  connectivity stack and `tcs` is a `T x 90 x N` ROI time series, with
  T = 350 for fetuses (TR = 2.2 s) and T = 2300 for neonates
  (TR = 0.392 s). Cohort-specific fields used by the scripts include
  `age` / `pma` (in weeks), `index_term` (logical, 1 for term),
  `tsnr_z_fet` / `tsnr_z_neo` (fMRI motion proxy), and `pct_slice` /
  `shard_outlier` (dMRI motion).

- **BrainSpan.** Download the "RNA-seq Gencode v10 summarized to genes"
  matrix from https://www.brainspan.org/. `suppletable13.csv` reproduces
  Supplementary Table S13 from Kang et al. 2011 (Nature). The
  AAL <-> BrainSpan region mapping in `trsp_brain_regions.mat`
  follows the convention used in our prior work.

- **AHBA.** Allen Human Brain Atlas microarray data, parcellated to the
  infant AAL90 atlas. Any standard pipeline (e.g. abagen) can produce
  `AHBA_adult_expression.csv`; the first column is the region ID and
  subsequent columns are genes.

- **Adult CMRglc.** Group-averaged FDG-PET map across 28 healthy adults
  from Shokri-Kojori et al. 2019, https://github.com/eshoko/COMET,
  registered to the AAL90 atlas.

- **Spin permutations.** Pre-computed 50 000 spherical-rotation
  permutations of the 90-node AAL parcellation, used for `spin_corr`.

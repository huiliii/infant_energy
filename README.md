# Energy characteristics of brain connectomes before and after birth

MATLAB code accompanying the manuscript on perinatal development of
network-control-theory **control energy** and its coupling to the
**amplitude of low-frequency fluctuations (ALFF)** in fetuses and
infants from the Developing Human Connectome Project (dHCP).

The pipeline:

1. Estimates per-subject regional control energy from individual
   structural connectomes (linear continuous-time optimal control).
2. Computes per-subject regional ALFF from preprocessed resting-state
   fMRI time series.
3. Tests how the spatial coupling between control energy and ALFF
   develops across the perinatal period.
4. Contextualises regional deviations in this coupling with
   developmental (BrainSpan) and adult (AHBA) gene-expression data and
   with adult cerebral glucose metabolism (FDG-PET; Shokri-Kojori et
   al. 2019).

## Repository layout

```
.
├── README.md
├── code/
│   ├── config.m                       # set local paths once, here
│   ├── s01_compute_ce_alff.m          # CE + ALFF per subject + group corr
│   ├── s02_ce_alff_coupling.m    # Fig. 2 analyses
│   ├── s03_cmrglc_alignment.m    # Fig. 5 analyses
│   ├── s04_transcriptomic_analysis.m  # Fig. 4 (BrainSpan + AHBA)
│   └── helpers/
│       ├── compute_subject_ce_alff.m  # CE + ALFF for one subject
│       └── compute_cohort_ce_alff.m   # loop wrapper over a cohort
└── data/
    └── README.md                      # data sources, expected layout
```

Scripts are numbered in the order they should be run; downstream
scripts depend on `.mat` files produced by upstream ones.

## Requirements

- MATLAB R2021a or later
- Statistics and Machine Learning Toolbox (`mafdr`, `partialcorr`,
  `ttest`, `ttest2`)
- Bioinformatics Toolbox (`mafdr` with `'BHFDR'`)
- A network control theory implementation providing `optim_fun`
  (see https://github.com/BassettLab/control_package)
- Helper functions distributed in `code/helpers/`:
  - `compute_subject_ce_alff.m`, `compute_cohort_ce_alff.m`
  - `cal_alff.m` — power-spectrum ALFF, 0.01-0.06 Hz (not included;
    standard implementation, e.g. from REST or in-house utilities)
  - `spin_corr.m` — spherical rotation test
    (see Alexander-Bloch et al. 2018, NeuroImage)
  - `scatter_corr_plot.m` — small plotting helper

## Quick start

1. Copy `code/config.m` and edit the paths at the top to point at your
   local copies of the data and the NCT toolbox.
2. Open MATLAB, `cd` into `code/`, and run:

   ```matlab
   run s01_compute_ce_alff.m         % builds dhcp_*_ce_alff.mat
   run s02_ce_alff_coupling.m   % Fig. 2 statistics
   run s03_cmrglc_alignment.m   % Fig. 5 statistics
   run s04_transcriptomic_analysis.m % Fig. 4 statistics
   ```

`s01_compute_ce_alff.m` has a `do_compute` flag; set it to `false` to
skip the optimisation loop and just load the cached `*_ce_alff.mat`
files.

## Data

Imaging and transcriptomic data are not redistributed in this
repository. See `data/README.md` for sources and the expected file
layout. Briefly, you will need:

- dHCP fetal and neonatal structural connectomes and ROI time series
  (request from https://www.developingconnectome.org/data-release/)
- BrainSpan gene-expression matrix and metadata
- An AAL-parcellated AHBA adult-expression CSV
- The adult CMRglc map from Shokri-Kojori et al. 2019, registered to
  the AAL90 atlas
- Spin-test permutations for the AAL90 atlas

## Pipeline overview

| Step | Script | Outputs |
| ---- | ------ | ------- |
| 1 | `s01_compute_ce_alff.m` | `dhcp_fet_ce_alff.mat`, `dhcp_neo_ce_alff.mat`, `group_corr.mat` |
| 2 | `s02_ce_alff_coupling.m` | whole-brain CE ~ age, cross-cohort CE/ALFF spatial similarity, group differences in r(CE, ALFF), r(CE, ALFF) ~ age (with and without motion/sex covariates) |
| 3 | `s03_cmrglc_alignment.m` | per-subject r(CE, adult CMRglc), r ~ age |
| 4 | `s04_transcriptomic_analysis.m` | BrainSpan gene-set tests, AHBA per-gene t-tests, ranked gene lists |

## Notes on the implementation

- The structural matrix is normalised as
  `A_norm = A / (max_eig + 1) - I` to keep the continuous-time system
  stable, following standard NCT conventions.
- ALFF is computed on the 0.01-0.06 Hz band; this is constrained by
  the fetal Nyquist (`1/(2 * TR) = 0.227 Hz`) but kept identical
  across cohorts for comparability.
- Subjects are kept in cohort-native order so that `index_term` and
  age vectors align with the connectivity stacks.
- `s04` exposes alternative definitions of the regions of interest
  (subcortical-network mask, 0→1 transition in residual sign,
  threshold on Δ(ALFF residual)) as commented blocks.

## Citation


## License

All rights reserved.

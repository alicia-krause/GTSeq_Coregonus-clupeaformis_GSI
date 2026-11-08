# GTSeq_Coregonus-clupeaformis_GSI

Scripts used to assess genetic stock identification and spatiotemporal stock mixing of lake whitefish (*Coregonus clupeaformis*) in Lake Michigan. This research was conducted at the University of Wisconsin-Stevens Point’s USGS Fisheries Co-op Molecular Conservation Genetics Lab and preformed by A.L. Krause, P.D. Barry, Y. Shi, W.A. Larson, D.A. Iserman, S.P. Hansen, T.O. Brenden, J.A. VanDeHey, and J.J. Homola.

Genomic data was generated with a new genotyping-in-thousands (GT-seq) panel. Genotyping and quality control was competed with [GTscore](https://github.com/gjmckinney/GTscore?tab=readme-ov-file) v1.3. Genetic stock identification was completed with [rubias](https://github.com/eriqande/rubias) v0.3.


## Contents in this Repository

| Script  | Purpose |
| ------------- | ------------- |
| 01_GTscore_genotyping.R  | PolyGen genotyping, genepop export, sample & locus summaries  |
| 02_Build_Mixture_Rubias.R  | merge runs, 30% missingness filters, metadata, mixture rubias files  |
| 03_Baseline_Development.R  | 465 loci → 295 loci, STRUCTURE-based migrant/admixture removal, PCA outliers, final baseline  |
| 04_Baseline_Evaluation.R  | self-assignment, 100%, and realistic simulations  |
| 05_Mixture_GSI.R  |  rubias `infer_mixture` (PB), 6 chains, Gelman-Rubin, CIs, P=0  |
| 06_SampleSize_Analysis.R  | Fig 6 subsampling and RMSE vs mixture size  |
| 07_Figures.R  | Figs 2, 4, 5, Suppl Fig 5  |

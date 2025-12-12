# sctransform revised by jliu678 to enable flexible analytic pearson residule


## Intro of Original Version
R package for normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression.

'The sctransform package was initially developed by [Hafemeister and Satija, Genome Biology 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1). Its version two improved significantly  [(Choudhary and Satija, Genome Biology, 2022)](https://doi.org/10.1186/s13059-021-02584-9).
Core functionality of this package has been integrated into [Seurat](https://satijalab.org/seurat/)'

## Quick start
I add a new function as below to allow the theta_given can be effetive together with `vst.flavor = "v2"` to use the stabler glmgampoi (implemented in sctransform v2) to run the analytic pearson residual. Check my blog to learn the logic flow of the original `sctransform()`
```
fit_glmGamPoi_offset_fixed_theta <- function(
    umi,
    model_str,
    data,
    theta_given = NULL,
    allow_inf_theta = FALSE
)
```
example: 
```
SCTransform(
  object,
  vst.flavor = "v2",
  theta_given = 100
)
```

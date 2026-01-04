# sctransform revised by jliu678 to enable flexible analytic pearson residule


## Intro of Original Version
R package for normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression.

'The sctransform package was initially developed by [Hafemeister and Satija, Genome Biology 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1). Its version two improved significantly  [(Choudhary and Satija, Genome Biology, 2022)](https://doi.org/10.1186/s13059-021-02584-9).
Core functionality of this package has been integrated into [Seurat](https://satijalab.org/seurat/)'

## Quick start
I add new functions including below to allow the theta_given to be effetive together with `vst.flavor = "v2"`, thus to use the stabler glmgampoi (implemented in sctransform v2) to run the analytic pearson residual. Check [my blog](https://jiyuanliu.netlify.app/project/improved_scrnaseq_normalization/) to learn the logic flow of the original `sctransform()` and the improved UMI normalization model. 

I also added another parameter `nb_mle_theta` which need be specified manualy (e.g. using fitted theta values obtained when theta is not fixed) when above theta_given is set. Because fixed theta cannot represent the theta_mle that is used in the condition (theta_mm/theta_mle < 1e-3) to determine whether a gene is poisson to exclude. Its default value is NULL


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
cur_theta <- 100
cur_as <- "RNA_scDesign"

n_gene <- nrow(seuratobj)
gene_names <- rownames(seuratobj)

theta_fitted <- seuratobj[[glue({cur_as}'_SCT')]]@SCTModel.list$model1@feature.attributes$theta
names(theta_fitted) <- rownames(seuratobj[[glue({cur_as}'_SCT')]]@SCTModel.list$model1@feature.attributes)
summary(theta_fitted)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.003994 0.033570 0.095337      Inf 0.288633      Inf 

# below try to avoid error by making sure theta_fitted has save length with theta
theta_fitted<-c(theta_fitted,
                setNames(rep(Inf,nrow(seuratobj[["RNA_scDesign"]])-length(theta_fitted)),
                         setdiff(rownames(seuratobj[["RNA_scDesign"]]),names(theta_fitted))))

theta_fitted <- data.frame('theta'=theta_fitted,row.names = names(theta_fitted))

SCTransform(
  seuratobj,
  verbose = T,
  variable.features.n = 3000, #our input has 1k genes, so make sure that the output includes them all
  vst.flavor = "v2",
  assay = cur_as,
  theta_given = setNames(rep(cur_theta,n_gene),nm=gene_names),
  nb_mle_theta = theta_fitted,
  new.assay.name = glue('{cur_as}_SCT_theta{as.integer(cur_theta)}')
)

```

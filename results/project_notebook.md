---
title: Single cell prediction
subtitle: Project notebook
---

# 1/May/2017 - 16/May/2017 

## 8/May/2017

### Description 

*Quan* provided a script via Slack to make predictions for the **Cardiodiff data**:

```
https://computationalgenomics.slack.com/files/quan/F59DM8QV7/lasso_extract_cvfit_parallel.r 
```

This script is designed to predict cell cluster-identity using gene expression data. This script is fixed to implement lasso only. It also uses the foreach package to parralelize a bootstrapping step.

### TODO

- Adapt script to run other models besides lasso
- Make a comparison of prediction accuracy across different regularization models such as lasso, elastic net (with different alpha parameters) and ridge

# 15/May/2017

## Description 

A comparison between lasso, ridge and elastic net was performed using `lasso_extract_cvfit_parallel.r` script. However, current implementation doest permit using the same training and test datasets if run for different models. 

Cardiodiff datasets used:

- **Cluster cell labels**: my.clusters_0.45_Day2.RDS
- **Expression data**: Exprs_DCVLnorm_unlog_minus1_pos_Day2.RDS
- **Differential expressed genes**: DEseq_Cluster1_vs_OtherClusters_Day2.txt_filtered_pAdjusted_sorted.txt

Data features:

- Day: 2

## Results

The following results corresponds to predicting cell cluster-identity (cluster 1 vs all other clusters). As seen in figure below, elastic net method with an alpha value of `0.1` tends to perform better than all remaining methods. Lasso accuracy is lower than all elastic net methods. Ridge is the worst method, probably because it uses all features to make the predictions (and not all differentially expressed genes may contribute to explain cluster identity). 


![](project_notebook_img/model_comparison_diff-rand-samples.png)


![](project_notebook_img/model_comparison_diff-rand-samples_density.png)


## Goals

- Adapt script to take same random samples for training and test data to make a valid comparison between regularization methods
- Modularize code

## Perspectives

- Possible R package implementation?


# 16/May/2017

## Description

Conceptual workflow was designed to restructure `lasso_extract_cvfit_parallel.r` in order to use different models to make predictions. See flowchart below.

## Results

![](project_notebook_img/single_cell_prediction_script_pipeline.png)

Coding process started.

# 17/May/2017

## Description

`lasso_extract_cvfit_parallel.r` was restructured into `prediction.Rmd`. Modularization of `Lit_New_Lasso` and accuracy extraction was done. 

`prediction.Rmd` was run to predict cells corresponding to cluster **1** in the day 2 for the **Cardiodiff** data. 100 interations were set. An error occured:

```
Error in { : task 100 failed - "NA/NaN argument"
Calls: %dopar% -> <Anonymous>
Execution halted
```

The 100th iteration failed. One possible explanation could be in the sampling process:

```R
  # Random sampling for training data --------------------------------------
  # Take random sample of cells corresponding to cluster of interest
  cluster.select.indx <- sample(cluster.select, 
                                size = round(length(cluster.select)/2), 
                                replace = FALSE)
  # Take random sample of cells corresponding to cluster to be compared
  cluster.compare.indx <- sample(cluster.compare, 
                                 size = round(length(cluster.compare)/2), 
                                 replace = FALSE)
  
  # Build predictor matrix -------------------------------------------------
  # Prepare predictor matrix containing both clutering classes
  predictor <- expression.data[features, c(cluster.select.indx, cluster.compare.indx)]
  
  # Generate categorical response ------------------------------------------
  # Set all values to cluster of interest
  response <- rep(clusterID.select, ncol(predictor))
  cluster.compare.names <- colnames(expression.data[,cluster.compare])
  sub.clustercompare.indx <- which(colnames(predictor) %in% cluster.compare.names)
```


# 18/May/2017

## Description

`prediction.Rmd` now includes three functions:

- `BuildTrainTest`: Creates Training and testing datasets
- `FitRegModel`: Runs regularization models using the output from `BuildTranTest`
- `ExtractResults`: Extracts results from *glmnet* object

It also uses the same training and testing datasets for all models.

## Results

`prediction.Rmd` was run to predict cells corresponding to cluster **1** in the day 2 for the **Cardiodiff** data. 50 bootstrapping interations were performed.

The following plots show that the **elastic net** algorithm performs better than lasso and ridge methods. The best model in terms of accuracy was using an alpha parameter of `0.1`.

![](project_notebook_img/model_comparison_18-05-2017.png)
![](project_notebook_img/model_comparison_density_18-05-2017.png)

@Joseph commented:

> great! However, we will need to have a justification for the alpha level - but the results shows that it’s better (than lasso) for all values! 


@Quan notice the following:
> it looks great, do you have the number of genes for each model?
> 
> usually biologists want to select a small set of genes that can predict well the cell groups

Considering the previous, the following plot was generated. 
![](project_notebook_img/number_genes_model_18-05-2017.png)


Elastic net method with and alpha equals 0.1 is the best. However, it uses **539 genes** in average while **lasso only uses 192 genes**.

Five-number + mean of number of genes included in each model:

|       |elastic.net.0.1  |elastic.net.0.5 |elastic.net.0.9 |lasso           | ridge        |
|:------|:---------------:|:--------------:|:--------------:|:--------------:|:------------:|
|Min    | 331.0           |102.0           |87.0            |82.0            |2817          |
|1st Qu.| 424.2           |183.2           |155.5           |116.0           |2817         |
|Median | 534.0           |214.5           |206.0           |174.0           |2817         |
|Mean   | 539.3           |239.4           |205.5           |192.3           |2817         |
|3rd Qu.| 646.0           |290.8           |253.8           |251.5           |2817         |
|Max.   | 767.0           |448.0           |386.0           |412.0           |2817         |

Considering @Quan's and @Joseph's comments on the alpha value, we should think about considering the trade off between accuracy and number of genes incorporated(the fewer genes are included for the regressions in each model, the lower the accuracy). The amount of genes in this case depends on the value of alpha.

## TODO

- Fix error associated to sampling `"NA/NaN argument"`
- Run regularization methods allowing interaction between variables
  + To reduce the number of combinations between coefficients when fitting the regression models, consider taking into account TF networks to decide which interaction between genes (features) are pertinent.

---
title: Single cell prediction
subtitle: Project notebook
author: Jose Alquicira Hernandez
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

**See commit** [ef27d15](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/commit/ef27d152ea7fb15c443d8231f7c4ade49ce89bd2)

## TODO

- Fix error associated to sampling `"NA/NaN argument"`
- Run regularization methods allowing interaction between variables
  + To reduce the number of combinations between coefficients when fitting the regression models, consider taking into account TF networks to decide which interaction between genes (features) are pertinent.

# 19/05/2017

# Description

`prediction.Rmd` was modified to find out training and testing dataset features that are problematic and causes the error `"NA/NaN argument"`. Random seeds were tracked to reproduce errors.

Desktop analysis (in laptop) of a small dataset with 120 genes and 5905 cells produced a error in the first iteration:

```R
##############################################################
Iteration: 1 
---------- 
Seed for cluster select: 30308
Seed for cluster compare: 29664

##############################################################
Iteration: 2 
---------- 
Seed for cluster select: 30288
Seed for cluster compare: 15632
---------- 
test    : 4347 20 
training: 2952 20 
response: 2952 
##############################################################
---------- 
test    : 4347 20 
training: 2952 20 
response: 2952 
##############################################################

##############################################################
Iteration: 3 
---------- 
Seed for cluster select: 4703
Seed for cluster compare: 11372
---------- 
test    : 4334 20 
training: 2952 20 
response: 2952 
##############################################################

##############################################################
Iteration: 4 
---------- 
Seed for cluster select: 20287
Seed for cluster compare: 15691
---------- 
test    : 4365 20 
training: 2952 20 
response: 2952 
##############################################################

##############################################################
Iteration: 5 
---------- 
Seed for cluster select: 16195
Seed for cluster compare: 6986
---------- 
test    : 4351 20 
training: 2952 20 
response: 2952 
##############################################################

##############################################################
Iteration: 6 
---------- 
Seed for cluster select: 16235
Seed for cluster compare: 2033
---------- 
test    : 4335 20 
training: 2952 20 
response: 2952 
##############################################################

##############################################################
Iteration: 7 
---------- 
Seed for cluster select: 12599
Seed for cluster compare: 10577
---------- 
test    : 4335 20 
training: 2952 20 
response: 2952 
##############################################################

##############################################################
Iteration: 8 
---------- 
Seed for cluster select: 26503
Seed for cluster compare: 30302
---------- 
test    : 4371 20 
training: 2952 20 
response: 2952 
##############################################################

##############################################################
Iteration: 9 
---------- 
Seed for cluster select: 10426
Seed for cluster compare: 14164
---------- 
test    : 4329 20 
training: 2952 20 
response: 2952 
##############################################################

##############################################################
Iteration: 10 
---------- 
Seed for cluster select: 7977
Seed for cluster compare: 16071
---------- 
test    : 4338 20 
training: 2952 20 
response: 2952 
##############################################################
Error in { : task 1 failed - "NA/NaN argument"
```

In the debugging process, the following error was found for iteration 1 using the same random seeds inside the `foreach` parallel loop:

```R
> predict.marker <- FitRegModel(test, training, response, family = "binomial", alpha = 0.1)
Error in 1:cvfit.dev.lambda.idx[1] : NA/NaN argument
```

Error was tracked manually using the random seeds form iteration 1. The error origin is the following:

```R
cvfit.dev.lambda.idx <- which(round(cvfit.dev$lambda,digit = 3) == round(cvfit$lambda.min,digits = 3))
```

Sometimes `cvfit.dev$lambda.min` is not found in `cvfit.dev$lambda`


```R
> round(cvfit.dev$lambda,digit = 3)
  [1] 1.817 1.655 1.508 1.374 1.252 1.141 1.040 0.947 0.863 0.786 0.716 0.653 0.595 0.542 0.494 0.450 0.410 0.374 0.340 0.310 0.283 0.258 0.235 0.214 0.195
 [26] 0.178 0.162 0.147 0.134 0.122 0.112 0.102 0.093 0.084 0.077 0.070 0.064 0.058 0.053 0.048 0.044 0.040 0.036 0.033 0.030 0.028 0.025 0.023 0.021 0.019
 [51] 0.017 0.016 0.014 0.013 0.012 0.011 0.010 0.009 0.008 0.008 0.007 0.006 0.006 0.005 0.005 0.004 0.004 0.004 0.003 0.003 0.003 0.002 0.002 0.002 0.002
 [76] 0.002 0.002 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
```

```R
round(cvfit$lambda.min,digits = 3)
[1] 0.111
```

**See commit** [c31d926](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/commit/c31d926e3a61ab85e82b7c95eefccc22baaf0b91)

From previous results, we can see that **0.112 ~ 0.111**. `0.112` is the closest value to `cvfit.dev$lambda`.


The following algorithm in [stat.ethz.ch](https://stat.ethz.ch/pipermail/r-help/2012-February/302827.html) was implemented to find the nearest lambda value to `cvfit.dev$lambda.min` in vector `cvfit.dev$lambda`:

```R
cvfit.dev.lambda.idx <- which.min(abs(cvfit.dev$lambda -cvfit$lambda.min))
cvfit.dev <- cvfit.dev[1:cvfit.dev.lambda.idx,]
```

See commit [da91950](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/commit/da91950303a123398375794aba44341739efbc95)

## Results

The following comparison between lasso, elastic net and ridge was performed for 100 bootstrap replicates using the new version of `prediction.Rmd`. No errors were generated due to selecting the minimum alpha value. We observe the same tendency: elastic net still performs better (particulary using an alpha value of `0.1`).


![](project_notebook_img/model_comparison_19-05-2017.png)
![](project_notebook_img/model_comparison_density_19-05-2017.png)

|         |elastic.net.0.1 |elastic.net.0.5 |elastic.net.0.9 |    lasso     |    ridge    |
|:--      |---------------:|---------------:|---------------:|-------------:|------------:|
|Min      |        301.0   |         84.0   |         73.0   |         64.0 |        2817 |
|1st Qu.  |        459.2   |        179.8   |        134.0   |        135.5 |        2817 |
|Median   |        543.5   |        227.5   |        171.5   |        184.0 |        2817 |
|Mean     |        543.5   |        247.3   |        190.0   |        186.4 |        2817 |
|3rd Qu.  |        623.8   |        317.8   |        240.5   |        230.2 |        2817 |
|Max.     |        775.0   |        491.0   |        379.0   |        393.0 |        2817 |

If the number of genes is considered to decide which model is better, a weighted accuracy (using the number of genes included in each model) may be useful measure. The following plot shows the weighted accuracy of each model (notice that ridge has a weighted accuracy of zero since it uses all genes for the prediction).

![](project_notebook_img/model_comparison_weighted_accuracy_19-05-2017.png)

See commit [4ae9b50](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/commit/4ae9b505a8d5eec34391d0c1cccdb9f2601a3529)

Results may be found on delta in the results directory

```
results/2017-05-19_CardioDiffAnalysis/prediction_cardiodiff_day2
```
## TODO

- Run 100 replicates using the IPSC single cell data


# 22/05/2017

## Description

hiPSC data was copied from @Quan's directory to `/shares/common/users/j.alquicira/SingleCell_Prediction/data/2017-05-22_HiPSC` on Delta.

```
cp /shares/common/groups/Group-Powell/shares-data/powell/quan/Expression_data_HiPSC_5day0Samples.RDS .
```

Differential gene expression was provided by @Quan via [Slack](https://files.slack.com/files-pri/T0F3BL6HX-F5E1FVDCL/download/significant_degenes_hipsc_day0.tar.gz)


## Results

Accuracy across models are similar. These results agree with the estimations reported in [previous analysis](http://biorxiv.org/content/early/2017/03/22/119255). Only 99 genes were used as features.

![](project_notebook_img/model_comparison_hipsc_1_vs_234_22-05-2017.png)

The number of genes included in all models (except ridge) is not variable across models as well.

![](project_notebook_img/number_genes_model_hipsc_1_vs_234_22-05-2017.png)

See commit [8fc3d60](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/commit/8fc3d60ac564a1933a92109c73e8834f5e80f35d)

> Note: style figures were modified. See commit [945839b](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/commit/945839b11d5f622b1e5e558dd6de7db618cccc78)

## TODO

- Run analysis comparing cluster 4 versus others
- Analyse pathways and identified co-regulation (in order to add variable interaction in models)


# 23/05/2017

## Results

Deviance explained was plotted for all models. Note that this results are similar to the ones [reported already](http://biorxiv.org/content/early/2017/03/22/119255). In the following figure the deviance explanation is shown. The cluster comparison was 1 vs 2,3 and 4.

![](project_notebook_img/model_hipsc_deviance_1_vs_234_23-05-2017.png)


# 24/05/2017

## Results

Model accuracy for lasso is better than the one [reported previously](http://biorxiv.org/content/early/2017/03/22/119255) (around 90% versus 99% - 100%). 

![](project_notebook_img/model_comparison_hipsc_4_vs_123_24-05-2017.png)

In the density plot below we can see a bimodal distribution of accuracy across all models.

![](project_notebook_img/model_comparison_hipsc_density_4_vs_123_24-05-2017.png)

Lasso and elastic net 0.9 chose a smaller number of genes for making the predictions.

![](project_notebook_img/number_genes_model_hipsc_4_vs_123_24-05-2017.png)


The performance is also better when comparing the deviance explained by the models.

![](project_notebook_img/model_hipsc_deviance_4_vs_123_24-05-2017.png)

See commit [945839b](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/commit/945839b11d5f622b1e5e558dd6de7db618cccc78)

## TODO

In a meeting with Joseph, the following points were discussed:

- Expand the methods to be evaluated: 
  + Multivariate regression
  + Stepwise regression
  + Random forest
  + SVM
  + Multinomial naive bayes
- Check these papers:
  + [Additive Genetic Variability and the Bayesian Alphabet](http://www.genetics.org/content/183/1/347)
  + [Extension of the bayesian alphabet for genomic selection](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-186)
- Keep using the HiPSC data to compare the models
- Look for R packages/code by *Ben Hayes* for Bayesian models

# 27/06/2017

Code has been restructured into and R package. See [scPrediction](https://github.com/IMB-Computational-Genomics-Lab/scPrediction)

The following methods have been incorporated for the prediction tasks:

- Generalized linear model
- Lasso *
- Elastic Net *
- Ridge *
- MARS (Multivariate adaptive regression splines
- Random Forests
- Naïve Bayes
- SVM (linear kernel)

\* *Already incorporated*

The R package [splatter](https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html) is being evaluated to make simulations and perform predictions using simulated gene expression datasets. 

We could consider the following simulations tuning specific parameters:

- **Number of genes**: Evaluate a range of genes (from 100 to 2000).
  + Number of cells: 100
  + Number of groups: 2 (50 in each group)
  + Number of bootstrap replicates (50 per dataset)

A script called `prediction_simulation.Rmd` has been [added](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/master/bin/prediction_simulation.Rmd).

# 29/06/2017

`prediction_simulation.Rmd` now permits run predictions from multiple simulated datasets.

`prediction_simulation.Rmd` was run under the version included in this commit [a561572](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/commit/a56157285acbfe78f3f93263bde0eb3169d2789a).

Results may be found in `results/2017-06-28_simulation_variable_genenumber`

The following figure shows the accuracy obtained by different machine learning methods. Note that the **generalized linear model** (glm) performed worse than other models. **Random forests** mantained the same accuracy along the different datasets (as well as elastic net with an alpha parameter of 0.5 and 0.1).

![](project_notebook_img/simulation_100_cells_ngenes-accuracy_per_model.png)

This same analysis will be performed using a range of genes from 1 to 3000.

![](project_notebook_img/simulation_100_cells_ngenes-100.png)
![](project_notebook_img/simulation_100_cells_ngenes-500.png)
![](project_notebook_img/simulation_100_cells_ngenes-1000.png)
![](project_notebook_img/simulation_100_cells_ngenes-2000.png)


# 12/07/2017

`prediction_simulation.Rmd` (commit [e18d026](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/e18d026f4aabbf91f31cadfe3132fa59969c5c40/bin/prediction_simulation.Rmd)).

Some of the updates so far in this version include:

- Use of new `TrainPredict()` function from the `scPrediction` package to optimize code and reduce redundancy
- Centers and scales the data before training model to improve accuracy. As noted by [ Jason Brownlee](http://machinelearningmastery.com/pre-process-your-dataset-in-r/)

> It is an easy step to forget or skip over and often has a huge impact on the accuracy of your final models

- Avoids complex call of `foreach()` function by creating output list within loop
- Adds new models such as:
  + **cart** (Classification and regression tree)
  + **bayesglm** (Bayesian generalized model)
  + **nnet** (Neural network)
  + **dnn** (Stacked AutoEncoder Deep Neural Network)
- Runs **elastic net** algorithm using the `caret` package. This way only the best model for a **alpha between 0.1 and 0.9** is selected
- `Ridge` and `lasso` methods are run outside `caret` using the `glmnet` package
- Creates accuracy plot per dataset and method. Plots boxplots in specified order.
- Saves accuracy summaries (to reproduce plots in case something goes wrong)

Notes:

> Least angle regression was not included in the analysis since it is not implemented for classification tasks in R. See [Stack Exchange post](https://stats.stackexchange.com/questions/48651/least-angle-regression-packages-for-r-or-matlab)

This version was run with the following parameters:

`Splatter`:

```R
nGenes <- c(2, 3, 4, 5, 10, 20, 30, 50, 100, 200, 500, 1000, 1500, 2000)
groupCells <- c(50, 50)
method <-  "groups"
seed <-  10
de.prob <- 1
```

`de.prob` refers to 

> Probability that a gene is differentially expressed in a group

Setting `de.prob` to 1 forces all genes to be differentially expressed.

After running the pipeline the following warnings stopped the process.

From the MARS method:

```R
earth glm Group2: did not converge after 25 iterations
Something is wrong; all the Accuracy metric values are missing:
```

And from `caret` function `preProcess()`.

```R
Warning in preProcess.default(thresh = 0.95, k = 5, freqCut = 19, uniqueCut = 10,  :
  These variables have zero variances: Gene92, Gene127, Gene193, Gene279, Gene343, Gene457, Gene532, Gene543, Gene589, Gene593, Gene694, Gene739, Gene773, Gene842, Gene906, Gene910, Gene934, Gene976
```

This last warning happens when trying to center and scale the data prior to training the model.

As noted by [Thiago G. Martins](https://tgmstat.wordpress.com/2014/03/06/near-zero-variance-predictors/), this error may be due to zero and near-zero predictors:

> Constant and almost constant predictors across samples (called zero and near-zero variance predictors in [1], respectively) happens quite often. One reason is because we usually break a categorical variable with many categories into several dummy variables. Hence, when one of the categories have zero observations, it becomes a dummy variable full of zeroes.

As some predictors after center and scaled have zero variances:

> This kind of predictor is not only non-informative, it can break some models you may want to fit to your data


As this error seemed to appeared when using 1000 genes, we will reduce the number of genes to 500 and see if we get the zero predictors.

Two explanations are possible:

1. This error did not happen before maybe due to the differential expression parameter (`de.prob`). As all genes are differentially expressed, the variance in the dataset may be reduced due to constant predictors with similar values
2. Centering and scaling affects estimations in some ML methods such as `glm`


# 13/07/2017

`prediction_simulation.Rmd` (commit [c187835](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/c187835dc1d6eae284cb2291e503ace58011db1d/bin/prediction_simulation.Rmd)) was run changing the number of genes to be analyzed.

```R
nGenes <- c(2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500)
```

![](project_notebook_img/simulation_100_cells_variable-nGenes_2017-07-13.png)


As shown in the following plot, in general all curves (number of genes vs. accuracy) **tend to asymptote with four genes** to make the predictions.

- **glm** (generalized linear model) tends to perform as good as other models. However, when the number of genes increases (from 50) the accuracy decreases.
- **rpart** (Classification and regression tree) presents more variability whern compared to other models as the number of genes increases
- **rf** (random forests) performs worse eith a low number of genes but the accuracy is more stable as more genes are added
- **dnn** (Stacked AutoEncoder Deep Neural Network) presents an oscillatory behaviour

![](project_notebook_img/simulation_100_cells_ngenes-accuracy_per_model_fixed_2017-07-13.png)

![](2017-07-13_simulation_variable_genenumber/simulation_100_cells_ngenes-2.png)
![](2017-07-13_simulation_variable_genenumber/simulation_100_cells_ngenes-3.png)
![](2017-07-13_simulation_variable_genenumber/simulation_100_cells_ngenes-4.png)

# 14/07/2017

In a meeting with @joseph, we proposed the following tasks to do:

## TODO

### Simulations

- Reduce the number of differentially expresed genes to evaluate to 100
- Test the accuracy as the number of cells increases
- Test the accuracy changing the number of cells in each cluster (ratio)

### Real datasets

- Build a predictor for the [Atlas of human blood dentritic cells and monocytes](https://portals.broadinstitute.org/single_cell/study/atlas-of-human-blood-dendritic-cells-and-monocytes)
- Build a predictor using bulk RNA-seq data to predict tissue origin of single cells


https://stackoverflow.com/questions/29995184/glmnet-error-for-logistic-regression-binomial

```R
one multinomial or binomial class has 1 or 0 observations; not allowed
```

https://stats.stackexchange.com/questions/232228/cross-validation-for-uneven-groups-using-cv-glmnet

```R
one multinomial or binomial class has fewer than 8  observations; dangerous ground
```

https://stats.stackexchange.com/questions/232228/cross-validation-for-uneven-groups-using-cv-glmnet

# 18/07/17


`eda_blood_atlas.Rmd` see commit ([7acc1c1](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/commit/7acc1c11ecdb944d4583a82fb05a80bd0d89ed2d?diff=unified)), was created to explore the data from the [Atlas of human blood dentritic cells and monocytes](https://portals.broadinstitute.org/single_cell/study/atlas-of-human-blood-dendritic-cells-and-monocytes). The following table shows the cell types and the number of cells in this dataset.

|Cell type | Frequency|
|:---------|---------:|
|DC1       |       165|
|DC2       |        94|
|DC3       |       107|
|DC4       |       173|
|DC5       |        30|
|DC6       |       173|
|Mono1     |       163|
|Mono2     |       122|
|Mono3     |        31|
|Mono4     |        20|

- **Number of genes**: 26593 
- **Number of cells**: 1078


The authors also provide a list of *discriminant genes* to differentiate between cell types. The following table shows the number of discriminant genes for each cell type.

|Cell type | Number of discriminant genes|
|:---------|----------------------------:|
|DC1       |                          112|
|DC2       |                           31|
|DC3       |                           59|
|DC4       |                          342|
|DC5       |                           85|
|DC6       |                          390|
|Mono1     |                          104|
|Mono2     |                           22|
|Mono3     |                          186|
|Mono4     |                          149|


A principal component analysis was performed using all genes included in the dataset. The following plot shows the first 3 components. Notice that PC2 and PC3 explain the separation of cell types and they account for 2.07 and 1.09 percent of the variance. This may suggest that only a small subset of genes may be necessary to predict cell type (e.g. discriminant genes proposed by the authors).

![](./2017-07-17_blood_atlas_prediction/blood_atlas_pca.png)

For this exploratory analysis, only the top 10 of discriminant genes for each cell type will be considered. The genes with the best **p-values** were selected as the top 10. In the following plot, the gene expression distribution for each discriminant gene is shown for every cell type. Some genes do explain only one cell type (e.g. **IL2RB** which is preferentially expressed in monocites 4). Other genes seem to separete monocites from dentritic cells such as **BC013828** and **LAIR2**

![](./2017-07-17_blood_atlas_prediction/top10_total_90_discriminant_genes.png)

A PCA was performed again only using the top 10 discriminant genes of each cell type (in total **90 genes** as 10 are shared between cell types). Cell type information is still preserved in this components.
![](./2017-07-17_blood_atlas_prediction/pca_top10_total_90_discriminant_genes.png)

# 19/07/17

`prediction_simulation_ratio.Rmd` (see commit ([7d380f4](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/7d380f497b20b52a1a08ea9841d321b6e2851055/bin/prediction_simulation_ratio.Rmd)) was created to see the behaviour of prediction models depending on the proportion of cell between two groups. This script was run using the followin number of genes and sample sizes:

```R
n.genes <- c(2, 3, 4, 5, 10, 20, 30, 40, 50, 100)
n.cells <- 10000
cell.ratio <- seq(0.1, 0.4, 0.1)
```

PCA plots of simulated data may be found [here](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/master/results/2017-07-14_simulation_variable_ratio_and_genes/PCA_simulations.pdf).


In the following plot, the accuracy is shown as a fuction of the number of genes. For each model, 4 subpanels are displayed indicating the proportion of cells between the simulated group 1 and 2. No appreciable differencre in performance is observed across all models. Some predictions of `rpart` fall to 70% of accuracy when the number of genes is 3 and the proportion of cells is 30/70.

![](2017-07-14_simulation_variable_ratio_and_genes/simulation-accuracy_per_model.png)

The following plot shows the same results but the scale is **independent** for each model. A **drop in accuracy is observable when 2 or 3 genes are used** as features.

![](2017-07-14_simulation_variable_ratio_and_genes/simulation-accuracy_per_model_scale_free.png)

Plots for all simulations may be found [here](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/tree/master/results/2017-07-14_simulation_variable_ratio_and_genes).

# 26/07/2017

`prediction_human_cell_atlas.Rmd` (commit [a865cfd](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/a865cfd2dcce48385afcbb9f80531c25b4e93555/bin/prediction_human_cell_atlas.Rmd)) was run to predict all 10 cell types from the **human blood atlas of monocytes and dentritic cells**.


The discriminant genes were selected based in the following criteria.

For each cell type:

- Genes were selected if the reported AUC value was greater than `0.7`. From those genes, the top `20` were considered for the prediction step.

```R
disc.genes %>%
  group_by(cluster.id) %>%
  filter(auc.value > 0.7) %>% 
  top_n(n = 20, auc.value) -> top.n
```

The table below shows the number of genes considered as features:

|Cell type  | Number genes    |
|:----------|----------------:|
|DC1        | 22              |
|DC2        | 19              |
|DC3        | 20              |
|DC4        | 20              |
|DC5        | 21              |
|DC6        | 22              |
|Mono1      | 20              |
|Mono2      | 20              |
|Mono3      | 20              |
|Mono4      | 20              |

The following plots how the acuraccy results obtained by prediction models for each cell type. 50 bootstrap replicates for each model for every cell type prediction were performed.


![](../results/2017-07-25_blood_atlas_prediction/accuracy_per_model_DC1.png)
![](../results/2017-07-25_blood_atlas_prediction/accuracy_per_model_DC2.png)
![](../results/2017-07-25_blood_atlas_prediction/accuracy_per_model_DC3.png)
![](../results/2017-07-25_blood_atlas_prediction/accuracy_per_model_DC4.png)
![](../results/2017-07-25_blood_atlas_prediction/accuracy_per_model_DC5.png)
![](../results/2017-07-25_blood_atlas_prediction/accuracy_per_model_DC6.png)
![](../results/2017-07-25_blood_atlas_prediction/accuracy_per_model_Mono1.png)
![](../results/2017-07-25_blood_atlas_prediction/accuracy_per_model_Mono2.png)
![](../results/2017-07-25_blood_atlas_prediction/accuracy_per_model_Mono3.png)
![](../results/2017-07-25_blood_atlas_prediction/accuracy_per_model_Mono4.png)

- In general, `svmPoly` (support vector machines with a polynomial kernel) performed better or roughly equal than the other models for all predictions
- For the toughest classification tasks (DC2 and DC3 prediction), `dnn` (Stacked AutoEncoder Deep Neural Network) performs significantly better than for easier classifications when all models perform similar.


# 27/07/2017

`prediction_simulation_ncells.Rmd` (commit [a6e8085](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/a6e808508441e438559e217874a1a06a882b15f8/bin/prediction_simulation_ncells.Rmd)) was run with the following parameters:

```R
n.genes <- c(2, 3, 4, 5, 10, 20, 30, 40, 50, 100)
n.cells <- c(100, 500, 1000, 5000, 10000, 50000)
```

The following plot shows the accuraccy performance of all tested models as the number of genes included in each simulated dataset changes (omitting the number of cells).


![](../results/2017-07-14_simulation_variable_cells_and_genes/simulation-accuracy_per_model.png)

# 28/07/2017

`prediction_hipsc_scprediction.Rmd` (commit [8ee70f8](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/8ee70f815448cb671cfa6cd1203f472bb4c6318f/bin/prediction_hipsc_scprediction.Rmd)) was run to predict cluster 1.

The following plot shows the results obtained when predicting cluster 1 versus other clusters (2, 3 and 4). The accuracy is almost the same across all models.

![](../results/2017-07-26_hiPSC_Analysis/accuracy_per_model_1.png)

# 31/07/2017

In a meeting with Joseph the following approach for prediction cell clusters was proposed:


1. Classify data into two clusters: cluster of interest and other clusters (e.g. cluster 1 vs. cluster 2, 3 and 4)
1. Perform a principal component analysis (PCA) using all genes
2. For all eigenvectors, test if there is a significant difference between the values for the cluster of interest and the other clusters using a Mann-Whitney test
3. Adjust p-values using a bonferroni correction
4. Keep those eigenvectors with adjusted p-values below `0.05`
5. Use significant eigenvectors as features for prediction


`prediction_human_blood_atlas_eigenvectors.Rmd` (commit [62039f6](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/62039f609993ada1779702ac4da0f08ca9f90f14/bin/prediction_human_blood_atlas_eigenvectors.Rmd)) was run using the previous methodology.


The following plots show the first three components for each cell type.

## PCA

### DC1
![](../results/2017-07-31_blood_atlas_prediction/pca_DC1.png)

### DC2
![](../results/2017-07-31_blood_atlas_prediction/pca_DC2.png)

### DC3
![](../results/2017-07-31_blood_atlas_prediction/pca_DC3.png)

### DC4
![](../results/2017-07-31_blood_atlas_prediction/pca_DC4.png)

### DC5
![](../results/2017-07-31_blood_atlas_prediction/pca_DC5.png)

### DC6
![](../results/2017-07-31_blood_atlas_prediction/pca_DC6.png)

### Mono1
![](../results/2017-07-31_blood_atlas_prediction/pca_Mono1.png)

### Mono2
![](../results/2017-07-31_blood_atlas_prediction/pca_Mono2.png)

### Mono3
![](../results/2017-07-31_blood_atlas_prediction/pca_Mono3.png)

### Mono4
![](../results/2017-07-31_blood_atlas_prediction/pca_Mono4.png)


## Prediction

The following plots show the accuracy results obtained using the significant components as features. In the right side, a table is shown with the selected **eigenvectors**, the **adjusted p-value** and the **percentage of variance** explained. 

- 10 bootstrap replicates
- 5 cross-validations per replicate

![](../results/2017-07-31_blood_atlas_prediction/accuracy_per_model_DC1.png)
![](../results/2017-07-31_blood_atlas_prediction/accuracy_per_model_DC2.png)
![](../results/2017-07-31_blood_atlas_prediction/accuracy_per_model_DC3.png)
![](../results/2017-07-31_blood_atlas_prediction/accuracy_per_model_DC4.png)
![](../results/2017-07-31_blood_atlas_prediction/accuracy_per_model_DC5.png)
![](../results/2017-07-31_blood_atlas_prediction/accuracy_per_model_DC6.png)
![](../results/2017-07-31_blood_atlas_prediction/accuracy_per_model_Mono1.png)
![](../results/2017-07-31_blood_atlas_prediction/accuracy_per_model_Mono2.png)
![](../results/2017-07-31_blood_atlas_prediction/accuracy_per_model_Mono3.png)
![](../results/2017-07-31_blood_atlas_prediction/accuracy_per_model_Mono4.png)

# 01/08/2017

Using the eigenvectors as features as in the previous analysis, the same approach was applied to the hiPSC dataset.

`prediction_hipsc_eigenvectors.Rmd` (commit [a3bb4ab](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/a3bb4abefb1f0cf0d11b6e93b1e96e57c0ce408c/bin/prediction_hipsc_eigenvectors.Rmd)) was run.



## PCA

### Cluster 1
![](../results/2017-08-01_hipsc_eigenvectors/pca_1.png)

### Cluster 2
![](../results/2017-08-01_hipsc_eigenvectors/pca_2.png)

### Cluster 3
![](../results/2017-08-01_hipsc_eigenvectors/pca_3.png)

### Cluster 4
![](../results/2017-08-01_hipsc_eigenvectors/pca_4.png)

## Prediction

The following plots show the accuracy results obtained using the significant components as features. In the right side, a table is shown with the selected **eigenvectors**, the **adjusted p-value** and the **percentage of variance** explained. 

- 10 bootstrap replicates
- 5 cross-validations per replicate

### Cluster 1
![](../results/2017-08-01_hipsc_eigenvectors/accuracy_per_model_1.png)

### Cluster 2
![](../results/2017-08-01_hipsc_eigenvectors/accuracy_per_model_2.png)

### Cluster 3
![](../results/2017-08-01_hipsc_eigenvectors/accuracy_per_model_3.png)

### Cluster 4
![](../results/2017-08-01_hipsc_eigenvectors/accuracy_per_model_4.png)


# 04/08/2017

`hipsc_PCA_cumulative.Rmd` (commit [276aa90](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/276aa9039faec7de1c9200a665d15a6e29919e55/bin/hipsc_PCA_cumulative.Rmd)) was run to determine if the prediction accuracy "plateus" at a certain number of significant components used as features. Below four plots with the accuracy results of the HiPSC dataset are shown.


![](../results/2017-08-04_hipsc_eigenvectors/accuracy_per_model_1.png)
![](../results/2017-08-04_hipsc_eigenvectors/accuracy_per_model_2.png)
![](../results/2017-08-04_hipsc_eigenvectors/accuracy_per_model_3.png)
![](../results/2017-08-04_hipsc_eigenvectors/accuracy_per_model_4.png)

# 08/08/2017

`hipsc_PCA_cumulative.Rmd` (commit [16d23b3](https://github.com/IMB-Computational-Genomics-Lab/SingleCell_Prediction/blob/16d23b30ccaf145922bd7bb0817777c03158f1af/bin/hipsc_PCA_cumulative.Rmd)) was run to observe how the accuracy flattens when the top 25 PCs -ranked by their p-value obtained from Mann-Whitney test comparing the PC-value distributions of the cluster evaluated (e.g cluster 1) and the remaining clusters (e.g. clusters 2, 3 and 4)-. In this analysis, the sensitivity and specificity were also included. Note that non-significant PCs are also included in the plots below as the top 25 were selected without taking into account an alpha value threshold.

![](../results/2017-08-08_hipsc_eigenvectors/metrics_cluster_1.png)
![](../results/2017-08-08_hipsc_eigenvectors/metrics_cluster_2.png)
![](../results/2017-08-08_hipsc_eigenvectors/metrics_cluster_3.png)
![](../results/2017-08-08_hipsc_eigenvectors/metrics_cluster_4.png)

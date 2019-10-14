# tested on R 3.4.3
# setup required libraries and constants
library(scater)  # tested on version 1.6.3,    install from Bioconductor: source("https://bioconductor.org/biocLite.R"); biocLite("scater")
library(xgboost) # tested on version 0.6.4.1, install from CRAN: install.packages("xgboost")
library(igraph)  # tested on version 1.2.1,   install from CRAN: install.packages("igraph")

library(tidyverse)

BREAKS=c(-1, 0, 1, 6, Inf)
nFeatures = 100

# 1. Load datasets in scater format: loaded files expected to contain "Large SingleCellExperiment" object
source = readRDS("data/2018-04-15_pancreas_baron/baron-human.rds")
target = readRDS("data/2018-04-15_pancreas_xin/xin.rds")


ds1 = t(exprs(source)) 
ds2 = t(exprs(target))
colnames(ds2) <- colnames(ds2) %>% str_remove("__.*$") -> rownames(target)
sourceCellTypes = colData(source)[,"cell_type1"]
sourceCellTypes <- as.factor(ifelse(sourceCellTypes %in% c("alpha", "beta", "delta", "gamma"), 
                                    as.character(sourceCellTypes), "other"))

# 2. Unify sets, excluding low expressed genes
source_n_cells_counts = apply(exprs(source), 1, function(x) { sum(x > 0) } )
target_n_cells_counts = apply(exprs(target), 2, function(x) { sum(x > 0) } )
common_genes = intersect( rownames(source)[source_n_cells_counts>10], 
                          rownames(target)[target_n_cells_counts>10]
)
remove(source_n_cells_counts, target_n_cells_counts)
ds1 = ds1[, colnames(ds1) %in% common_genes]
ds2 = ds2[, colnames(ds2) %in% common_genes]
ds = rbind(ds1[,common_genes], ds2[,common_genes])
isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
remove(ds1, ds2)

# 3. Highest mean in both source and target
topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]

# 4. Highest mutual information in source
topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),sourceCellTypes,method = "nmi") }), decreasing = T))

# 5. Top n genes that appear in both mi and avg
selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )

# 6. remove correlated features
tmp = cor(ds[,selectedFeatures], method = "pearson")
tmp[!lower.tri(tmp)] = 0
selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
remove(tmp)

# 7,8. Convert data from continous to binned dummy vars
# break datasets to bins
dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
# use only bins with more than one value
nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
# convert to dummy vars
ds = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
remove(dsBins, nUniq)

# 9. Classify
train = runif(nrow(ds[isSource,]))<0.8
# slightly different setup for multiclass and binary classification
if (length(unique(sourceCellTypes)) > 2) {
  xg=xgboost(data=ds[isSource,][train, ] , 
             label=as.numeric(sourceCellTypes[train])-1,
             objective="multi:softmax", num_class=length(unique(sourceCellTypes)),
             eta=0.7 , nthread=5, nround=20, verbose=0,
             gamma=0.001, max_depth=5, min_child_weight=10)
} else {
  xg=xgboost(data=ds[isSource,][train, ] , 
             label=as.numeric(sourceCellTypes[train])-1,
             eta=0.7 , nthread=5, nround=20, verbose=0,
             gamma=0.001, max_depth=5, min_child_weight=10)
}

# 10. Predict
predictedClasses = predict(xg, ds[!isSource, ])
predictedClasses <- levels(sourceCellTypes)[predictedClasses + 1]

targetCellTypes = colData(target)[,"cell_type1"]

props <- table(predictedClasses, targetCellTypes) %>% 
  as.data.frame() %>% 
  spread(key = "targetCellTypes", value = "Freq") %>% 
  column_to_rownames("predictedClasses")


mapply(function(x,y) x/y, props, colSums(props)) %>% 
  `rownames<-`(rownames(props)) %>% 
  round(2)

       alpha alpha.contaminated beta beta.contaminated delta delta.contaminated
alpha  0.97               0.30 0.51              0.19  0.71               0.22
beta   0.00               0.02 0.44              0.10  0.06               0.00
delta  0.00               0.00 0.00              0.00  0.10               0.11
gamma  0.00               0.00 0.00              0.00  0.00               0.00
other  0.03               0.68 0.05              0.71  0.12               0.67
       gamma gamma.contaminated
alpha  0.61               0.38
beta   0.00               0.00
delta  0.01               0.00
gamma  0.00               0.00
other  0.38               0.62


       alpha alpha.contaminated beta beta.contaminated delta delta.contaminated
alpha   855                 18  239                 6    35                  2
beta      1                  1  209                 3     3                  0
delta     0                  0    0                 0     5                  1
gamma     0                  0    0                 0     0                  0
other    30                 41   24                22     6                  6
       gamma gamma.contaminated
alpha    52                  3
beta      0                  0
delta     1                  0
gamma     0                  0
other    32                  5


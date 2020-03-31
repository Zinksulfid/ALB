#neue normalisierungsfunktion

normalize_vsn_Anna <- function (se, spikeins) {
  #spikeins = 100:200
  spfit = vsn2(se[spikeins,], lts.quantile=1)
  se_vsn = predict(spfit, newdata=se)
  return(se_vsn)
}

#DEP script zum normalisieren

library("DEP")
library("vsn")
library("ggplot2")
library("dplyr")
Interactome <- read.csv("Fly_Interactome_curated.csv", header=TRUE, sep=",")
expdesign <- read.table("NFS1_expdesign.txt", header = TRUE, sep="\t")
proteinGroups <- read.table("proteinGroups.txt", header = TRUE, sep="\t")

data <- filter(proteinGroups, Reverse != "+", Potential.contaminant != "+")
dim(data)
colnames(data)
data$Gene.names %>% duplicated() %>% any()
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
data$name %>% duplicated() %>% any()

LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se <- make_se_parse(data_unique, LFQ_columns)
#fehlerbehebung
LFQ_columns
dim(data_unique)
data_unique[102:117]

plot_frequency(data_se)
data_filt <- filter_missval(data_se, thr = 0)
plot_numbers(data_filt)
plot_coverage(data_filt)

data_norm <- normalize_vsn_Anna(data_filt, c(498, 631, 286, 156, 232 )) #ohne wt
  #data_norm <- normalize_vsn_Anna(data_filt, c(590, 1084, 461, 541, 491)) #mit wt
  #data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm)

plot_missval(data_filt)
plot_detect(data_filt)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
plot_imputation(data_norm, data_imp)

data_diff <- test_diff(data_imp, type = "control", control = "wt")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

plot_pca(dep, x = 1, y = 2, n = 374, point_size = 4)


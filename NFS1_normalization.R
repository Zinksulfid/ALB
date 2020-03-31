#neue normalisierungsfunktion

normalize_vsn_anna <- function(se, spikeins) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))  # Variance stabilization transformation on assay data
  se_vsn <- se
  spfit = vsn::vsnMatrix(2 ^ assay(se_vsn)[spikeins,], minDataPointsPerStratum = 1)
  assay(se_vsn) <- vsn::predict(spfit, 2 ^ assay(se_vsn))
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

LFQ_columns <- grep("LFQ.", colnames(data_unique))# get LFQ column numbers
data_se <- make_se_parse(data_unique, LFQ_columns)
#fehlerbehebung
LFQ_columns
dim(data_unique)
data_unique[102:117]

plot_frequency(data_se)
data_filt <- filter_missval(data_se, thr = 0)
plot_numbers(data_filt)
plot_coverage(data_filt)


#data_norm <- normalize_vsn_anna(se=data_filt, spikeins=c("CG30022", "Jafrac1", "RpS19a", "Hsp60", "alpha-Spec" )) #ohne wt
#data_norm <- normalize_vsn_anna(data_filt, c( "skap", "Tina-1", "EG:87B1.3", "Q8MT18", "Q7PLE6")) #mit wt
data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm)

plot_missval(data_filt)
plot_detect(data_filt)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
plot_imputation(data_norm, data_imp)

data_diff <- test_diff(data_imp, type = "control", control = "wt")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

plot_pca(dep, x = 1, y = 2, n = 374, point_size = 4)
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
plot_volcano(dep, contrast = "F_vs_wt", label_size = 2, add_names = TRUE)

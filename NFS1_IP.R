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
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)
plot_volcano(dep, contrast = "F_vs_wt", label_size = 2, add_names = TRUE)

#barplots
plot_single(dep, proteins = "CG12264", type = "centered")
plot_single(data_se, proteins = "CG12264", type = "centered")
install.packages("ggplot2")
library(ggplot2)
install.packages("stats")
library(stats)
#daten auslesen
LFQ_columns <- grep("LFQ.", colnames(data_unique))
data_unique %>%
  select(Gene.names) -> data_lfq
data <- cbind(data_unique[LFQ_columns], data_lfq)
data_final <- subset(data, Gene.names == "CG12264")
#mittelwerte berechnen
data_final %>%
  select(starts_with("LFQ.intensity.wt")) -> wt 
  mean_wt <-rowMeans(wt)
  sd_wt <-apply(wt, 1, sd)
data_final %>%
  select(starts_with("LFQ.intensity.FLAG")) -> FLAG
  mean_FLAG <-rowMeans(FLAG)
  sd_FLAG <-apply(FLAG, 1, sd)
data_final %>%
  select(starts_with("LFQ.intensity.K")) -> K 
  mean_K <-rowMeans(K)
  sd_K <-apply(K, 1, sd)
data_final %>%
  select(starts_with("LFQ.intensity.wt")) -> F 
  mean_F <-rowMeans(F)
  sd_F <-apply(F,1,sd)

genotype<- c("wt", "FLAG", "K", "F")
mean <- c(mean_wt, mean_FLAG, mean_K, mean_F)
sd <- c(sd_wt, sd_FLAG, sd_K, sd_F)
data_plot_vektor <- cbind(genotype, mean, sd)
data_plot<-data.frame(data_plot_vektor)
d$mean <- as.numeric(as.character(d$mean))
d$genotype <- factor(d$genotype, level = c("wt", "Flag", "K", "F"))
d$sd <- as.numeric(as.character(d$sd))

#plotten
ggplot(data_plot, aes(x = genotype, y = mean)) +
  
  geom_bar(stat="identity", position=position_dodge())+
  
  scale_fill_brewer(type = "seq", palette = 3, direction = 1, aesthetics = "fill" ) + 
  
  geom_errorbar(aes(ymax = mean + sd, ymin= mean - sd), position = position_dodge(), width=0.4) 
  

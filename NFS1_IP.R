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
#install.packages("ggplot2")
#library(ggplot2)
#install.packages("stats")
#library(stats)
#daten auslesen
LFQ_columns <- grep("LFQ.", colnames(data_unique))
data_unique %>%
  select(Gene.names) -> data_lfq
data <- cbind(data_unique[LFQ_columns], data_lfq)
data_final <- subset(data, Gene.names == "CG12264")
#mittelwerte berechnen
data_final %>%
  select(starts_with("LFQ.intensity.wt")) -> wt 
  mean_wt <-log2(rowMeans(wt))
  sd_wt <-log2(apply(wt, 1, sd))
data_final %>%
  select(starts_with("LFQ.intensity.FLAG")) -> FLAG
  mean_FLAG <-log2(rowMeans(FLAG))
  sd_FLAG <-log2(apply(FLAG, 1, sd))
data_final %>%
  select(starts_with("LFQ.intensity.K")) -> K
  mean_K <-log2(rowMeans(K))
  sd_K <-log2(apply(K, 1, sd))
data_final %>%
  select(starts_with("LFQ.intensity.F")) -> F 
  mean_F <-log2(rowMeans(F))
  sd_F <-log2(apply(F,1,sd))

genotype<- c("wt", "FLAG", "K", "F")
mean <- c(mean_wt, mean_FLAG, mean_K, mean_F)
sd <- c(sd_wt, sd_FLAG, sd_K, sd_F)
data_plot_vektor <- cbind(genotype, mean, sd)
data_plot<-data.frame(data_plot_vektor)
data_plot$mean <- as.numeric(as.character(data_plot$mean))
data_plot$genotype <- factor(data_plot$genotype, level = c("wt", "FLAG", "K", "F"))
data_plot$sd <- as.numeric(as.character(data_plot$sd))

#plotten
ggplot(data_plot, aes(x = genotype, y = mean)) +

  
  geom_bar(stat="identity", position=position_dodge())+
  
  scale_fill_brewer(type = "seq", palette = 3, direction = 1, aesthetics = "fill" ) + 
  
  geom_errorbar(aes(ymax = mean + sd, ymin= mean - sd), position = position_dodge(), width=0.4) 

#plotten für LFQ gegen intensitäten
data_final %>%
  select(starts_with("LFQ.intensity.wt")) -> wt_LFQ
data_final %>%
  select(starts_with("intensity.wt"))  -> wt
data_final %>%
  select(starts_with("LFQ.intensity.K")) -> K_LFQ
data_final %>%
  select(starts_with("intensity.K"))  -> K
data_final %>%
  select(starts_with("LFQ.intensity.FLAG")) -> FLAG_LFQ
data_final %>%
  select(starts_with("intensity.FLAG"))  -> FLAG 
# hier manuell übertragen!
x_data <- c(6266100000,        3552700000,        7829100000,        5912300000, 
            7865600000,       8293500000,       2.1979e+10,       1.8967e+10,
            9054900000, 8279800000, 1.6293e+10, 2.1106e+10)
y_data <- c(6.906e+09,    4070100000,     7.908e+09,    4686900000,
            1.1925e+10,    8.662e+09,   1.9115e+10,   1.7939e+10,
            7455500000, 7222200000, 1.5766e+10, 1.5468e+10)
genotyp <-c(rep("wt", length(x_data)/3),rep("K", length(x_data)/3),rep("FLAG", length(x_data)/3))
data <- cbind(x_data, y_data, genotyp)
data_plot <-data.frame(data)
colnames(data_plot)<- c("FLAG_LFQ_Data", "FLAG_Data")
ggplot(data_plot, aes(x = "FLAG_LFQ_Data", y = "FLAG_Data"), fill=genotyp)+ 
  geom_jitter()


# alle werte plotten mit facet_wrap
library("DEP")
library("dplyr")
library("ggplot2")
library("stringr")
proteinGroups <- read.table("proteinGroups.txt", header = TRUE, sep="\t")
data <- filter(proteinGroups, Reverse != "+", Potential.contaminant != "+")

data$Gene.names %>% duplicated() %>% any()
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
data$name %>% duplicated() %>% any()
LFQ <- grep("LFQ.", colnames(data_unique))
data_LFQ <- cbind(data_unique[LFQ], data_unique["Protein.IDs"])
data_LFQ %>%
  gather( "S1", "LFQ", 1:16) -> dLFQ
#replace <- 
#dLFQ <- str_replace(dLFQ[S1],)
dLFQ %>%
  mutate(ID.unique = paste(Protein.IDs, S1)) -> dL

intensity <- grep("Intensity.", colnames(data_unique))
data_intense <- cbind(data_unique[intensity], data_unique["Protein.IDs"])
data_intense %>%
  gather( "sample", "intense", 1:16) -> dint
#dint <- str_replace(dLFQ[S1],)
dint %>%
  mutate(ID.unique = paste(Protein.IDs,sample)) -> dI

data_plot <- full_join(dL, dI, by="ID.unique")

ggplot(data_plot, aes(x = LFQ, y = intense)) +
  geom_point()+
  geom_abline()+
  facet_wrap( ~ sample, ncol=4, scales="free") 
  

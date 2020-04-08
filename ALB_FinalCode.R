#open libraries 
library("dplyr")
library("ggplot2")
library("DEP")
library("tidyr")
library("stringr")
library("org.Dm.eg.db")
library("ggrepel")
#necessary functions
#functions for evaluating spikeins
#criteria R^2
Regression <- function (prot) {
  #data for each genotype sorted after the prot (=protein of a row)-> ends up with a file with intensities for the protein (=prot), NFS1 and corresponding name
  d_LFQ %>% 
    select(matches(".FLAG.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs == prot | Majority.protein.IDs=="Q9VKD3") %>%
    gather("prot", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_FLAG
  colnames(data_FLAG) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches(".K.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==prot | Majority.protein.IDs=="Q9VKD3") %>%
    gather("prot", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_K
  colnames(data_K) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches("F.$") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==prot | Majority.protein.IDs=="Q9VKD3") %>%
    gather("prot", "CG12264", 1:4)%>%
    spread("Majority.protein.IDs", "CG12264") -> data_F
  colnames(data_F) <- c("LFQ", "NFS1", "Protein")
  
  #linear regression model for prot against NFS1 and calculating R^2-value
  model_FLAG <- lm(NFS1 ~ Protein, data_FLAG)
  R_FLAG<-summary(model_FLAG)$r.squared
  model_K <- lm(NFS1 ~ Protein, data_K)
  R_K<-summary(model_K)$r.squared
  model_F <- lm(NFS1 ~ Protein, data_F)
  R_F<-summary(model_F)$r.squared
  
  #calculating R^2-value of all genotypes
  vec <- c( R_FLAG, R_K, R_F)
  mean_R <- mean(vec)
  
  
  return(mean_R)
}

#criteria slope of linear regression
slope <- function (prot) {
  #data for each genotype sorted after the prot (=protein of a row)-> ends up with a file with intensities for the protein (=prot), NFS1 and corresponding name
  d_LFQ %>% 
    select(matches(".FLAG.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs == prot | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Prot", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_FLAG
  colnames(data_FLAG) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches(".K.") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==prot | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Prot", "CG12264", 1:4) %>%
    spread("Majority.protein.IDs", "CG12264") -> data_K
  colnames(data_K) <- c("LFQ", "NFS1", "Protein")
  d_LFQ %>% 
    select(matches("F.$") | ends_with("Majority.protein.IDs") ) %>%
    subset(Majority.protein.IDs ==prot | Majority.protein.IDs=="Q9VKD3") %>%
    gather("Prot", "CG12264", 1:4)%>%
    spread("Majority.protein.IDs", "CG12264") -> data_F
  colnames(data_F) <- c("LFQ", "NFS1", "Protein")
  
  
  #linear regression model for prot against NFS1 and calculating the slope
  model_FLAG <- lm(NFS1 ~ Protein, data_FLAG)
  S_FLAG<-summary(model_FLAG)$coefficients[2] 
  model_K <- lm(NFS1 ~ Protein, data_K)
  S_K<-summary(model_K)$coefficients[2] 
  model_F <- lm(NFS1 ~ Protein, data_F)
  S_F<-summary(model_F)$coefficients[2] 
  
  #returing a boolean, depending on the slope -> 1 negative slope = FALSE,  all positive slope = TRUE
  if(S_F > 0 #& S_wt >0 
     & S_K > 0 & S_FLAG >0) {
    S=TRUE
  }
  else{
    S=FALSE
  }
  return(S)
}

#normalization
normalize_vsn_anna <- function(se, spikeins) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))# Variance stabilization transformation on assay data
  se_vsn <- se
  spfit = vsn::vsnMatrix(2 ^ assay(se_vsn)[spikeins,], minDataPointsPerStratum = 1)
  assay(se_vsn) <- vsn::predict(spfit, 2 ^ assay(se_vsn))
  return(se_vsn)
}
#read files
Interactome <- read.csv("Fly_Interactome_curated.csv", header=TRUE, sep=",")
data_int <- Interactome [,1]
proteinGroups <- read.table("proteinGroups.txt", header = TRUE, sep="\t")
data <- filter(proteinGroups, Reverse != "+", Potential.contaminant != "+")

#sort data 
#no more duplicates, and exclude rows with intensity-values = 0
#data$Gene.names %>% duplicated() %>% any()
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#data$name %>% duplicated() %>% any()
#exclude rows with LFQ.intensity-values = 0 for normalization
LFQ <- grep("LFQ.", colnames(data_unique))
data_LFQ <- cbind(data_unique["Majority.protein.IDs"], data_unique[LFQ])
data_LFQ[data_LFQ == 0] <- NA
d_LFQ <- drop_na(data_LFQ)

#get the spikeins for normalization
#iterate through the rows of d_LFQ and apply Regression-Funktion -> collect all R^2-values
R_values <- data.frame(1,2)
for(i in c(1:dim(d_LFQ)[1])){
  d<-d_LFQ[i,]
  prot <- d[,"Majority.protein.IDs"]
  #make sure, that NFS1 is not aligened against itself
  if(prot == "Q9VKD3"){
    R_values <- rbind(R_values, c(NA,NA))
  }
  else {
    mean_R <- Regression (prot) # Muss in Variable gespeichert werden
    #print(c(as.character(zeile), mean_R))
    R_values <- rbind(R_values, c(as.character(zeile), mean_R)) 
  }
}
colnames(R_values) <- c("zeilen", "R")

#iterate through the rows of d_LFQ and apply slope-Funktion -> collect all slope-values
S_values <- data.frame(1,2)
for(i in c(1:dim(d_LFQ)[1])){
  d<-d_LFQ[i,]
  prot <- d[,"Majority.protein.IDs"]
  #make sure, that NFS1 is not aligened against itself
  if(prot == "Q9VKD3"){
    S_values <- rbind(S_values, c(NA,NA))
  }
  else {
    #print(zeile)
    S <- slope(prot)
    # Muss in Variable gespeichert werden
    #print(c(as.character(zeile), mean_R))
    S_values <- rbind(S_values, c(as.character(zeile), S)) 
  }
}
colnames(S_values) <- c("zeilen", "S")

#one table with R^2-values and slope-values; 
data_S_R <- full_join (R_values, S_values, by = "zeilen")
#take only the ones with s= TRUE into account
data_S_R %>%
  subset(S == TRUE) -> data_sub
#sort the data.frame with a descending R^2 value and take only the ones with R^2 > 0.75 into account -> the names of proteins = spikeins
data_final <- arrange(data_sub, desc(R))
data_final$name <- data_unique[match(data_final$zeilen, data_unique$Majority.protein.IDs),]$name
d_F <- subset(data_final, R >0.75)
spikeins <- c(d_F [,4])

#prepare data for normalization
LFQ_columns <- grep("LFQ.", colnames(data_unique))
data_se <- make_se_parse(data_unique, LFQ_columns)
#plot_frequency(data_se)
data_filt <- filter_missval(data_se, thr = 0)
#plot_numbers(data_filt)
#plot_coverage(data_filt)

#normalization
data_norm <- normalize_vsn_anna(data_filt, spikeins) 
#plot_normalization(data_filt, data_norm)

#prepare data for plotting
#plot_missval(data_filt)
#plot_detect(data_filt)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
#plot_imputation(data_norm, data_imp)
data_diff <- test_diff(data_imp, type = "control", control = "wt")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

#differentplots by DEP
#plot_pca(dep, x = 1, y = 2, n = 374, point_size = 4)
#plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
#plot_heatmap(dep, type = "centered", kmeans = TRUE, 
#k = 6, col_limit = 4, show_row_names = FALSE,
#indicate = c("condition", "replicate"))
#plot_volcano(dep, contrast = "F_vs_wt", label_size = 2, add_names = TRUE)
#plot_single(dep, proteins = "CG12264", type = "centered")


#prepare data for volcano-Plot

data_results <- get_results(dep)
data_results$Gene.ID <- mapIds(org.Dm.eg.db, keys=str_replace_all(data_results$ID, ";.*", ""),
                               column="ENSEMBL", keytype="UNIPROT",  multiVals="first")
#get ratio (x-axis), p-val (y-axis) and gene-ID
#replace replicate number for better facet_wrap
#get information, if protein is mitochondrial, by comparing with interactome
data_results %>%
  dplyr::select(matches("name|ID|Gene.ID|.ratio|p.adj")) -> data_ratio 
  tidyr::gather("measure", "ratio", 3:8 ) %>%
  mutate(type = str_replace(measure, ".*wt_", "")) %>%
  mutate(measure = str_replace(measure, "_.*", "")) %>%
  spread(type, ratio) %>%
  mutate(mito = Gene.ID %in% Interactome$Gene.ID) -> data_rf
#plot data in a PDF 
# colour correspond to significant level, shape correspond to mitochondrial protein
#new x/y-label, horizontale linie bei x=0 und facet-wrap zur aufteilung nach genotype
cairo_pdf("volcano_plot.pdf", 15,5)
ggplot(data_rf, aes(x = ratio, y = -log10(p.adj),color=p.adj <0.05,  shape=mito ))+# shape= mitochondrial protein)) +
  geom_point()+#color = ifelse(data_plot$significant == TRUE, "orange", "grey50"))+#, shape= ifelse(data_plot$Mitochondrial_protein == TRUE, 15, 16))+
  facet_wrap(~ measure, ncol=3, scales="free")+
  xlab("log2 Fold change")+
  ylab("-log10 p-value")+
  geom_line(aes(x=0), colour="#990000")+
  theme_bw()+
  scale_color_manual(values=c("#999999", "#E69F00"))+
  geom_text_repel(data = subset(data_rf, p.adj <0.05 & mito==TRUE), aes(label=name),   segment.size  = 0.2,
                  segment.color = "grey50", size=1.5)
dev.off()



#datatable with significant proteins for STRING analysis
data_sigF <- subset(data_rf, p.adj <0.05 & measure == "F" ) 
data_sig <- dplyr::select(data_sigF, name) 
rownames(data_sig) <- c()
write.table(file = "sig.proteine_F.xlsx", data_sig, sep = "\t", quote = FALSE)

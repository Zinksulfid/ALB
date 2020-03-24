install.packages("ggplot2")
library("ggplot2")
# erster Versuch 

a <- c("wt", "OE", "Flag", "K", "F")
b <- c(1.000, 115.162, 71.477, 55.136, 66.766)
d <- rbind(a,b)
d <- data.frame(d)
colnames(d) <- a
rownames(d) <- c( "genotype", "transcriptlevel")
pp <- ggplot(d, aes(x = "genotype", y = "transcriptlevel")) +
  geom_bar()
pp

# (1) Ich hab rausgenommen, dass der ggplot in pp gespeichert wird und man das separat wieder auslesen muss. Das kann man später machen, wenn der Code funktioniert, gerade ist das zu umständlich
# (2) Wie du siehst, sind die Genotypen nicht in der richtigen Reihenfolge. Das kann man theoretisch mit factor(d$genotype, levels = d$genotype) lösen
# (3) Auch die Zahlen sind nicht fließend abgebildet, sondern scheinbar als Faktor eingelesen. Das kannst du überprüfen, indem du str(d) laufen lässt. Aha, d$transcriptlevel ist tatsächlich ein Faktor. Mach daraus eine Zahlenreihe mit as.numeric(as.character(x))

# Daraus ergibt sich folgender Lösungsvorschlag
a <- c("wt", "OE", "Flag", "K", "F")
b <- c(1.000, 115.162, 71.477, 55.136, 66.766)
d <- cbind(a,b) # Beachte, hier cbind statt rbind, weil ggplot in Spalten und nicht Zeilen arbeitet.
# colnames(d) <- a # Brauchts nicht unbedingt
colnames(d) <- c( "genotype", "transcriptlevel") # Hier dementsprechend geändert
d <- data.frame(d)
str(d)
d$transcriptlevel <- as.numeric(as.character(d$transcriptlevel))
d$genotype <- factor(d$genotype, level = c("wt", "OE", "Flag", "K", "F"))

ggplot(d, aes(x = genotype, y = transcriptlevel)) +
  geom_bar(stat="identity", position=position_dodge())

# weiterarbeit an farbe und errorbalken
a <- c("wt", "OE", "Flag", "K", "F")
b <- c(1.000, 115.162, 71.477, 55.136, 66.766)
e <- c(0.084344263, 2.142360052, 13.36486933, 2.244103961, 1.595136296)
d <- cbind(a,b,e) # Beachte, hier cbind statt rbind, weil ggplot in Spalten und nicht Zeilen arbeitet.
# colnames(d) <- a # Brauchts nicht unbedingt
colnames(d) <- c( "genotype", "transcriptlevel", "se") # Hier dementsprechend geändert
d <- data.frame(d)
str(d)
d$transcriptlevel <- as.numeric(as.character(d$transcriptlevel))
d$genotype <- factor(d$genotype, level = c("wt", "OE", "Flag", "K", "F"))
d$se <- factor(d$se, level = c(0.084344263, 2.142360052, 13.36486933, 2.244103961, 1.595136296))


ggplot(d, aes(x = genotype, y = transcriptlevel, fill=genotype)) +
  geom_bar(stat="identity", position=position_dodge())+
 scale_fill_brewer( type = "seq", palette = 1, direction = 1, aesthetics = "fill" )+
  geom_errorbar(aes(ymax = transcriptlevel + se, ymin= transcriptlevel - se), position = position_dodge())

#analog für iron-test (gleiche funktion mit anderen Werten)
a <- c("wt", "OE", "Flag", "K", "F")
b <- c(1,	1.113059025,	0.987477587,	0.924865842,	0.513058988)
e <- c(0.163055447,	0.110694525,	0.140618712,	0.273389083,	0.249686538)
d <- cbind(a,b,e) # Beachte, hier cbind statt rbind, weil ggplot in Spalten und nicht Zeilen arbeitet.
# colnames(d) <- a # Brauchts nicht unbedingt
colnames(d) <- c( "genotype", "transcriptlevel", "se") # Hier dementsprechend geändert
d <- data.frame(d)
str(d)
d$transcriptlevel <- as.numeric(as.character(d$transcriptlevel))
d$genotype <- factor(d$genotype, level = c("wt", "OE", "Flag", "K", "F"))
d$se <- factor(d$se, level = c(0.163055447,	0.110694525,	0.140618712,	0.273389083,	0.249686538))


ggplot(d, aes(x = genotype, y = transcriptlevel, fill=genotype)) +
  geom_bar(stat="identity", position=position_dodge())+
 scale_fill_brewer( type = "seq", palette = 1, direction = 1, aesthetics = "fill" )

# lösung für nfs1
a <- c("wt", "OE", "Flag", "K", "F")
b <- c(1.000, 115.162, 71.477, 55.136, 66.766)
e <- c(0.084344263, 2.142360052, 13.36486933, 2.244103961, 1.595136296)
d <- cbind(a,b,e) # Beachte, hier cbind statt rbind, weil ggplot in Spalten und nicht Zeilen arbeitet.
# colnames(d) <- a # Brauchts nicht unbedingt
colnames(d) <- c( "genotype", "transcriptlevel", "se") # Hier dementsprechend geändert
d <- data.frame(d)
str(d)
d$transcriptlevel <- as.numeric(as.character(d$transcriptlevel))
d$genotype <- factor(d$genotype, level = c("wt", "OE", "Flag", "K", "F"))
d$se <- as.numeric(as.character(d$se))


ggplot(d, aes(x = genotype, y = transcriptlevel, fill=genotype)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(
    type = "seq",
    palette = 3,
    direction = 1,
    aesthetics = "fill" ) + 
  geom_errorbar(aes(ymax = transcriptlevel + se, ymin= transcriptlevel - se), position = position_dodge())
  
# lösung für iron

a <- c("wt", "OE", "Flag", "K", "F")
b <- c(1,	1.113059025,	0.987477587,	0.924865842,	0.513058988)
e <- c(0.163055447,	0.110694525,	0.140618712,	0.273389083,	0.249686538)
f <- cbind(a,b)
d <- cbind(f,e)# Beachte, hier cbind statt rbind, weil ggplot in Spalten und nicht Zeilen arbeitet.
# colnames(d) <- a # Brauchts nicht unbedingt
colnames(d) <- c( "genotype", "transcriptlevel", "se") # Hier dementsprechend geändert
d <- data.frame(d)
str(d)
d$transcriptlevel <- as.numeric(as.character(d$transcriptlevel))
d$genotype <- factor(d$genotype, level = c("wt", "OE", "Flag", "K", "F"))
d$se <- d$se <- as.numeric(as.character(d$se))


ggplot(d, aes(x = genotype, y = transcriptlevel, fill=genotype)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(
    type = "seq",
    palette = 3,
    direction = 1,
    aesthetics = "fill" ) + 
  geom_errorbar(aes(ymax = transcriptlevel + se, ymin= transcriptlevel - se), position = position_dodge())

# plot für andere primer (Qtzl, LYRMA4, ISCU, NFS1)
genotype <- c("wt", "OE", "Flag", "K", "F", "wt", "OE", "Flag", "K", "F", "wt", "OE", "Flag", "K", "F", "wt", "OE", "Flag", "K", "F")
transcriptlevel <- c(1.000, 0.865, 0.838, 0.956, 1.350, 1.000, 45.601, 37.844, 0.809, 0.955, 1.000, 0,898
                     0.698, 1.041, 1.027, 1.000, 115.162, 71.477, 55.136, 66.766)
se <- c(0.17972206, 0.046230542, 0.387379952, 0.183992176, 0.133365007, 0.203818635,
       3.312987884, 11.61584221, 0.003496234, 0.131513385, 0.179068215, 0.060578249, 0.281034778,
       0.138713883, 0.111072547, 0.084344263, 2.142360052, 13.36486933, 2.244103961, 1.595136296)
primer <- c("ISCU", "ISCU", "ISCU", "ISCU", "ISCU", "Qtzl", "Qtzl", "Qtzl", "Qtzl", "Qtzl", "LYRMA4", "LYRMA4", "LYRMA4", "LYRMA4", "LYRMA4", "NFS1", "NFS1", "NFS1", "NFS1", "NFS1")
f <- cbind(genotype,primer)
h<- cbind(f,transcriptlevel)
d <- cbind(h,se)# Beachte, hier cbind statt rbind, weil ggplot in Spalten und nicht Zeilen arbeitet.
# colnames(d) <- a # Brauchts nicht unbedingt
colnames(d) <- c( "genotype", "transcriptlevel", "se") # Hier dementsprechend geändert
d <- data.frame(d)
str(d)
d$transcriptlevel <- as.numeric(as.character(d$transcriptlevel))
d$genotype <- factor(d$genotype, level = c("wt", "OE", "Flag", "K", "F"))
d$primer <- factor(d$primer, level = c("ISCU", "ISCU", "ISCU", "ISCU", "ISCU", "QTZL", "Qtzl", "Qtzl", "Qtzl", "Qtzl"))
d$se <- d$se <- as.numeric(as.character(d$se))

ggplot(d, aes(x = genotype, y = transcriptlevel, fill=genotype)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(
    type = "seq",
    palette = 3,
    direction = 1,
    aesthetics = "fill" ) + 
  geom_errorbar(aes(ymax = transcriptlevel + se, ymin= transcriptlevel - se), position = position_dodge())+
facet_wrap( ~ primer, ncol=2, scales="free")+
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#gleiches script mit einlesen 

d <-read.csv("daten_zsf.csv", header = TRUE, sep = ";")
d <- data.frame(d)
str(d)
d$transcriptlevel <- as.numeric(as.character(d$transcriptlevel))
d$genotype <- factor(d$genotype, level = c("wt", "OE", "Flag", "K", "F"))
d$primer <- factor(d$primer, level = c("ISCU", "ISCU", "ISCU", "ISCU", "ISCU", "QTZL", "Qtzl", "Qtzl", "Qtzl", "Qtzl"))
d$se <- d$se <- as.numeric(as.character(d$se))

cairo_pdf("qPCR_primer.pdf", 12, 7)
ggplot(d, aes(x = genotype, y = transcriptlevel, fill=genotype)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(
    type = "seq",
    palette = 3,
    direction = 1,
    aesthetics = "fill" ) + 
  geom_errorbar(aes(ymax = transcriptlevel + se, ymin= transcriptlevel - se), position = position_dodge())+
  facet_wrap( ~ primer, ncol=2, scales="free") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

#relvel-function zum sortieren
d <-read.csv("daten_zsf.csv", header = TRUE, sep = ";")
d <- data.frame(d)
data %>% 
  dplyr::mutate(genotype = fct_relevel(genotype, "wt", "OE", "Flag", "K", "F")) %>%

cairo_pdf("qPCR_primer.pdf", 12, 7)
ggplot(d, aes(x = genotype, y = transcriptlevel, fill=genotype)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(type = "seq", palette = 3, direction = 1, aesthetics = "fill" ) + 
  geom_errorbar(aes(ymax = transcriptlevel + se, ymin= transcriptlevel - se), position = position_dodge())+
  facet_wrap( ~ primer, ncol=2, scales="free") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()


# template
ddf <- data.frame(
  time  = c("t1", "t2", "t1", "t2"),
  group = c("f", "f", "m", "m"),
  mean  = c(1.6, 1.8, 3.5, 4.1),
  se    = c(0.1, 0.3, 0.3, 0.2)
)
# package ggplot2 has to be loaded
require("ggplot2")

# create a plot object, define dataframe to use, add gender to differ in colour already in base plot
pplot <- ggplot(ddf, x=group, y=mean, aes(group, mean, fill = time))

# create a bar plot
pplot +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.7), width=0.2) +
  #geom_point(stat="identity", shape=21, size=5, position=position_dodge(width=0.7), width=0.2)
  geom_point(stat="identity", shape=21, size=5, position=position_dodge(width=0.7))


#according to template
ddf <- data.frame(
  time = c("NFS1", "NFS1", "NFS1", "NFS1", "NFS1", "ISCU", "ISCU", "ISCU", "ISCU", "ISCU"),
  group = c("wt", "OE", "Flag", "K", "F", "wt", "OE", "Flag", "K", "F"),
  mean  = c(1.000, 115.162, 71.477, 55.136, 66.766, 1.000, 0.865, 0.838, 0.956, 1.350),
  se    = c(0.084344263, 2.142360052, 13.36486933, 2.244103961, 1.595136296, 0.17972206, 0.046230542, 0,387379952
            0.183992176, 0.133365007)
  # package ggplot2 has to be loaded
  # create a plot object, define dataframe to use, add gender to differ in colour already in base plot
  pplot <- ggplot(ddf, x=group, y=mean, aes(group, mean, fill = time)
                  
                  # create a bar plot
                  pplot +
                    geom_bar(stat="identity", position=position_dodge()) +
                    geom_errorbar(aes(ymax = mean + se, ymin= mean - se), position = position_dodge(width=0.7), width=0.2) +
                    #geom_point(stat="identity", shape=21, size=5, position=position_dodge(width=0.7), width=0.2)
                    geom_point(stat="identity", shape=21, size=5, position=position_dodge(width=0.7))
                  
                  

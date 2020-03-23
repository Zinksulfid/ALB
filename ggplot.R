install.packages("ggplot2")
library("ggplot2")
# erster Versuch 
a <- c("wt", "OE", "Flag", "K", "F")
b <- c(1.000, 115.162, 71.477, 55.136, 66.766)
d <- rbind(a,b)
colnames(d) <- a
rownames(d) <- c( "genotype", "transcriptlevel")
d <- t(d)  ## Die Tabelle hat die verkehrte Orientierung für ggplot, aes bezieht sich auf Spalten, hier hast du aber Zeilen. Deswegen transformiere ich mit t(d), das ändert die Orientierung
d <- data.frame(d)

ggplot(d, aes(x = genotype, y = transcriptlevel)) +
  geom_bar(stat="identity", position=position_dodge())

# (1) Ich hab rausgenommen, dass der ggplot in pp gespeichert wird und man das separat wieder auslesen muss. Das kann man später machen, wenn der Code funktioniert, gerade ist das zu umständlich
# (2) Wie du siehst, sind die Genotypen nicht in der richtigen Reihenfolge. Das kann man theoretisch mit as.vector(d$genotype, levels = d$genotype) lösen
# (3) Auch die Zahlen sind nicht fließend abgebildet, sondern scheinbar als Faktor eingelesen. Das kannst du überprüfen, indem du str(d) laufen lässt. Aha, d$transcriptlevel ist tatsächlich ein Faktor. Mach daraus eine Zahlenreihe mit as.numeric(as.character(x))





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
                  
                  
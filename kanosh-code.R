#clean routine
rm (list=ls())

#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("gridExtra")
#install.packages("grid")

library(vegan) 
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(tidyverse)

###Load the Data####
setwd("//naturkundemuseum-berlin.de/MuseumDFSRoot/Benutzer/richard.hofmann/Desktop/kanosh")
#setwd("//yourlocalpath")

#Data <- read.csv(url("https://raw.githubusercontent.com/fossilrich/Kanosh/master/kanosh.csv"),
#                 header = TRUE,  sep=";", check.names=FALSE)

Data <-  read.csv("//naturkundemuseum-berlin.de/MuseumDFSRoot/Benutzer/richard.hofmann/Desktop/kanosh/kanosh.csv", 
                     header = TRUE, sep=";", check.names=FALSE)


####Multivariate Analysis####

#Trilo's and Brachs only and Collection n>10
#data.vetted <- read.csv("//naturkundemuseum-berlin.de/MuseumDFSRoot/Benutzer/richard.hofmann/Desktop/Kanosh/kanosh-vet.csv", header = TRUE, row.names = 1, sep=";", check.names=FALSE)
data.vetted <- read.csv("//naturkundemuseum-berlin.de/MuseumDFSRoot/Benutzer/richard.hofmann/Desktop/kanosh/kanosh-vet.csv", header = TRUE, row.names = 1, sep=";", check.names=FALSE)


kanosh.multi <- data.vetted [,-c(1:6)] 


#kanosh.multi.brach <- kanosh.multi[,-c(1:6)] 
#kanosh.multi.tril <- - kanosh.multi[,-c(7:11)] 


####Base Figure for Cluster####
#dat.vet <- read.csv("//naturkundemuseum-berlin.de/MuseumDFSRoot/Benutzer/richard.hofmann/Desktop/Kanosh/kanosh-vet.csv", header = TRUE, sep=";", check.names=FALSE)
dat.vet <- read.csv("kanosh-vet.csv", header = TRUE, sep=";", check.names=FALSE)
data.fig <- melt(dat.vet, id.vars = c(1:7)) 
colnames(data.fig)[8] <- "species"
colnames(data.fig)[9] <- "abundance"
data.fig[data.fig == 0] <- NA

ggplot(data.fig, aes(y=Sample, x=species)) +
  geom_point(aes(size = abundance))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "black"),
        # panel.background = element_rect(fill = "white", colour = "black"),
        #panel.grid.major.x = element_blank(),
        #panel.grid.minor.x = element_blank(),
        #   axis.title.y = element_text(size=10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(size = 9,angle = 45, hjust = 1),
        #   axis.text.y = element_text(size = 9),
        #axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank())
#ggsave(("occurences.eps"), plot = last_plot(), width=12, height=19, units = "cm")


####Clusteranalysis####

distances.kanosh.raup <- vegdist(kanosh.multi, method = "raup")
#try different method =""
#"manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", 
#"jaccard", "gower", "altGower", "morisita", "horn", "mountford",
#"raup", "binomial", "chao", "cao" or "mahalanobis".
kanosh.clust.raup <- hclust(distances.kanosh.raup, method = "average")
#try also "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
plot (kanosh.clust.raup)#raup appears to yield well resolved cluster



cluster.raup <- as.data.frame (cutree(kanosh.clust.raup, h = 0.2)) #choose "h" in accordance with desired number of clusters or, respectively, at the distance level which should serve as threshold to delineate cluster 
#cluster.bray <- as.data.frame (cutree(kanosh.clust.bray, k = 3)) # "k" indicates the number of desired cluster
#cluster.jacc <- as.data.frame (cutree(kanosh.clust.jacc, k = 3)) # "k" indicates the number of desired cluster

names(cluster.raup)[1]<-"groups"
#names(cluster.bray)[1]<-"groups"
#names(cluster.jacc)[1]<-"groups"

####NDMS####
mdsKanosh <- metaMDS(kanosh.multi, distance = "jaccard", k = 3, trymax = 20, plot = TRUE)
stressplot(mdsKanosh)
ordiplot (mdsKanosh, display = "sites", type = 't')
#orditorp(mdsKanosh, display="sites")
ordihull(mdsKanosh,groups=cluster.raup$groups, draw="polygon",col="grey96",label=T)

####Built the dataframe for associations 
kanosh.multi.cluster <- cbind(kanosh.multi, cluster.raup)
kanosh.multi.ass <- aggregate(. ~  groups, data = kanosh.multi.cluster, sum)
kanosh.multi.ass.melt <- melt(kanosh.multi.ass, id.vars = 1)
colnames(kanosh.multi.ass.melt)[2] <- "species"
colnames(kanosh.multi.ass.melt)[3] <- "abundance"
is.na(kanosh.multi.ass.melt) <- !kanosh.multi.ass.melt

association1 <- subset(kanosh.multi.ass.melt, groups == 1)
association1 <- na.omit(association1)
association1 <- droplevels(association1)
sum(association1$abundance)

association2 <- subset(kanosh.multi.ass.melt, groups == 2)
association2 <- na.omit(association2)
association2 <- droplevels(association2)
sum(association2$abundance)

association3 <- subset(kanosh.multi.ass.melt, groups == 3)
association3 <- na.omit(association3)
association3 <- droplevels(association3)
sum(association3$abundance)

#Getting percentages
association1$percent <- association1$abundance/sum(association1$abundance) * 100
association2$percent <- association2$abundance/sum(association2$abundance) * 100
association3$percent <- association3$abundance/sum(association3$abundance) * 100

####associations####
#plot.1 <- 
ggplot(data=association1,  aes(x= reorder(species, -percent), y=percent))+
  geom_bar(stat = "identity")+
  labs(    x = "species",    y = "abundance [%]", title = "A")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "black"), 
        # panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        #   axis.title.y = element_text(size=10),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        #axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank())
#ggsave(("assoc1.eps"), plot = last_plot(), width=8, height=10, units = "cm")

#plot.2  <- 
ggplot(data=association2,  aes(x= reorder(species, -percent), y=percent))+
  geom_bar(stat = "identity")+
  labs(    x = "species",    y = "abundance [%]", title = "B")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "black"), 
        # panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_text(size=10),
        axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        # axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.3),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        #axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank())
#ggsave(("assoc2.eps"), plot = last_plot(), width=8, height=10, units = "cm")
#plot.3  <- 

ggplot(data=association3, aes(x= reorder(species, -percent), y=percent))+
  geom_bar(stat = "identity")+
  labs(    x = "species",    y = "abundance [%]", title = "C")+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "black"), 
        # panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        #   axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10), 
        # axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.3),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        #axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank())
#ggsave(("assoc3.eps"), plot = last_plot(), width=8, height=10, units = "cm")



####rarefaction####
rare.ass1 <- cbind.data.frame(association1$species, association1$abundance)
colnames(rare.ass1) <- c("species","association1")
rare.ass2 <- cbind.data.frame(association2$species, association2$abundance)
colnames(rare.ass2) <- c("species","association2")
rare.ass3 <- cbind.data.frame(association3$species, association3$abundance)
colnames(rare.ass3) <- c("species","association3")

rarefac <- merge(rare.ass1, rare.ass2, all = TRUE)
rarefac <- merge(rarefac, rare.ass3, all = TRUE)
rarefac[is.na(rarefac)] <- 0

rownames(rarefac) <- rarefac$species
rarefac <- rarefac[,-1]

#how to transpose a dataframe in tidyverse
test <- rarefac %>% 
  rownames_to_column %>%
  gather(variable, value, -rowname) %>% 
  spread(rowname, value)

rownames(test) <- test$variable
test <- test[,-1]

S <- specnumber(test)
Srare <- rarefy(test, se = TRUE)
plot(Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(test, step = 4,  col = "blue", cex = 0.6)


#####data preparation for stratigraphic distribution####
Samplenames <- rownames(cluster.raup)
#rownames(cluster.raup) <- NULL
cluster.raup2 <- cbind(Samplenames,cluster.raup)

Data2 <- merge(cluster.raup2, Data, by.x = 1, by.y = 1, all.y = TRUE)


#take the groups from the clusters 
data.melt <- melt(Data2, id.vars = c(1:8)) 
colnames(data.melt)[3] <- "composite"
colnames(data.melt)[9] <- "species"
colnames(data.melt)[10] <- "abundance"
data.melt[data.melt == 0] <- NA

#add a new column where we copy the species'
data.melt$group <- data.melt$species

#now major groups will be added fpr later use as a factor based on species name
levels(data.melt$group)[levels(data.melt$group)%in%c("cf. Bathyurus sp.","Punka cf. nitida","Pseudoolenoides dilectus",
                                                     "Pseudoolenoides pogonipensis","Pseudoolenoides ludificatus","Pseudomera barrandei","Kanoshia sp.","Bathyurellus pogonipensis","Kanoshia kanoshensis")] <- "Trilobite"
levels(data.melt$group)[levels(data.melt$group)%in%c("Shoshonorthis michaelis","Anomalorthis  lonensis","Anomalorthis utahensis",
                                                     "Desmorthis nevadensis","Hesperonomiella minor","cf. Wahwahlingula sp.")] <- "Brachiopod"
levels(data.melt$group)[levels(data.melt$group)%in%c("Murchisonia sp.", "Lophospira perangulata","Malayispira sp.")] <- "Gastropod"
levels(data.melt$group)[levels(data.melt$group)%in%c("Receptaculites","leperditicopids")] <- "other"
data.melt$groups <- as.factor(data.melt$groups)

####Plot Species (abundance) alongside Sections####
ggplot(data.melt, aes(y=composite, x=species)) +
  geom_point(aes(size = abundance, colour = groups))+
  #facet_wrap( ~ Section, ncol =3,  scales = "fixed", drop=TRUE)+
 #facet_grid(cols = vars(Section), space = "free_x")+
  #find a way to get rid of unused levels within facet,  ->drop = TRUE in facet_wrap doesnt work dumbass
   geom_hline(yintercept=c(33, 47, 90, 117, 140), colour= "black")+
  #geom_hline(yintercept= c(2.6, 9.8,11.2,15.9,24.5,28,35, 35.5,36,38.9, 41.5, 
  #                         44, 46,48,50,51,69,70.5,73,76,80,85,86, 90.2,104.4,
  #                         109,119, 121.5, 126, 135, 161,167.8,188.4), colour= "grey", linetype = "dashed")+
  #geom_hline(yintercept= lines, colour= "pink")+
  scale_x_discrete(position = "top") +
  scale_y_continuous("Stratigraphic Unit", sec.axis = sec_axis(~. *1), 
                     breaks = c(15,40,68.5, 104.5,128.5, 165), 
                     labels = c("Lower Shale","Silty Limestone","Upper Shale", 
                                "Sandstone", "Calcisiltite", "Lehman Formation"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(colour = "black", fill = "lightgrey"),
        panel.grid.major = element_line(colour = "black"),
        # panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #   axis.title.y = element_text(size=10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(size = 9,angle = 55, hjust = 0, vjust = 0.5),
        #   axis.text.y = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.ticks.y.right = element_blank())
#ggsave(("fossilranges.pdf"), plot = last_plot(), width=15, height=25, units = "cm")
#ggsave(("fossilranges.png"), plot = last_plot(), width=15, height=27, units = "cm", dpi = 1000)


#####Stratigraphic Distribution####
 
ggplot(data.melt, aes(y=composite, x=species)) +
  geom_point(aes(size = abundance, colour = groups))+
  facet_wrap( ~ Section, ncol =3,  scales = "fixed", drop=TRUE)+
  #facet_grid(cols = vars(Section), space = "free_x")+
  #find a way to get rid of unused levels within facet,  ->drop = TRUE in facet_wrap doesn't work dumbass
  # geom_hline(yintercept=c(33, 47, 90, 117, 140), colour= "black")+
  geom_hline(yintercept= c(2.6, 9.8,11.2,15.9,24.5,28,35, 35.5,36,38.9, 41.5, 
                           44, 46,48,50,51,69,70.5,73,76,80,85,86, 90.2,104.4,
                           109,119, 121.5, 126, 135, 161,167.8,188.4), colour= "grey", linetype = "dashed")+
  #geom_hline(yintercept= lines, colour= "pink")+
  scale_x_discrete(position = "top") +
  scale_y_continuous("Stratigraphic Unit", sec.axis = sec_axis(~. *1), 
                     breaks = c(15,40,68.5, 104.5,128.5, 165), 
                     labels = c("Lower Shale","Silty Limestone","Upper Shale", 
                                "Sandstone", "Calcisiltite", "Lehman Formation"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(colour = "black", fill = "lightgrey"),
        panel.grid.major = element_line(colour = "black"),
        # panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        #   axis.title.y = element_text(size=10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(size = 9,angle = 55, hjust = 0, vjust = 0.5),
        #   axis.text.y = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.ticks.y.right = element_blank())
#ggsave(("fossilranges_wrap.pdf"), plot = last_plot(), width=20, height=27, units = "cm")
#ggsave(("fossilranges_wrap.eps"), plot = last_plot(), width=20, height=27, units = "cm")























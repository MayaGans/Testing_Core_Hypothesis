require(gplots)
require(tidyverse)
require(matrixStats)
require(factoextra)
require(gridExtra)
require(mctoolsr)
require(parallelDist)
require(hexbin)
require(qpcR)
require(cowplot)

# import OTU tables: OTUs as rows, sample replicates as columns
human <- readRDS("human.RDS")
arabadopsis <- readRDS("arabadopsis.RDS")

# by reassigning the imported dataset to the variable OTU
# all code below can be easily executed
otu <- human
print(otu[1:4,1:4])


# Data Clean Up
# Remove all non-numeric data and make OTU IDs into row names
otu1 <- otu[,-1]
rownames(otu1) <- otu[,1]
otu <- otu1

# Get number of sites for the experiment
# both datasets here are single treatments
# this will need to be modified if the data includes multiple treatments
N <- ncol(otu)

# Get the total number of reads for the entire experiment
# this number would also need to be the total number of reads per treatment
s <- sum(colSums(otu))

# visualize the distribution of reads
All_OTUs <- rownames(otu)


# Mean, Variance, Covariance

# In this section we calculate the mean, variance and CV of each OTU.  
# This will be used downstream to see if there is correlation between CV/Mean/Variance and the 
# OTUs defined as core.  The end product of this section is a dataframe with the Mean,
# Variance and CV for each OTU. 

#transpose matrix so we can use the two functions colMeans and colVars
dat<-as.matrix(t(otu))
OTU_M_V_CV<-as.data.frame(colMeans(dat))
OTU_M_V_CV$Mean<-OTU_M_V_CV$`colMeans(dat)`
OTU_M_V_CV$`colMeans(dat)`= NULL
OTU_M_V_CV$Variance<-colVars(dat)
OTU_M_V_CV$CV <- (OTU_M_V_CV$Variance/OTU_M_V_CV$Mean)
OTU_M_V_CV$OTU<-rownames(OTU_M_V_CV)
OTU_M_V_CV.2<-OTU_M_V_CV
OTU_M_V_CV.2$OTU=NULL

# Proportion of replicates method

# This method assigns taxa to the core based upon the number of sites it is present in. 
# In this example we assign core membership when the taxa's abundance is atleast 10 fold 
# the number of sites. This method accounts for abundance as a function of sites. 

#taxa rows and sites columns
focaltaxa<-otu
focaltaxa$numofsites <- apply(focaltaxa, 1, function(x) sum(x>0))
focaltaxa<-as.data.frame(subset(focaltaxa, focaltaxa$numofsites >=N*0.5))
focaltaxa<-as.character(rownames(focaltaxa))



prop$numofsites <- apply(prop, 1, function(x) sum(x>0))
prop<-as.data.frame(subset(prop, prop$numofsites >=N*0.5))

# Hard cut offs

# This method assigns taxa to the core if they are present in more than a pre-determined number of 
# sites and have a total abudance greater than a pre-determined number of reads. In our example we
# set the minimum number of sites to 5 and the minimum number of reads to 25. Here we use the hard 
# cut off described in Lundberg (2012), but realize this is any threshhold. 

M5_25<- otu
#M5_25 <- as.data.frame(subset(df5_25, (rowSums(df5_25) >= 25)))
#from maya
#sum(rownsums(replicatesSim>=25))>=5)- ncore
M5_25$numofsites <- apply(M5_25, 1, function(x) sum(x>25))
M5_25<-as.data.frame(subset(M5_25, M5_25$numofsites >=5))
M5_25<-as.character(rownames(M5_25))
length(M5_25)

# Proportion of reads and replicates

# This method assigns taxa to the core if they account for some proportion of the total
# reads for the sequencing run and if they are preseant in atleast x% of the total number of 
# replicates. In this example, a core taxa must account for 0.01% of the total reads for the entire 
# otu table and be present in at least 50% of sites.

prop <- as.data.frame(subset(otu, (rowSums(otu) >= (s/1000))))
prop$numofsites <- apply(prop, 1, function(x) sum(x>0))
prop<-as.data.frame(subset(prop, prop$numofsites >=N*0.5))
prop<- as.character(rownames(prop))

# Proportion of reads
# This method assigns taxa to the core if they are in the top X% of reads. 
# Taxa are ranked in abudnace and the cumulative sum is recoreded. 
# Any taxa which appears before some cutoff percentage is included in the core. 
# In this example, a taxa will be assigned to the core if they account for the first 75% of the reads

top_75_percent<-as.data.frame(otu)
top_75_percent$otuappearance<- rowSums(top_75_percent)
sortedtop_75_percent<-top_75_percent[order(-top_75_percent$otuappearance),]
sortedtop_75_percent$prop<-sortedtop_75_percent$otuappearance/s
sortedtop_75_percent$cum_sum<-cumsum(sortedtop_75_percent$prop)
sortedtop_75_percent<-as.data.frame(subset(sortedtop_75_percent, sortedtop_75_percent$cum_sum <=0.75))
sortedtop_75_percent<-as.character(rownames(sortedtop_75_percent))

# Combine Methods

# Make a dataframe with all observed taxa, their inclusion to the core by method 
# (deliniated as a 1 or 0), the mean, variance, and coefficient of variation.

combined_data  <- qpcR:::cbind.na(focaltaxa, prop, sortedtop_75_percent, M5_25, All_OTUs )
combined_data <- as.data.frame(combined_data)
rownames(combined_data) <- combined_data$All_OTUs

combined_data2  <- qpcR:::cbind.na(focaltaxa, prop, sortedtop_75_percent, M5_25 )
combined_data2 <- as.data.frame(combined_data2)

Core_community_Ids<-unique(unlist(combined_data[,]))
Core_community_Ids<-as.character(Core_community_Ids)

combined_data1<- data.frame(OTU = Core_community_Ids, 
focaltaxa = Core_community_Ids %in% combined_data2$focaltaxa, 
prop = Core_community_Ids %in% combined_data2$prop,  
sortedtop_75_percent = Core_community_Ids %in% combined_data2$sortedtop_75_percent, 
M5_25 = Core_community_Ids %in% combined_data2$M5_25)


#convert true false into 1 and 0
index <- c("TRUE", "FALSE")
values <- c("1", "0")

combined_data1$focaltaxa<- values[match(combined_data1$focaltaxa, index)]
combined_data1$prop<- values[match(combined_data1$prop, index)]
combined_data1$M5_25<- values[match(combined_data1$M5_25, index)]
combined_data1$sortedtop_75_percent<- values[match(combined_data1$sortedtop_75_percent, index)]

MergedOTU_Stats_human <- na.omit(left_join(combined_data1, OTU_M_V_CV))


#################
################

test <- MergedOTU_Stats_human %>%
  gather(method, TF, 2:5) %>%
  mutate(method = factor(method, levels = c("focaltaxa",
                                            "prop",
                                            "sortedtop_75_percent",
                                            "M5_25"), 
                         labels = c("Proportion of Replicates",
                                    "Proportion of Sequence Reads and Replicates", 
                                    "Proportion of Sequence Reads", 
                                    "Hard Cutoff")))


test$method <- factor(test$method, 
                           levels = c("Proportion of Sequence Reads", 
                                      "Proportion of Replicates", 
                                      "Proportion of Sequence Reads and Replicates", 
                                      "Hard Cutoff"))


p1 <- ggplot(test, aes(x = log(Mean), y = CV, colour = TF)) +
  geom_hex(bins = 30) +
  facet_wrap(.~method, nrow = 1) +
  scale_fill_gradient(low = "#EAE6f3", high = "#432976") +
  scale_color_manual(values = c("lightgray", "black")) +
  theme_minimal() +
  guides(color = FALSE) +
  labs(fill = "HMP \n  Taxa Count") +
  theme(text = element_text(size = 20)) +
  ylab("Coefficent of Variance") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


##########################################
# Arabidopsis
#########################################

otu <- arabadopsis
print(otu[1:4,1:4])

otu1 <- otu[,-1]
rownames(otu1) <- otu[,1]
otu <- otu1

N <- ncol(otu)
s <- sum(colSums(otu))
All_OTUs <- rownames(otu)

#transpose matrix so we can use the two functions colMeans and colVars
dat<-as.matrix(t(otu))
OTU_M_V_CV<-as.data.frame(colMeans(dat))
OTU_M_V_CV$Mean<-OTU_M_V_CV$`colMeans(dat)`
OTU_M_V_CV$`colMeans(dat)`= NULL
OTU_M_V_CV$Variance<-colVars(dat)
OTU_M_V_CV$CV <- (OTU_M_V_CV$Variance/OTU_M_V_CV$Mean)
OTU_M_V_CV$OTU<-rownames(OTU_M_V_CV)
OTU_M_V_CV.2<-OTU_M_V_CV
OTU_M_V_CV.2$OTU=NULL

# Proportion of replicates method
#taxa rows and sites columns
focaltaxa<-otu
focaltaxa$numofsites <- apply(focaltaxa, 1, function(x) sum(x>0))
focaltaxa<-as.data.frame(subset(focaltaxa, focaltaxa$numofsites >=N*0.5))
focaltaxa<-as.character(rownames(focaltaxa))

# Hard cut offs
M5_25<- otu
#M5_25 <- as.data.frame(subset(df5_25, (rowSums(df5_25) >= 25)))
#from maya
#sum(rownsums(replicatesSim>=25))>=5)- ncore
M5_25$numofsites <- apply(M5_25, 1, function(x) sum(x>25))
M5_25<-as.data.frame(subset(M5_25, M5_25$numofsites >=5))
M5_25<-as.character(rownames(M5_25))
length(M5_25)

# Proportion of reads and replicates
prop <- as.data.frame(subset(otu, (rowSums(otu) >= (s/1000))))
prop$numofsites <- apply(prop, 1, function(x) sum(x>0))
prop<-as.data.frame(subset(prop, prop$numofsites >=N*0.5))
prop<- as.character(rownames(prop))

# Proportion of reads
top_75_percent<-as.data.frame(otu)
top_75_percent$otuappearance<- rowSums(top_75_percent)
sortedtop_75_percent<-top_75_percent[order(-top_75_percent$otuappearance),]
sortedtop_75_percent$prop<-sortedtop_75_percent$otuappearance/s
sortedtop_75_percent$cum_sum<-cumsum(sortedtop_75_percent$prop)
sortedtop_75_percent<-as.data.frame(subset(sortedtop_75_percent, sortedtop_75_percent$cum_sum <=0.75))
sortedtop_75_percent<-as.character(rownames(sortedtop_75_percent))

# Combine Methods
combined_data  <- qpcR:::cbind.na(focaltaxa, prop, sortedtop_75_percent, M5_25, All_OTUs )
combined_data <- as.data.frame(combined_data)
rownames(combined_data) <- combined_data$All_OTUs
combined_data2  <- qpcR:::cbind.na(focaltaxa, prop, sortedtop_75_percent, M5_25 )
combined_data2 <- as.data.frame(combined_data2)
Core_community_Ids<-unique(unlist(combined_data[,]))
Core_community_Ids<-as.character(Core_community_Ids)
combined_data1<- data.frame(OTU = Core_community_Ids, 
                            focaltaxa = Core_community_Ids %in% combined_data2$focaltaxa, 
                            prop = Core_community_Ids %in% combined_data2$prop,  
                            sortedtop_75_percent = Core_community_Ids %in% combined_data2$sortedtop_75_percent, 
                            M5_25 = Core_community_Ids %in% combined_data2$M5_25)
index <- c("TRUE", "FALSE")
values <- c("1", "0")
combined_data1$focaltaxa<- values[match(combined_data1$focaltaxa, index)]
combined_data1$prop<- values[match(combined_data1$prop, index)]
combined_data1$M5_25<- values[match(combined_data1$M5_25, index)]
combined_data1$sortedtop_75_percent<- values[match(combined_data1$sortedtop_75_percent, index)]

MergedOTU_Stats_arab <- na.omit(left_join(combined_data1, OTU_M_V_CV))
#focal taxa  Proportion of Replicates
table(MergedOTU_Stats_arab[,2])
#prop Proportion of Sequence Reads and Replicates
table(MergedOTU_Stats_arab[,3])
#sorted top 75 prop reads Proportion of Sequence Reads
table(MergedOTU_Stats_arab[,4])
# m5_25 Hard Cutoff
table(MergedOTU_Stats_arab[,5])


#focal taxa  Proportion of Replicates
table(MergedOTU_Stats_human[,2])
#prop Proportion of Sequence Reads and Replicates
table(MergedOTU_Stats_human[,3])
#sorted top 75 prop reads Proportion of Sequence Reads
table(MergedOTU_Stats_human[,4])
# m5_25 Hard Cutoff
table(MergedOTU_Stats_human[,5])

################
# Plot
################

test <- MergedOTU_Stats_arab %>%
  gather(method, TF, 2:5) %>%
  mutate(method = factor(method, levels = c("focaltaxa",
                                            "prop",
                                            "sortedtop_75_percent",
                                            "M5_25"), 
                         labels = c("Proportion of Replicates",
                                    "Proportion of Sequence Reads and Replicates", 
                                    "Proportion of Sequence Reads", 
                                    "Hard Cutoff")))

test$method <- factor(test$method, 
                      levels = c("Proportion of Sequence Reads", 
                                 "Proportion of Replicates", 
                                 "Proportion of Sequence Reads and Replicates", 
                                 "Hard Cutoff"))


p2 <- ggplot(test, aes(x = log(Mean), y = CV, colour = TF)) +
  geom_hex(bins = 30) +
  facet_wrap(.~method, nrow = 1) +
  scale_fill_gradient(low = "#E6ECF1", high = "#284B5E") +
  scale_color_manual(values = c("lightgray", "black")) +
  theme_minimal() +
  guides(color = FALSE) +
  labs(fill = "Arabidopsis \n Taxa Count") +
  theme(text = element_text(size = 20)) +
  ylab("Coefficent of Variance") +
  theme(strip.text = element_text(color="white"))


# g <- gridExtra::arrangeGrob(p1, p2, nrow = 2)
g <- cowplot::plot_grid(p1, p2, labels=c("A", "B"), ncol = 1, nrow = 2)
ggsave("/Users/mayagans/Desktop/Core_Hypothesis/Figure3.pdf", g, width = 25, height = 12)
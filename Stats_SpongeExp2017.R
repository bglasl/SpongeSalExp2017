####### Statistical Analysis ###########
# Glasl et al. (2018) Exploring the diversity-stability paradigm using sponge microbial communities. Scientific Reports 8:8425
# doi: DOI:10.1038/s41598-018-26641-9

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

#### load PHYLOSEQ_TABLE ####
load(file="PHYLOSEQ_TABLE.RData")
METADATA<-read.csv("Metadata_2.csv", header = TRUE, sep = ",")

##########################
##### ALPHA DIVERSITY ####
##########################
RICHNESS<-estimate_richness(PHYLOSEQ_TABLE, split=TRUE, measures=c("Shannon","Observed"))
RICHNESS<-as.data.frame(RICHNESS)
RICHNESS<-tibble::rownames_to_column(as.data.frame(RICHNESS), var="Sample_ID")
Richness_table<-right_join(METADATA, RICHNESS, by="Sample_ID")

MEAN_Richness<-Richness_table %>% dplyr::group_by(Code_short) %>% dplyr:: summarise(mean(Observed),sd(Observed))
write.csv(MEAN_Richness, file="Richness.csv")

MEAN_Richness<-Richness_table %>% dplyr::group_by(Host) %>% dplyr:: summarise(mean(Observed),sd(Observed))
write.csv(MEAN_Richness, file="Richness_Host.csv")

MEAN_Shannon<-Richness_table %>% dplyr::group_by(Code_short) %>%dplyr:: summarise(mean(Shannon),sd(Shannon))
write.csv(MEAN_Shannon, file="Shannon.csv")

MEAN_Shannon<-Richness_table %>% dplyr::group_by(Host_Treatment) %>% dplyr::summarise(mean(Shannon),sd(Shannon))
write.csv(MEAN_Shannon, file="Shannon_Treatment.csv")

MEAN_Shannon<-Richness_table %>% dplyr::group_by(Host) %>% dplyr:: summarise(mean(Shannon),sd(Shannon))
write.csv(MEAN_Shannon, file="Shannon_Host.csv")

library(vegan)
Evenness<-as.data.frame((Richness_table$Shannon)/(log(Richness_table$Observed)))
Evenness_table<-bind_cols(Richness_table, Evenness)
Evenness_table<-Evenness_table %>% dplyr::rename(Evenness=`(Richness_table$Shannon)/(log(Richness_table$Observed))`)
MEAN_Evenness<-Evenness_table %>% dplyr::group_by(Code_short) %>% dplyr::summarise(mean(Evenness),sd(Evenness))
write.csv(MEAN_Evenness, file="Evenness.csv")

#### plot Alpha Div ####
library(nlme)
library(car)

Richness_table$SamplingTimepoint<-factor(Richness_table$SamplingTimepoint, levels=c("before disturbance", "24h after disturbance", "168h after disturbance"))

# function for Standard Error
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

#calculate SE
Richness_summary<-summarySE(Richness_table, measurevar = "Shannon", groupvars = c("Host","SamplingTimepoint","Treatment"))
Richness_summary

#plot Shannon index for all Hosts
pd<-position_dodge(0.3)

ShannonPlot<-ggplot(Richness_summary, aes(x=SamplingTimepoint, y=Shannon, colour=Treatment, group=Treatment))+
  facet_grid(.~Host)+
  geom_errorbar(aes(ymin=Shannon-se, ymax=Shannon+se), width=.1, position = pd)+
  geom_line(position = pd)+
  geom_point(position = pd, size=2)+
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10),
        panel.grid =element_blank())+
  scale_x_discrete(breaks=c("before disturbance","24h after disturbance","168h after disturbance"),
                   labels=c("1", "11", "17"), name="sampling day")+
  scale_color_manual(breaks=c("Control","Pulse"),
                     labels=c("Control", "Disturbance"), values=c(Control="#7fcdbb",Pulse="#2c7fb8"))

ShannonPlot

pdf('Shannon_Graph.pdf', width=6, height=6)
print(ShannonPlot)
graphics.off()

setEPS()
postscript("Shannon_Graph.eps", width=6, height=6)
ShannonPlot
dev.off()

### ANOVA AlphaDiv ####
# ANOVA with interactions of SamplingTimepoint*Treatment*Host
model1<-aov(Shannon~SamplingTimepoint*Treatment*Host, data=Richness_table)
summary(model1) # Host sig ***
capture.output(summary(model1),file="ANOVA_Shannon.doc")

model1<-aov(Shannon~Host, data=Richness_table)
summary(model1)
TukeyHSD(model1)
capture.output(TukeyHSD(model1),file="TukeyANOVA_Shannon.doc")

# repeated measures ANOVA - check if Genotype randome effects has an influence on the dataset
library(lme4)
model2<-lmer(Shannon~SamplingTimepoint*Treatment*Host+(1|Genotype), data=Richness_table)
summary(model2) # Genotype random effects explains for 0.04280 of the Variance in the dataset - this is very low and therefore Genotype can be ignored

library(nlme)
model3<-lme(Shannon~SamplingTimepoint*Treatment*Host, random =~1|Genotype, data=Richness_table) #same as model2 but with different package
summary(model3)

#check which model is the best
library(MuMIn)
AICc(model1, model2, model3) #model one has the lowest AICc

#### plot Richness  ####
library(nlme)
library(car)

head(Richness_table)
Richness_table$SamplingTimepoint<-factor(Richness_table$SamplingTimepoint, levels=c("before disturbance", "24h after disturbance", "168h after disturbance"))

# function for Standard Error
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

#calculate SE
Richness_summary<-summarySE(Richness_table, measurevar = "Observed", groupvars = c("Host","SamplingTimepoint","Treatment"))
Richness_summary

#plot Observed index for all Hosts
pd<-position_dodge(0.3)

ObservedPlot<-ggplot(Richness_summary, aes(x=SamplingTimepoint, y=Observed, colour=Treatment, group=Treatment))+
  facet_grid(.~Host)+
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se), width=.1, position = pd)+
  geom_line(position = pd)+
  geom_point(position = pd, size=2)+
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10),
        panel.grid =element_blank())+
  scale_x_discrete(breaks=c("before disturbance","24h after disturbance","168h after disturbance"),
                   labels=c("1", "11", "17"), name="sampling day")+
  scale_color_manual(breaks=c("Control","Pulse"),
                     labels=c("Control", "Disturbance"), values=c(Control="#7fcdbb",Pulse="#2c7fb8"))

ObservedPlot

pdf('Richness_Graph.pdf', width=6, height=6)
print(ObservedPlot)
graphics.off()

model1<-aov(Observed~SamplingTimepoint*Treatment*Host, data=Richness_table)
summary(model1) 


##################
#### Beta Div ####
##################

# transform sample counts to relative abundance
PHYLOSEQ_TABLE = transform_sample_counts(PHYLOSEQ_TABLE, function(x)100*x/sum(x))
otu_table(PHYLOSEQ_TABLE)
PHYLOSEQ_TABLE

#### NMDS ####
set.seed(1000)
ORD<-ordinate(PHYLOSEQ_TABLE,"NMDS","bray", k=2, trymax=1000)
library(vegan)
library(data.table)
library(plyr)
test<-as.data.frame(scores(ORD))
test<-setDT(test, keep.rownames=TRUE)

detach("package:plyr", unload=TRUE)
test<-test %>% rename(Sample_ID=rn) #plyr masks dplyr::rename
ALLHosts<-METADATA %>% select(Host_Treatment, Sample_ID, SamplingTimepoint)

test<- inner_join(test,ALLHosts,by="Sample_ID")
test<- test%>%select(NMDS1,NMDS2,Host_Treatment,SamplingTimepoint)
head(test)

library(plyr)
find_hull <- function(test) test[chull(test$NMDS1, test$NMDS2), ]
hulls <- ddply(test, c("Host_Treatment"), find_hull)
hulls
detach("package:plyr", unload=TRUE)

plot<-plot_ordination(PHYLOSEQ_TABLE, ORD, type="samples",color = "Host_Treatment")
plot$layers=plot$layers[-1]
plot<-plot +
  geom_point(aes(shape=SamplingTimepoint), size=3)+ 
  theme_classic()+
  geom_polygon(data=hulls, alpha=0, show.legend = FALSE, linetype=1)+
  scale_colour_manual(values=c(AQ_Control="#73d09f",AQ_Treatment="#50916f",
                               CO_Control="#aaaaaa",CO_Treatment="#37677E",
                               CY_Control="#a95382",CY_Treatment="#45416b",
                               IB_Control="#f9e07c",IB_Treatment="#ddba2d",
                               IR_Control="#ae1c1c",IR_Treatment="#6b0202",
                               ST_Control="#fab74f",ST_Treatment="#ff8533"),
                      breaks=c("AQ_Control","AQ_Treatment",
                               "CO_Control","CO_Treatment",
                               "CY_Control","CY_Treatment",
                               "IB_Control","IB_Treatment",
                               "IR_Control","IR_Treatment",
                               "ST_Control","ST_Treatment"),
                      labels=c("AQ Control","AQ Disturbance",
                               "CO Control","CO Disturbance",
                               "CY Control","CY Disturbance",
                               "IB Control","IB Disturbance",
                               "IR Control","IR Disturbance",
                               "ST Control","ST Disturbance"))+
  annotate("text", x=-3, y=-1.7, label= "stress = 0.1717499", hjust=0)+
  annotate("text", x=-3, y=-1.9, label="k = 2", hjust=0)+
  guides(shape=guide_legend(title="sampling day"))+
  guides(color=guide_legend(title="host and treatments"))+
  scale_shape_discrete(breaks=c("before disturbance","24h after disturbance","168h after disturbance"),
                       labels=c("1", "11", "17"))

plot

pdf('NMDS_HostTreatment.pdf', width=8, height=6)
print(plot)
graphics.off()

setEPS()
postscript('NMDS_HostTreatment.eps', width=8, height=6)
plot
dev.off()

#### ANOSIM ####
# An R value close to "1.0" suggests dissimilarity between groups 
# while an R value close to "0" suggests an even distribution of high and low ranks within and between groups. 
# R values below "0" suggest that dissimilarities are greater within groups than between groups.

##all Hosts ##
Host_group<-get_variable(PHYLOSEQ_TABLE, "Host")
HOSTANOSIM<-anosim(distance(PHYLOSEQ_TABLE,"bray"), Host_group)
HOSTANOSIM$signif #p=0.001
HOSTANOSIM$statistic #R=0.979312 --> Each sponge host is associated with a "unique" microbiome

Treatment_group<-get_variable(PHYLOSEQ_TABLE, "Treatment")
HOSTANOSIM<-anosim(distance(PHYLOSEQ_TABLE,"bray"), Treatment_group, strata=Host_group)
HOSTANOSIM$signif # p=0.0.27
HOSTANOSIM$statistic # R=-0.0070019 --> even distribution between treatments within host species

Genotype_group<-get_variable(PHYLOSEQ_TABLE, "Genotype")
HOSTANOSIM<-anosim(distance(PHYLOSEQ_TABLE,"bray"), Genotype_group, strata=Host_group)
HOSTANOSIM$signif # p=0.001
HOSTANOSIM$statistic # R=0.9427445 

MicrobiomeDiversity_group<-get_variable(PHYLOSEQ_TABLE, "MicrobiomeDiversity")
HOSTANOSIM<-anosim(distance(PHYLOSEQ_TABLE,"bray"), MicrobiomeDiversity_group)
HOSTANOSIM$signif # p=0.001
HOSTANOSIM$statistic # p=0.4109492

Tank_group<-get_variable(PHYLOSEQ_TABLE, "TankNr")
HOSTANOSIM<-anosim(distance(PHYLOSEQ_TABLE,"bray"), Tank_group)
HOSTANOSIM$signif # no tank effect p=0.789

TreatmentComb_group<-get_variable(PHYLOSEQ_TABLE, "TreatmentComb")
HOSTANOSIM<-anosim(distance(PHYLOSEQ_TABLE,"bray"), TreatmentComb_group)
HOSTANOSIM$signif # not sign!

HT_group<-get_variable(PHYLOSEQ_TABLE, "Host_Treatment")
HOSTANOSIM<-anosim(distance(PHYLOSEQ_TABLE,"bray"), HT_group)
HOSTANOSIM$signif # sign!
HOSTANOSIM$statistic # p=0.901555 -> this result might be biased because Host has a significant effect 

Time_group<-get_variable(PHYLOSEQ_TABLE, "SamplingTimepoint")
HOSTANOSIM<-anosim(distance(PHYLOSEQ_TABLE,"bray"), Time_group, strata=HT_group)
HOSTANOSIM$signif #p=0.003
HOSTANOSIM$statistic #R=-0.01117194

#### ADONIS & BETA DISPERSION ####
# Adonis = permutational Multivariate Analysis of variance using distance matrices
### using adonis2 ### - this is the way to go!
df=as(sample_data(PHYLOSEQ_TABLE), 'data.frame')
df=df%>%dplyr::mutate(Genotype_SamplingTimepoint=paste(Genotype,SamplingTimepoint,sep="_"))
d=phyloseq::distance(PHYLOSEQ_TABLE,'bray')

perm<-how(nperm=10000)
setBlocks(perm)<-with(df, Host_Treatment)
ADONIS2<-adonis2(d~SamplingTimepoint, data=df, permutations = perm, method = "bray")
ADONIS2 #not sign -> no shift in the microbial community over time irrespective of treatment
capture.output(ADONIS2,file="Adonis2_TimeeffectHOSTTREATMENT.doc")

perm<-how(nperm=10000)
setBlocks(perm)<-with(df, Genotype)
ADONIS2<-adonis2(d~SamplingTimepoint, data=df, permutations = perm, method = "bray")
ADONIS2 #0.0104 *

perm<-how(nperm=10000)
setBlocks(perm)<-with(df, Host_Treatment)
ADONIS2<-adonis2(d~Genotype, data=df, permutations = perm, method = "bray")
ADONIS2 # 9.999e-05 ***

perm<-how(nperm=10000)
setBlocks(perm)<-with(df, Host)
ADONIS2<-adonis2(d~Genotype, data=df, permutations = perm, method = "bray")
ADONIS2 #p=9.999e-05 *** microbial communities of different genotypes of the same species differ in their composition
capture.output(ADONIS2,file="Adonis2_GenotypeeffectHOST.doc")

library(RVAideMemoire)
PAIRWISE<-pairwise.perm.manova(d,df$Host, nperm=10000)
PAIRWISE # sign. difference between all species 

dfC<-df%>%filter(Treatment=="Control")
PRUNED<-phyloseq::subset_samples(PHYLOSEQ_TABLE, Treatment=="Control")
d=phyloseq::distance(PRUNED,'bray')
perm<-how(nperm=10000)
setBlocks(perm)<-with(dfC, Genotype)
ADONIS2<-adonis2(d~SamplingTimepoint, data=dfC, permutations = perm, method = "bray")
ADONIS2 #0.0019 *

dfC<-df%>%filter(Treatment=="Pulse")
PRUNED<-phyloseq::subset_samples(PHYLOSEQ_TABLE, Treatment=="Pulse")
d=phyloseq::distance(PRUNED,'bray')
perm<-how(nperm=10000)
setBlocks(perm)<-with(dfC, Genotype)
ADONIS2<-adonis2(d~SamplingTimepoint, data=dfC, permutations = perm, method = "bray")
ADONIS2 #0.2724

# effect of genotype per host species
#AQ
dfC<-df%>%filter(Host=="AQ")
PRUNED<-phyloseq::subset_samples(PHYLOSEQ_TABLE, Host=="AQ")
d=phyloseq::distance(PRUNED,'bray')
perm<-how(nperm=10000)
setBlocks(perm)<-with(dfC, Genotype)
ADONIS2<-adonis2(d~SamplingTimepoint, data=dfC, permutations = perm, method = "bray")
ADONIS2 #0.3705

Genotype<-get_variable(PRUNED, "Genotype")
HOSTANOSIM<-anosim(distance(PRUNED,"bray"), Genotype)
HOSTANOSIM$signif # p=0.001
HOSTANOSIM$statistic #R=0.6156379

#CO
dfC<-df%>%filter(Host=="CO")
PRUNED<-phyloseq::subset_samples(PHYLOSEQ_TABLE, Host=="CO")
d=phyloseq::distance(PRUNED,'bray')
perm<-how(nperm=10000)
setBlocks(perm)<-with(dfC, Genotype)
ADONIS2<-adonis2(d~SamplingTimepoint, data=dfC, permutations = perm, method = "bray")
ADONIS2 # 0.07699 .

Genotype<-get_variable(PRUNED, "Genotype")
HOSTANOSIM<-anosim(distance(PRUNED,"bray"), Genotype)
HOSTANOSIM$signif # p=0.001
HOSTANOSIM$statistic #R=0.8559671

#CY
dfC<-df%>%filter(Host=="CY")
PRUNED<-phyloseq::subset_samples(PHYLOSEQ_TABLE, Host=="CY")
d=phyloseq::distance(PRUNED,'bray')
perm<-how(nperm=10000)
setBlocks(perm)<-with(dfC, Genotype)
ADONIS2<-adonis2(d~SamplingTimepoint, data=dfC, permutations = perm, method = "bray")
ADONIS2 # 0.6689

Genotype<-get_variable(PRUNED, "Genotype")
HOSTANOSIM<-anosim(distance(PRUNED,"bray"), Genotype)
HOSTANOSIM$signif # p=0.001
HOSTANOSIM$statistic #R=0.8049383

#IB
dfC<-df%>%filter(Host=="IB")
PRUNED<-phyloseq::subset_samples(PHYLOSEQ_TABLE, Host=="IB")
d=phyloseq::distance(PRUNED,'bray')
perm<-how(nperm=10000)
setBlocks(perm)<-with(dfC, Genotype)
ADONIS2<-adonis2(d~SamplingTimepoint, data=dfC, permutations = perm, method = "bray")
ADONIS2 # 0.6266

Genotype<-get_variable(PRUNED, "Genotype")
HOSTANOSIM<-anosim(distance(PRUNED,"bray"), Genotype)
HOSTANOSIM$signif # p=0.001
HOSTANOSIM$statistic #R=0.5522634

#IR
dfC<-df%>%filter(Host=="IR")
PRUNED<-phyloseq::subset_samples(PHYLOSEQ_TABLE, Host=="IR")
d=phyloseq::distance(PRUNED,'bray')
perm<-how(nperm=10000)
setBlocks(perm)<-with(dfC, Genotype)
ADONIS2<-adonis2(d~SamplingTimepoint, data=dfC, permutations = perm, method = "bray")
ADONIS2 # 0.8046

Genotype<-get_variable(PRUNED, "Genotype")
HOSTANOSIM<-anosim(distance(PRUNED,"bray"), Genotype)
HOSTANOSIM$signif # p=0.001
HOSTANOSIM$statistic #R=0.9168724

#ST
dfC<-df%>%filter(Host=="ST")
PRUNED<-phyloseq::subset_samples(PHYLOSEQ_TABLE, Host=="ST")
d=phyloseq::distance(PRUNED,'bray')
perm<-how(nperm=10000)
setBlocks(perm)<-with(dfC, Genotype)
ADONIS2<-adonis2(d~SamplingTimepoint, data=dfC, permutations = perm, method = "bray")
ADONIS2 # 0.0006999

Genotype<-get_variable(PRUNED, "Genotype")
HOSTANOSIM<-anosim(distance(PRUNED,"bray"), Genotype)
HOSTANOSIM$signif # p=0.001
HOSTANOSIM$statistic #R=0.9744856

#DISPERSION -  multivariate homogeneity of group dispersion (variances).
Host_Treatment_group<-get_variable(PHYLOSEQ_TABLE, "Host_Treatment")
BETADISP<-betadisper(d, Host_Treatment_group,type=c("centroid"))
anova(BETADISP) #significant <2.2e-16
capture.output(anova(BETADISP),file="BETADISP_HostTreatment.doc")
plot(BETADISP)
boxplot(BETADISP)
permutest(BETADISP, pairwise = TRUE, permutations = 10000) # significant (however-pairwise comparison shows no sign. differnce between Control und Treatment group within a Host)
capture.output(permutest(BETADISP, pairwise = TRUE, permutations = 10000),file="BETADISP_HostTreatment_pairwise.doc")
BETADISP_HSD<-TukeyHSD(BETADISP)
BETADISP_HSD
capture.output(BETADISP_HSD,file="BETADISP_TukeyHSD.doc")
plot(BETADISP_HSD)

Host_group<-get_variable(PHYLOSEQ_TABLE, "Host")
BETADISP<-betadisper(d, Host_group,type=c("centroid"))
anova(BETADISP) #significant <2.2e-16
capture.output(anova(BETADISP),file="BETADISP_Host.doc")
plot(BETADISP)
boxplot(BETADISP)
permutest(BETADISP, pairwise = TRUE, permutations = 10000) # for all Hosts significant 
BETADISP_HSD<-TukeyHSD(BETADISP)
BETADISP_HSD
plot(BETADISP_HSD)

Treatment_group<-get_variable(PHYLOSEQ_TABLE, "Treatment")
BETADISP<-betadisper(d, Treatment_group,type=c("centroid"))
anova(BETADISP) #not sign
plot(BETADISP)
boxplot(BETADISP)
permutest(BETADISP, pairwise = TRUE, permutations = 10000) # not significant
BETADISP_HSD<-TukeyHSD(BETADISP)
BETADISP_HSD
plot(BETADISP_HSD)

#  multivariate homogeneity of group dispersion (variances).
TreatmentComb_group<-get_variable(PHYLOSEQ_TABLE, "TreatmentComb")
BETADISP<-betadisper(d, TreatmentComb_group,type=c("centroid"))
anova(BETADISP)
boxplot(BETADISP)
permutest(BETADISP, pairwise = TRUE, permutations = 10000)
BETADISP_HSD<-TukeyHSD(BETADISP)
BETADISP_HSD
plot(BETADISP_HSD)

#  multivariate homogeneity of group dispersion (variances)
d=phyloseq::distance(PHYLOSEQ_TABLE,'bray')
TreatmentComb_group<-get_variable(PHYLOSEQ_TABLE, "Genotype")
BETADISP<-betadisper(d, TreatmentComb_group,type=c("centroid"))
anova(BETADISP)
boxplot(BETADISP)
(BETADISP)
permutest(BETADISP, pairwise = TRUE, permutations = 10000)
BETADISP_HSD<-TukeyHSD(BETADISP)
BETADISP_HSD
plot(BETADISP_HSD)

### graph Betadisper
Host_Treatment_group<-get_variable(PHYLOSEQ_TABLE, "Host_Treatment")
BETADISP<-betadisper(d, Host_Treatment_group,type=c("centroid"))

Disp<-as.data.frame(BETADISP$distances) #from Code_short
Sample_ID<-rownames(Disp)
Disp<-cbind(Sample_ID, Disp)
Disp_new<-right_join(METADATA, Disp, by="Sample_ID")
Disp_new<-Disp_new%>% group_by(Host) %>% arrange(desc(BETADISP$distances), .by_group=TRUE)
Disp_new$TreatmentComb<-factor(Disp_new$TreatmentComb, levels=c("Control", "before Pulse", "after Pulse"))

plot<-ggplot(Disp_new, aes(y=(BETADISP$distances),x=Host_Treatment))+
  geom_boxplot(aes(y=BETADISP$distances, x=Host_Treatment))+
  geom_point(aes(color=Host_Treatment, shape=SamplingTimepoint), size=3)+
  scale_colour_manual(values=c(AQ_Control="#73d09f",AQ_Treatment="#50916f",
                               CO_Control="#aaaaaa",CO_Treatment="#37677E",
                               CY_Control="#a95382",CY_Treatment="#45416b",
                               IB_Control="#f9e07c",IB_Treatment="#ddba2d",
                               IR_Control="#ae1c1c",IR_Treatment="#6b0202",
                               ST_Control="#fab74f",ST_Treatment="#ff8533"),
                      breaks=c("AQ_Control","AQ_Treatment",
                               "CO_Control","CO_Treatment",
                               "CY_Control","CY_Treatment",
                               "IB_Control","IB_Treatment",
                               "IR_Control","IR_Treatment",
                               "ST_Control","ST_Treatment"),
                      labels=c("AQ Control","AQ Disturbance",
                               "CO Control","CO Disturbance",
                               "CY Control","CY Disturbance",
                               "IB Control","IB Disturbance",
                               "IR Control","IR Disturbance",
                               "ST Control","ST Disturbance"))+
  theme_classic()+
  coord_flip()+
  guides(color=guide_legend(title="host and treatments"), shape=guide_legend(title = "sampling day"))+
  labs(x= "host and treatments", y = "distance to group centroids")+
  scale_shape_discrete(breaks=c("before disturbance","24h after disturbance","168h after disturbance"),
                       labels=c("1", "11", "17"))

plot

pdf('Betadisper_HostTreatment.pdf', width=10, height=6)
print(plot)
graphics.off()

setEPS()
postscript('Betadisper_HostTreatment.eps', width=14, height=10)
plot
dev.off()

### Phyla of dataset ####
top20otus=names(sort(taxa_sums(PHYLOSEQ_TABLE), TRUE)[1:1000])
top20otus
taxtab20=cbind(tax_table(PHYLOSEQ_TABLE), OTU1000=NA)
taxtab20[top20otus, "OTU1000"] <- as(tax_table(PHYLOSEQ_TABLE)[top20otus, "Phylum"], "character")
tax_table(PHYLOSEQ_TABLE)<-tax_table(taxtab20)

merged<-merge_samples(PHYLOSEQ_TABLE, "Host_Treatment")
sample_data(merged)$Host_Treatment<-levels(sample_data(PHYLOSEQ_TABLE)$Host_Treatment)
merged=transform_sample_counts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

pal <- function(col, border = "light gray", ...)
{
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}
# create palett 
tol21rainbow1= c("#771122", "#AA4455", "#DD7788","#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA","#AAAA44", "#DDDD77")
tol21rainbow<- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598" ,"#ABDDA4" ,"#66C2A5", "#498699",
                 "#35a4dc", "#2083c5", "#4554a4","#65499d", "#8f3e97",
                 "#466CA6", "#3F477D","#7A6E99", "#865599","#BF5178", "#993A65","#991A51")

pal(tol21rainbow)

plot<-plot_bar(merged20, x="Host_Treatment", fill="Phylum")+
  geom_bar(aes(fill=Phylum, color=Phylum), stat="identity", position="stack")+
  theme_classic()+
  scale_fill_manual(values = tol21rainbow)+
  scale_color_manual(values = tol21rainbow)+
  theme(legend.position="right", legend.direction="vertical", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
plot


pdf('Phyla_Barchart.pdf', width=8, height=6)
print(plot)
graphics.off()

setEPS()
postscript('Phyla_Barchart.eps')
plot
dev.off()

### Class of dataset ####
top20otus=names(sort(taxa_sums(PHYLOSEQ_TABLE), TRUE)[1:100])
top20otus
taxtab20=cbind(tax_table(PHYLOSEQ_TABLE), OTU1000=NA)
taxtab20[top20otus, "OTU1000"] <- as(tax_table(PHYLOSEQ_TABLE)[top20otus, "Class"], "character")
tax_table(PHYLOSEQ_TABLE)<-tax_table(taxtab20)

merged<-merge_samples(PHYLOSEQ_TABLE, "Host_Treatment")
sample_data(merged)$Host_Treatment<-levels(sample_data(PHYLOSEQ_TABLE)$Host_Treatment)
merged=transform_sample_counts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

pal <- function(col, border = "light gray", ...)
{
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}
# create palett 
tol21rainbow1= c("#771122", "#AA4455", "#DD7788","#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA","#AAAA44", "#DDDD77")
tol21rainbow<- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598" ,"#ABDDA4" ,"#66C2A5", "#498699",
                 "#35a4dc", "#2083c5", "#4554a4","#65499d", "#8f3e97",
                 "#466CA6", "#3F477D","#7A6E99", "#865599","#BF5178", "#993A65","#991A51")

pal(tol21rainbow)

plot<-plot_bar(merged20, x="Host_Treatment", fill="Class")+
  geom_bar(aes(fill=Class, color=Class), stat="identity", position="stack")+
  theme_classic()+
  scale_fill_manual(values = tol21rainbow)+
  scale_color_manual(values = tol21rainbow)+
  theme(legend.position="right", legend.direction="vertical", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
plot


pdf('Class_Barchart.pdf', width=8, height=6)
print(plot)
graphics.off()

#### data for alluvial graph in RAWGraphs ####
# https://rawgraphs.io

#### AQ ####
AQ=subset_samples(PHYLOSEQ_TABLE, Host=="AQ")
AQ
AQ<-filter_taxa(AQ, function(x) sum(x) > 0.000000000, prune=TRUE)
AQ

top20otus=names(sort(taxa_sums(AQ), TRUE)[1:10])
top20otus
taxtab20=cbind(tax_table(AQ), OTU10=NA)
taxtab20[top20otus, "OTU10"] <- as(tax_table(AQ)[top20otus, "Phylum"], "character")
tax_table(AQ)<-tax_table(taxtab20)

TEST<-as.data.frame(tax_table(merged20))
TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="AQ_TEST.csv")

merged<-merge_samples(AQ, "Host")
sample_data(merged)$Host<-levels(sample_data(AQ)$Host)
merged=transform_sample_counts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

TEST<-as.data.frame(tax_table(merged20))
TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="AQ_TEST.csv")

# Nitrosococcus rel abundance
NITRO<-subset_taxa(AQ, Genus=="D_5__Nitrosococcus")
NITRO
NITROmelt<-psmelt(NITRO)

Nitrosococcus_com<-get_variable(NITRO, "Genotype")
HOSTANOSIM<-anosim(distance(NITRO,"bray"), Nitrosococcus_com)
HOSTANOSIM$signif # 0.001 sign
HOSTANOSIM$statistic # R=0.7127572



scaleFUN <- function(x) sprintf("%.0f", x)

relabundplot<-ggplot(NITROmelt, aes(SamplingTimepoint,Abundance))+
  geom_col()+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  coord_cartesian(ylim=c(0,100))+
  labs(y = "rel.abundance [%]")

relabundplot

### Oligo 100 percent ####
top20otus=names(sort(taxa_sums(NITRO), TRUE)[1:4])
top20otus
taxtab20=cbind(tax_table(NITRO), OTU15=NA)
taxtab20[top20otus, "OTU15"] <- as(tax_table(NITRO)[top20otus, "OTU15"], "character")
tax_table(NITRO)<-tax_table(taxtab20)

sample_data(NITRO)

merged20 = prune_taxa(top20otus, NITRO)
phy <- transform_sample_counts(merged20,function(x)100*x/sum(x))

tax_table(phy)

# create palett (20 colours) - "#DD7788" can be used as 21 if needed
tol21rainbow<- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598" ,"#ABDDA4" ,"#66C2A5","#498699", "#466CA6", "#3F477D","#7A6E99", "#865599","#BF5178", "#993A65","#991A51")
tol10<-c("#69bf69","#5e4fa2","#3288bd","#66c2a5","#abdda4","#e6f598","#fee08b","#fdae61","#f46d43","#d53e4f","#9e0142")


oligoplot<-plot_bar(phy, x="SamplingTimepoint", fill="OTUs")+
  geom_bar(aes(fill=OTUs, color=OTUs), stat="identity", position="stack")+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  theme(strip.text.x = element_blank(), text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0), labels=c("before disturbance"="before","24h after disturbance"="24h after","168h after disturbance"="168h after")) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  labs(y = "rel.abundance [%]")+
  theme(legend.position="none")

oligoplot

relabundplot

library("cowplot")
combplot<-ggdraw() +
  draw_plot(relabundplot, x = 0, y = 0.55, width = 1, height = .45)+ 
  draw_plot(oligoplot, x = 0, y = 0, width = 1, height = .55) 
draw_plot_label(label = c("B", "C"), size = 15,
                x = c(0, 0), y = c(1, 0.7))

combplot

pdf('AQ_Nitrosococcus_Oligotyping_Genotype.pdf', width=7, height=6)
print(combplot)
graphics.off()

tiff(filename='AQ_Nitrosococcus_Oligotyping_Genotype.tiff', width=7, height=6, units="in",
     pointsize=16, compression="lzw", bg="white", res=600)
print(combplot)
dev.off()

setEPS()
postscript('AQ_Nitrosococcus_Oligotyping_Genotype.eps', width=7, height=6)
combplot
dev.off()


#### CO ####
CO=subset_samples(PHYLOSEQ_TABLE, Host=="CO")
CO
CO<-filter_taxa(CO, function(x) sum(x) > 0.000000000, prune=TRUE)
CO

top20otus=names(sort(taxa_sums(CO), TRUE)[1:10])
top20otus
taxtab20=cbind(tax_table(CO), OTU10=NA)
taxtab20[top20otus, "OTU10"] <- as(tax_table(CO)[top20otus, "Phylum"], "character")
tax_table(CO)<-tax_table(taxtab20)

merged<-merge_samples(CO, "Host")
sample_data(merged)$Host<-levels(sample_data(CO)$Host)
merged=transform_sample_counts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="CO_alluvial.csv")

# PAUC34f rel abundance
PAUC<-subset_taxa(CO, Phylum=="D_1__PAUC34f")
PAUC
PAUCmelt<-psmelt(PAUC)

scaleFUN <- function(x) sprintf("%.0f", x)

relabundplot<-ggplot(PAUCmelt, aes(SamplingTimepoint,Abundance))+
  geom_col()+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  coord_cartesian(ylim=c(0,100))+
  labs(y = "rel.abundance [%]")

relabundplot

### Oligo 100 percent ####
top20otus=names(sort(taxa_sums(PAUC), TRUE)[1:3])
top20otus
taxtab20=cbind(tax_table(PAUC), OTU15=NA)
taxtab20[top20otus, "OTU15"] <- as(tax_table(PAUC)[top20otus, "OTU15"], "character")
tax_table(PAUC)<-tax_table(taxtab20)

sample_data(PAUC)

merged20 = prune_taxa(top20otus, PAUC)
phy <- transform_sample_counts(merged20,function(x)100*x/sum(x))

tax_table(phy)

oligoplot<-plot_bar(phy, x="SamplingTimepoint", fill="OTUs")+
  geom_bar(aes(fill=OTUs, color=OTUs), stat="identity", position="stack")+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  theme(strip.text.x = element_blank(), text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0), labels=c("before disturbance"="before","24h after disturbance"="24h after","168h after disturbance"="168h after")) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  labs(y = "rel.abundance [%]")
theme(legend.position="none")

oligoplot

relabundplot

library("cowplot")
combplot<-ggdraw() +
  draw_plot(relabundplot, x = 0, y = 0.55, width = 1, height = .45)+ 
  draw_plot(oligoplot, x = 0, y = 0, width = 1, height = .55) 
draw_plot_label(label = c("B", "C"), size = 15,
                x = c(0, 0), y = c(1, 0.7))

combplot

pdf('CO_PAUC34f_Oligotyping_Genotype.pdf', width=7, height=6)
print(combplot)
graphics.off()

tiff(filename='CO_PAUC34f_Oligotyping_Genotype.tiff', width=7, height=6, units="in",
     pointsize=16, compression="lzw", bg="white", res=600)
print(combplot)
dev.off()

setEPS()
postscript('CO_PAUC34f_Oligotyping_Genotype.eps', width=7, height=6)
combplot
dev.off()


#### CY ####
CY=subset_samples(PHYLOSEQ_TABLE, Host=="CY")
CY
CY<-filter_taxa(CY, function(x) sum(x) > 0.000000000, prune=TRUE)
CY

top20otus=names(sort(taxa_sums(CY), TRUE)[1:10])
top20otus
taxtab20=cbind(tax_table(CY), OTU10=NA)
taxtab20[top20otus, "OTU10"] <- as(tax_table(CY)[top20otus, "Phylum"], "character")
tax_table(CY)<-tax_table(taxtab20)

merged<-merge_samples(CY, "Host")
sample_data(merged)$Host<-levels(sample_data(CY)$Host)
merged=transform_sample_CYunts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="CY_alluvial.csv")

# Cyano FamilyI rel abundance
Cyano<-subset_taxa(CY, Family=="D_4__FamilyI")
Cyano
Cyanomelt<-psmelt(Cyano)

scaleFUN <- function(x) sprintf("%.0f", x)

relabundplot<-ggplot(Cyanomelt, aes(SamplingTimepoint,Abundance))+
  geom_col()+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  coord_cartesian(ylim=c(0,100))+
  labs(y = "rel.abundance [%]")

relabundplot

### Oligo 100 percent ####
top20otus=names(sort(taxa_sums(Cyano), TRUE)[1:3])
top20otus
taxtab20=cbind(tax_table(Cyano), OTU15=NA)
taxtab20[top20otus, "OTU15"] <- as(tax_table(Cyano)[top20otus, "OTU15"], "character")
tax_table(Cyano)<-tax_table(taxtab20)

sample_data(Cyano)

merged20 = prune_taxa(top20otus, Cyano)
phy <- transform_sample_counts(merged20,function(x)100*x/sum(x))

tax_table(phy)

oligoplot<-plot_bar(phy, x="SamplingTimepoint", fill="OTUs")+
  geom_bar(aes(fill=OTUs, color=OTUs), stat="identity", position="stack")+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  theme(strip.text.x = element_blank(), text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0), labels=c("before disturbance"="before","24h after disturbance"="24h after","168h after disturbance"="168h after")) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  labs(y = "rel.abundance [%]")+
  theme(legend.position="none")

oligoplot

relabundplot

library("cowplot")
combplot<-ggdraw() +
  draw_plot(relabundplot, x = 0, y = 0.55, width = 1, height = .45)+ 
  draw_plot(oligoplot, x = 0, y = 0, width = 1, height = .55) 
draw_plot_label(label = c("B", "C"), size = 15,
                x = c(0, 0), y = c(1, 0.7))

combplot

pdf('CY_Cyano_Oligotyping_Genotype.pdf', width=7, height=6)
print(combplot)
graphics.off()

tiff(filename='CY_Cyano_Oligotyping_Genotype.tiff', width=7, height=6, units="in",
     pointsize=16, compression="lzw", bg="white", res=600)
print(combplot)
dev.off()

setEPS()
postscript('CY_Cyano_Oligotyping_Genotype.eps', width=7, height=6)
combplot
dev.off()

#### IB ####
IB=subset_samples(PHYLOSEQ_TABLE, Host=="IB")
IB
IB<-filter_taxa(IB, function(x) sum(x) > 0.000000000, prune=TRUE)
IB

top20otus=names(sort(taxa_sums(IB), TRUE)[1:10])
top20otus
taxtab20=cbind(tax_table(IB), OTU10=NA)
taxtab20[top20otus, "OTU10"] <- as(tax_table(IB)[top20otus, "Phylum"], "character")
tax_table(IB)<-tax_table(taxtab20)

merged<-merge_samples(IB, "Host")
sample_data(merged)$Host<-levels(sample_data(IB)$Host)
merged=transform_sample_IBunts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="IB_alluvial.csv")

# Gamma rel abundance
Gamma<-subset_taxa(IB, Class=="D_2__Gammaproteobacteria")
Gamma
Gammamelt<-psmelt(Gamma)

scaleFUN <- function(x) sprintf("%.0f", x)

relabundplot<-ggplot(Gammamelt, aes(SamplingTimepoint,Abundance))+
  geom_col()+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  coord_cartesian(ylim=c(0,100))+
  labs(y = "rel.abundance [%]")

relabundplot

### Oligo 100 percent ####
top20otus=names(sort(taxa_sums(Gamma), TRUE)[1:3])
top20otus
taxtab20=cbind(tax_table(Gamma), OTU15=NA)
taxtab20[top20otus, "OTU15"] <- as(tax_table(Gamma)[top20otus, "OTU15"], "character")
tax_table(Gamma)<-tax_table(taxtab20)

sample_data(Gamma)

merged20 = prune_taxa(top20otus, Gamma)
phy <- transform_sample_counts(merged20,function(x)100*x/sum(x))

tax_table(phy)

oligoplot<-plot_bar(phy, x="SamplingTimepoint", fill="OTUs")+
  geom_bar(aes(fill=OTUs, color=OTUs), stat="identity", position="stack")+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  theme(strip.text.x = element_blank(), text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0), labels=c("before disturbance"="before","24h after disturbance"="24h after","168h after disturbance"="168h after")) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  labs(y = "rel.abundance [%]")+
  theme(legend.position="none")

oligoplot

relabundplot

library("cowplot")
combplot<-ggdraw() +
  draw_plot(relabundplot, x = 0, y = 0.55, width = 1, height = .45)+ 
  draw_plot(oligoplot, x = 0, y = 0, width = 1, height = .55) 
draw_plot_label(label = c("B", "C"), size = 15,
                x = c(0, 0), y = c(1, 0.7))

combplot

pdf('IB_Gamma_Oligotyping_Genotype.pdf', width=7, height=6)
print(combplot)
graphics.off()

tiff(filename='IB_Gamma_Oligotyping_Genotype.tiff', width=7, height=6, units="in",
     pointsize=16, compression="lzw", bg="white", res=600)
print(combplot)
dev.off()

setEPS()
postscript('IB_Gamma_Oligotyping_Genotype.eps', width=7, height=6)
combplot
dev.off()

#### IR ####
IR=subset_samples(PHYLOSEQ_TABLE, Host=="IR")
IR
IR<-filter_taxa(IR, function(x) sum(x) > 0.000000000, prune=TRUE)
IR

top20otus=names(sort(taxa_sums(IR), TRUE)[1:10])
top20otus
taxtab20=cbind(tax_table(IR), OTU10=NA)
taxtab20[top20otus, "OTU10"] <- as(tax_table(IR)[top20otus, "Phylum"], "character")
tax_table(IR)<-tax_table(taxtab20)

merged<-merge_samples(IR, "Host")
sample_data(merged)$Host<-levels(sample_data(IR)$Host)
merged=transform_sample_IRunts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="IR_alluvial.csv")

# Rhodothermaceae rel abundance
Rhodo<-subset_taxa(IR, Family=="D_4__Rhodothermaceae")
Rhodo
Rhodomelt<-psmelt(Rhodo)

scaleFUN <- function(x) sprintf("%.0f", x)

relabundplot<-ggplot(Rhodomelt, aes(SamplingTimepoint,Abundance))+
  geom_col()+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  coord_cartesian(ylim=c(0,100))+
  labs(y = "rel.abundance [%]")

relabundplot

### Oligo 100 percent ####
top20otus=names(sort(taxa_sums(Rhodo), TRUE)[1:3])
top20otus
taxtab20=cbind(tax_table(Rhodo), OTU15=NA)
taxtab20[top20otus, "OTU15"] <- as(tax_table(Rhodo)[top20otus, "OTU15"], "character")
tax_table(Rhodo)<-tax_table(taxtab20)

sample_data(Rhodo)

merged20 = prune_taxa(top20otus, Rhodo)
phy <- transform_sample_counts(merged20,function(x)100*x/sum(x))

tax_table(phy)

oligoplot<-plot_bar(phy, x="SamplingTimepoint", fill="OTUs")+
  geom_bar(aes(fill=OTUs, color=OTUs), stat="identity", position="stack")+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  theme(strip.text.x = element_blank(), text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0), labels=c("before disturbance"="before","24h after disturbance"="24h after","168h after disturbance"="168h after")) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  labs(y = "rel.abundance [%]")+
  theme(legend.position="none")

oligoplot

relabundplot

library("cowplot")
combplot<-ggdraw() +
  draw_plot(relabundplot, x = 0, y = 0.55, width = 1, height = .45)+ 
  draw_plot(oligoplot, x = 0, y = 0, width = 1, height = .55) 
draw_plot_label(label = c("B", "C"), size = 15,
                x = c(0, 0), y = c(1, 0.7))

combplot

pdf('IR_Rhodo_Oligotyping_Genotype.pdf', width=7, height=6)
print(combplot)
graphics.off()

tiff(filename='IR_Rhodo_Oligotyping_Genotype.tiff', width=7, height=6, units="in",
     pointsize=16, compression="lzw", bg="white", res=600)
print(combplot)
dev.off()

setEPS()
postscript('IR_Rhodo_Oligotyping_Genotype.eps', width=7, height=6)
combplot
dev.off()





#### ST ####
ST=subset_samples(PHYLOSEQ_TABLE, Host=="ST")
ST
ST<-filter_taxa(ST, function(x) sum(x) > 0.000000000, prune=TRUE)
ST

top20otus=names(sort(taxa_sums(ST), TRUE)[1:10])
top20otus
taxtab20=cbind(tax_table(ST), OTU10=NA)
taxtab20[top20otus, "OTU10"] <- as(tax_table(ST)[top20otus, "Phylum"], "character")
tax_table(ST)<-tax_table(taxtab20)

merged<-merge_samples(ST, "Host")
sample_data(merged)$Host<-levels(sample_data(ST)$Host)
merged=transform_sample_STunts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

TESTmelt<-psmelt(merged20)
write.csv(TESTmelt, file="ST_alluvial.csv")

# Nitrospira rel abundance
SPIRA<-subset_taxa(ST, Genus=="D_5__Nitrospira")
SPIRA
Spiramelt<-psmelt(SPIRA)

scaleFUN <- function(x) sprintf("%.0f", x)

relabundplot<-ggplot(Spiramelt, aes(SamplingTimepoint,Abundance))+
  geom_col()+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  coord_cartesian(ylim=c(0,100))+
  labs(y = "rel.abundance [%]")

relabundplot

### Oligo 100 percent ####
top20otus=names(sort(taxa_sums(SPIRA), TRUE)[1:2])
top20otus
taxtab20=cbind(tax_table(SPIRA), OTU15=NA)
taxtab20[top20otus, "OTU15"] <- as(tax_table(SPIRA)[top20otus, "OTU15"], "character")
tax_table(SPIRA)<-tax_table(taxtab20)

sample_data(SPIRA)

merged20 = prune_taxa(top20otus, SPIRA)
phy <- transform_sample_counts(merged20,function(x)100*x/sum(x))

tax_table(phy)

oligoplot<-plot_bar(phy, x="SamplingTimepoint", fill="OTUs")+
  geom_bar(aes(fill=OTUs, color=OTUs), stat="identity", position="stack")+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")+
  theme_classic()+
  facet_wrap(Treatment~Genotype, ncol=6)+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  theme(strip.text.x = element_blank(), text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0), labels=c("before disturbance"="before","24h after disturbance"="24h after","168h after disturbance"="168h after")) + scale_y_continuous(labels=scaleFUN, expand = c(0, 0))+
  labs(y = "rel.abundance [%]")+
  theme(legend.position="none")

oligoplot

relabundplot

library("cowplot")
combplot<-ggdraw() +
  draw_plot(relabundplot, x = 0, y = 0.55, width = 1, height = .45)+ 
  draw_plot(oligoplot, x = 0, y = 0, width = 1, height = .55) 
draw_plot_label(label = c("B", "C"), size = 15,
                x = c(0, 0), y = c(1, 0.7))

combplot

pdf('ST_Nitrospira_Oligotyping_Genotype.pdf', width=7, height=6)
print(combplot)
graphics.off()

tiff(filename='ST_Nitrospira_Oligotyping_Genotype.tiff', width=7, height=6, units="in",
     pointsize=16, compression="lzw", bg="white", res=600)
print(combplot)
dev.off()

setEPS()
postscript('ST_Nitrospira_Oligotyping_Genotype.eps', width=7, height=6)
combplot
dev.off()

#### Indicator Value Analysis ####
#not included in paper
library(labdsv)
library(MASS)
library(vegan)
library(cluster)
library(indicspecies)
library(permute)
library(phyloseq)
library(data.table)

OTU_table_rel<-as.data.frame(otu_table(PHYLOSEQ_TABLE))
OTU_table_rel<-as.data.frame(t(OTU_table_rel))
OTU_table_rel<-setDT(OTU_table_rel, keep.rownames = TRUE)
OTU_table_rel<-dplyr::rename(OTU_table_rel, Sample_ID=rn)

OTU_table_rel<-right_join(METADATA, OTU_table_rel, by="Sample_ID")
head(OTU_table_rel)

#indval
(INDVAL_OTUs_species_GC=(as.data.frame(OTU_table_rel[,15:6054])))
(INDVAL_Groups_Origin_GC=(as.character(OTU_table_rel$Host)))
INDVAL_Origin_GC=multipatt(INDVAL_OTUs_species_GC, INDVAL_Groups_Origin_GC, func="IndVal.g", duleg=FALSE, control=how(nperm=100))
summary(INDVAL_Origin_GC, indvalcomp=TRUE)
options(max.print=1000000000)
(summary.multipatt(INDVAL_Origin_GC, alpha = 0.05, indvalcomp=TRUE, At=.85, Bt=.85)) 
# Component ‘A’ is the probability that the surveyed
# site belongs to the target site group given the fact that the species has been found. 
# This conditional probability is called the specificity or positive predictive value of the species as indicator of the site group. 
# Component ‘B’ is the probability of finding the species in sites belonging to the site group. 
# This second conditional probability is called the fidelity or sensitivity of the species as indicator of the target site group
options(max.print=1000000000)
capture.output(summary.multipatt(INDVAL_Origin_GC, alpha = 0.05,At=0.85, Bt=0.85, indvalcomp=TRUE),file="IndVal_Hostspecific085.csv")

## graphs IndVal ####
#load datasets & merge them so that only indicator otus are left
IndValHOST<-read.csv('IndVallist.csv', header = TRUE, sep = ",", strip.white = TRUE)
head(IndValHOST)
IndValHOST$OTUs<-as.character(IndValHOST$OTUs)
INDVAL_PHYLO<-phyloseq::prune_taxa(IndValHOST$OTUs, PHYLOSEQ_TABLE)
INDVAL_PHYLO
tax_table(INDVAL_PHYLO)

#### Genus of the 15 most abundant indicators in the dataset ####
top20otus=names(sort(taxa_sums(IndVal_PHYLO), TRUE)[1:15])
top20otus
taxtab20=cbind(tax_table(IndVal_PHYLO), OTU15=NA)
taxtab20[top20otus, "OTU15"] <- as(tax_table(IndVal_PHYLO)[top20otus, "Order"], "character")
tax_table(IndVal_PHYLO)<-tax_table(taxtab20)

merged<-merge_samples(IndVal_PHYLO, "Host_Treatment")
sample_data(merged)$Host_Treatment<-levels(sample_data(IndVal_PHYLO)$Host_Treatment)
merged=transform_sample_counts(merged,function(x)100*x/sum(x))

merged20 = prune_taxa(top20otus, merged)

pal <- function(col, border = "light gray", ...)
{
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}
# create palett 
tol21rainbow1= c("#771122", "#AA4455", "#DD7788","#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA","#AAAA44", "#DDDD77")
tol21rainbow<- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598" ,"#ABDDA4" ,"#66C2A5","#498699", "#466CA6", "#3F477D","#7A6E99", "#865599","#BF5178", "#993A65","#991A51")

pal(tol21rainbow)
library(ggplot2)
plot<-plot_bar(merged, x="Host_Treatment", fill="Phylum")+
  geom_bar(aes(fill=Phylum, color=Phylum), stat="identity", position="stack")
scale_fill_manual(values = tol21rainbow)+
  scale_color_manual(values = tol21rainbow)+
  theme_classic()+
  theme(legend.position="right", legend.direction="vertical", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
plot


merged<-merge_samples(IndVal_PHYLO, "Host_Treatment")
sample_data(merged)$Host_Treatment<-levels(sample_data(IndVal_PHYLO)$Host_Treatment)
merged=transform_sample_counts(merged,function(x)100*x/sum(x))
plot<-plot_bar(merged20, x="Host_Treatment", fill="Order")+
  geom_bar(aes(fill=Order, color=Order), stat="identity", position="stack")
scale_fill_manual(values = tol21rainbow)+
  scale_color_manual(values = tol21rainbow)+
  theme_classic()+
  theme(legend.position="right", legend.direction="vertical", legend.title = element_blank(), axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1))+
  scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
plot

plot_heatmap(INDVAL_PHYLO, "NMDS", "bray", "Code_short", "Genus")
#####Full analysis for: Wright et al. Adaptive divergence in Heliconius eyes#####

#load packages
library(MASS)
library(lme4)
library(ggplot2)
library(stats)
library(effects)
library(splines)
library(xlsx)
library(car)
library(lattice)
library(grid)
library(gridExtra)
library(pbkrtest)
library(readxl)
library(dplyr)

#import data
eye <- read.csv("~/Documents/LMU/Eye morphology/selection manuscript/eye_morphology_cp_mp_f1.csv")


#is the data properly structured?
str(eye)
eye$type <- as.factor(eye$type)
eye$sex <- as.factor(eye$sex)
eye$brood <- as.factor(eye$brood)
eye$wild_insectary <- as.factor(eye$wild_insectary)
eye$raised <- as.factor(eye$raised)
eye$observer <- as.factor(eye$observer)
str(eye) #now correct

######new composite variables of facet count & corneal area####
#when left eye was missing or damaged, right is used

eye$whole_count <- ifelse(is.na(eye$L_whole_count), eye$R_whole_count, eye$L_whole_count)
eye$area <- ifelse(is.na(eye$L_area), eye$R_area, eye$L_area)

#log-transformed variables
eye$area_log <- (log10(eye$area))
hist(eye$area_log)

eye$tibia_length_log <- (log10(eye$tibia_length))
hist((eye$tibia_length_log))

eye$whole_count_log <- (log10(eye$whole_count))
hist((eye$whole_count_log))

#some individuals are missing body size estimate (tibia length)
#remove from the dataset
eye <- eye[!(is.na(eye$tibia_length) | eye$tibia_length==""), ]

#how many per species & sex
t.first <- eye[match(unique(eye$ID), eye$ID),]
t.first %>%
  group_by(collection) %>%
  summarize(count = n())

#CP    female    19
#CP    male      15
#CPxMP female    15
#CPxMP male      15
#MP    female    15
#MP    male      18
#MPxCP female    15
#MPxCP male      15


#dataset without F1 hybrids
eye2 <- eye[-which(eye$type=="CPxMP"|
                  eye$type=="MPxCP"),]


#how many per species & sex
t.second <- eye2[match(unique(eye2$ID), eye2$ID),]
t.second %>%
  group_by(type,sex) %>%
  summarize(count = n())

#CP    female    19
#CP    male      15
#MP    female    15
#MP    male      18


####Correlations between left and right sides####
#using all processed samples (pure species + F1 hybrids)

cor.test(eye$L_whole_count,eye$R_whole_count, method="pearson")
cor.test(eye$L_area,eye$R_area, method="pearson")

#import dataset with L & R hind tibia lengths
legs <- read.csv("~/Documents/LMU/Eye morphology/selection manuscript/legs.csv")
cor.test(legs$L_length,legs$R_length, method="pearson")


####Correlations between body size, facet count, corneal area####
#without hybrids, onle pure species
cor.test(eye2$tibia_length,eye2$whole_count, method="pearson")
cor.test(eye2$area,eye2$whole_count, method="pearson")
cor.test(eye2$tibia_length,eye2$area, method="pearson")


####WHAT IS THE RELATIONSHIP BETWEEN FACET COUNT AND CORNEAL AREA?####
#pure species only

hist(eye2$area_log)

d1 <- lm(area_log~whole_count_log*type*sex, data=eye2)
plot(d1, which=1)
plot(d1, which=2) 
hist(resid(d1)) 

drop1(d1, test="Chisq") #3-way interaction 0.098, remove
d2 <- update(d1,.~.-whole_count_log:type:sex)
drop1(d2, test="Chisq") #all 2-way interactions ns
d3 <- update(d2,.~.-type:sex)
drop1(d3, test="Chisq")
d4 <- update(d3,.~.-whole_count_log:sex)
drop1(d4, test="Chisq")
d5 <- update(d4,.~.-whole_count_log:type)
drop1(d5, test="Chisq") #sex ns
d6 <- update(d5,.~.-sex)
drop1(d6, test="Chisq")

plot(allEffects(d6))

Anova(d6, test="Chisq")
#facet count p<0.001
#type p<0.001

Anova(d4, test="Chisq")
#facet count:species n.s. p=0.12


####ANALYSES OF FACET COUNT - C & M only####

#full statistical model with only C & M
hist(eye2$whole_count_log) 

m1 <- lm(whole_count_log~type*sex+tibia_length_log, data=eye2)
plot(m1, which=1)
plot(m1, which=2)
hist(resid(m1))

drop1(m1, test="Chisq") #type:sex interaction ns
m2 <- update(m1,.~.-type:sex) 
drop1(m2, test="Chisq") #all three significant

plot(allEffects(m2))
plot(m2, which = 2) 
hist(resid(m2))

Anova(m2, test.statistic="Chisq") 
#type p=0.0087
#sex p>0.0001
#tibia p>0.0001

Anova(m1, test.statistic="Chisq") 
#type:sex n.s. p=0.78

m2.1 <- update(m2,.~.-type)
PBmodcomp(m2,m2.1,nsim = 1000) #same result with parametric bootstrapping as with Anova


#####plot of facet count based on MAM####
#using tibia length

#theme settings for plots
require(grid)
pub<- theme_update(
  panel.grid.major=element_line(colour=NA),
  panel.grid.minor=element_line(colour=NA),
  panel.background = element_rect(colour = NA,fill=NA,linewidth = 0.5),
  panel.border = element_rect(linetype = "solid", colour = "black",fill=NA),
  axis.line.x = element_line(color="black"),
  axis.line.y = element_line(color="black"),
  axis.title.x=element_text(size=15,face="bold",hjust=0.5,vjust=0.5,angle=0),
  axis.title.y=element_text(size=15,face="bold",hjust=0.5,vjust=1,angle=90),
  axis.text.x=element_text(colour="black",angle=0,size=15),
  axis.text.y=element_text(colour="black",angle=0,size=15),
  axis.ticks=element_line(colour="black",linewidth=0.5))


#get estimated marginal means of the final model (m2)
######**species**####

library(tidyverse)
library(emmeans)
library(multcomp)
library(multcompView)

model_means <- emmeans(object = m2, specs = ~ type, type = "response") 

#p-value comparison table
pwpm(model_means, adjust="bonferroni",  diffs = F)

#effect size
eff_size(model_means, sigma = sigma(m2), edf = df.residual(m2))

#add letters to each mean to indicate difference
model_means_cld <- cld(object = model_means,
                       reversed=T,
                       adjust="bonferroni",
                       Letters = letters,
                       alpha = 0.05)

#set marginal means at a dataframe
facet <- as.data.frame(summary(model_means_cld))

#species

#color selection for consistenty with other papers
library(RColorBrewer)
display.brewer.pal(12, 'Paired')
brewer.pal(12, 'Paired')

facet$type <- factor(facet$type, c("CP", "MP"))

species <- ggplot(facet, aes(type,emmean))+
  geom_errorbar(data=facet,mapping=aes(x=type, ymin=upper.CL, ymax=lower.CL), color="grey45",
                width=0, linewidth=1, position = position_dodge(0.4))+
  geom_point(data=facet, aes(type,emmean, color=type), size=4,
             position = position_dodge(0.4))+
  geom_text(aes(label = .group, y = emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=4.5) +
  labs(x="",y="e.m.m. log10 (facet count)")+
  scale_y_continuous(breaks=seq(4.11, 4.17, 0.02), limits=c(4.10, 4.174))+
  scale_x_discrete(labels=c("", "", ""))+
  scale_color_manual(values = c("#1F78B4","#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(plot.margin = margin(0.4,0.5,0.3,0.75, "cm"))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(axis.text.x=element_blank())+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.title.y = element_text(vjust = 5))+
  ggtitle("species effect") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


######**sex**####
model_means2 <- emmeans(object = m2, specs = ~ sex, type = "response") 

#p-value comparison table
pwpm(model_means2, adjust="bonferroni",  diffs = F)

#effect size
eff_size(model_means2, sigma = sigma(m2), edf = df.residual(m2))

#add letters to each mean to indicate difference
model_means_cld2 <- cld(object = model_means2,
                        reversed = T,
                        adjust="bonferroni",
                        Letters = letters,
                        alpha = 0.05)

#set marginal means at a dataframe
facet2 <- as.data.frame(summary(model_means_cld2))

#plot marginal means
facet2$sex <- factor(facet2$sex, c("male", "female"))

#sex

sex <- ggplot(facet2, aes(sex,emmean))+
  geom_errorbar(data=facet2,mapping=aes(x=sex, ymin=upper.CL, ymax=lower.CL), color="grey45",
                width=0, linewidth=1)+
  geom_point(data=facet2, aes(sex,emmean), size=4)+
  geom_text(aes(label = .group, y = emmean),
            vjust=0.3, hjust=-0.3,size=4.5) +
  labs(x="",y="")+
  scale_y_continuous(breaks=seq(4.11, 4.17, 0.02), limits=c(4.10, 4.174),
                     position="right")+
  scale_x_discrete(labels = c("", ""))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(legend.position=c(0.88,0.89), legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(0.4,1.05,0.3,0.85, "cm"))+
  theme(axis.text.x=element_blank())+
  theme(text = element_text(family = "Times New Roman"))+
  ggtitle("sex effect") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


######**species:sex**####
model_means3 <- emmeans(object = m2, specs = ~ type:sex, type = "response") 

#p-value comparison table
pwpm(model_means3, adjust="bonferroni",  diffs = F)

#add letters to each mean to indicate difference
model_means_cld3 <- cld(object = model_means3,
                        adjust = "bonferroni",
                        reversed=T,
                        Letters = letters,
                        alpha = 0.05)

#set marginal means at a dataframe
facet3 <- as.data.frame(summary(model_means_cld3))

#plot marginal means
facet3$sex <- factor(facet3$sex, c("male", "female"))
facet3$type <- factor(facet3$type, c("CP", "MP"))

#species:sex plot
combo1<-ggplot(facet3, aes(type,emmean))+
  geom_errorbar(data=facet3,mapping=aes(x=type, ymin=upper.CL, ymax=lower.CL, color=sex),
                width=0, linewidth=1, position = position_dodge(width=0.5))+
  geom_point(data=facet3, aes(type,emmean, color=sex), size=3,
             position = position_dodge2(width=0.5))+
  geom_text(aes(label = .group, y = emmean),
            vjust=0.3, hjust=-0.3,size=4.5, position = position_dodge2(width=0.5)) +
  labs(x="",y="e.m.m. log10 (facet count)", color="")+
  scale_y_continuous(breaks=seq(4.09, 4.21, 0.03), limits=c(4.08, 4.21))+
  scale_x_discrete(labels=c("H. cydno", "H. melp"))+
  scale_color_manual(values=c("black","grey65"))+
  theme(legend.position="", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(1,0.5,0.5,0.5, "cm"))+
  theme(legend.text=element_text(size=11))+
  theme(legend.key=element_blank())+
  theme(axis.text.x=element_text(face = "bold"))+
  theme(text = element_text(family = "Times New Roman"))


####ANALYSES OF CORNEAL AREA - C & M ONLY#### 
hist(eye2$area_log)

x1 <- lm(area_log~type*sex+tibia_length_log, data=eye2)
plot(x1,which=1)
plot(x1,which=2)
hist(resid(x1))

drop1(x1, test="Chisq") #interaction weakly trending (0.08)

Anova(x1, test="Chisq") #interaction n.s. p=0.102, remove
plot(allEffects(x1))

x2 <- update(x1,.~.-type:sex)
drop1(x2, test="Chisq") 

plot(x2,which=2) 
hist(resid(x2))

Anova(x2, test="Chisq")
#type p=0.004
#sex p>0.001
#tibia p>0.001

plot(allEffects(x2))


#####plot of corneal area based on MAM####

#get estimated marginal means of the final model (x2)

######**species**#####

model_means_x <- emmeans(object = x2, specs = ~ type, type = "response") 

#p-value comparison table
pwpm(model_means_x, adjust="bonferroni",  diffs = F)

#effect size
eff_size(model_means_x, sigma = sigma(x2), edf = df.residual(x2))

#add letters to each mean to indicate difference
model_means_cld_x <- cld(object = model_means_x,
                         reversed=T,
                         adjust="bonferroni",
                         Letters = letters,
                         alpha = 0.05)

#set marginal means at a dataframe
area <- as.data.frame(summary(model_means_cld_x))

#species

area$type <- factor(area$type, c("CP", "MP"))

species2 <- ggplot(area, aes(type,emmean))+
  geom_errorbar(data=area,mapping=aes(x=type, ymin=upper.CL, ymax=lower.CL), color="grey45",
                width=0, linewidth=1, position = position_dodge(0.4))+
  geom_point(data=area, aes(type,emmean,color=type), size=4,
             position = position_dodge(0.4))+
  geom_text(aes(label = .group, y = emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=4.5) +
  labs(x="",y="e.m.m. log10 (corneal area)")+
  scale_y_continuous(breaks=seq(0.82, 0.92, 0.03), limits=c(0.81, 0.92))+
  scale_x_discrete(labels=c("H. cydno", "H. melp"))+
  scale_color_manual(values = c("#1F78B4","#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(0.3,0.5,0,0.75, "cm"))+
  theme(axis.text.x=element_text(face = "bold"))+
  theme(text = element_text(family = "Times New Roman"))+
  theme(axis.title.y = element_text(vjust = 5))+
  ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


######**sex**####
model_means2_x <- emmeans(object = x2, specs = ~ sex, type = "response") 

#p-value comparison table
pwpm(model_means2_x, adjust="bonferroni",  diffs = F)

#effect size
eff_size(model_means2_x, sigma = sigma(x2), edf = df.residual(x2))

#add letters to each mean to indicate difference
model_means_cld2_x <- cld(object = model_means2_x,
                          reversed = T,
                          adjust="bonferroni",
                          Letters = letters,
                          alpha = 0.05)

#set marginal means at a dataframe
area2 <- as.data.frame(summary(model_means_cld2_x))

#plot marginal means
area2$sex <- factor(area2$sex, c("male", "female"))

#sex
sex2<- ggplot(area2, aes(sex,emmean))+
  geom_errorbar(data=area2,mapping=aes(x=sex, ymin=upper.CL, ymax=lower.CL), color="grey45",
                width=0, linewidth=1)+
  geom_point(data=area2, aes(sex,emmean), size=4)+
  geom_text(aes(label = .group, y = emmean),
            vjust=0.3, hjust=-0.3,size=4.5) +
  labs(x="",y="")+
  scale_y_continuous(breaks=seq(0.82, 0.92, 0.03), limits=c(0.81, 0.92),
                     position="right")+
  scale_x_discrete(labels = c("male", "female"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(legend.position=c(0.88,0.89), legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(0.3,1.05,0,0.85, "cm"))+
  theme(axis.text.x=element_text(face = "bold"))+
  theme(text = element_text(family = "Times New Roman"))+
  ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


######**species:sex**####
model_means3_x <- emmeans(object = x2, specs = ~ type:sex, type = "response") 

#p-value comparison table
pwpm(model_means3_x, adjust="bonferroni",  diffs = F)


#add letters to each mean to indicate difference
model_means_cld3_x <- cld(object = model_means3_x,
                          adjust = "bonferroni",
                          reversed = T,
                          Letters = letters,
                          alpha = 0.05)

#set marginal means at a dataframe
area3 <- as.data.frame(summary(model_means_cld3_x))


#plot marginal means
area3$sex <- factor(area3$sex, c("male", "female"))

#species:sex plot
combo2<-ggplot(area3, aes(type,emmean))+
  geom_errorbar(data=area3,mapping=aes(x=type, ymin=upper.CL, ymax=lower.CL, color=sex),
                width=0, linewidth=1, position = position_dodge(width=0.5))+
  geom_point(data=area3, aes(type,emmean, color=sex), size=3,
             position = position_dodge2(width=0.5))+
  geom_text(aes(label = .group, y = emmean),
            vjust=0.3, hjust=-0.3,size=4.5, position = position_dodge2(width=0.5)) +
  labs(x="",y="e.m.m. log10 (corneal area)", color="")+
  scale_y_continuous(breaks=seq(0.83, 0.91, 0.02), limits=c(0.825, 0.915),
                     position="left")+
  scale_x_discrete(labels=c("H. cydno", "H. melp"))+
  scale_color_manual(values=c("black","grey65"))+
  theme(legend.position="", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(1,0.5,0.5,0.5, "cm"))+
  theme(legend.text=element_text(size=11))+
  theme(legend.key=element_blank())+
  theme(axis.text.x=element_text(face = "bold"))+
  theme(text = element_text(family = "Times New Roman"))


####MASTER PLOT OF ALL 4 EMMS - C & M ONLY####

library(ggpubr)

#####**species & sex separate**#####

ggarrange(species, sex, species2, sex2, 
                       ncol = 2, nrow = 2, widths = c(1, 1.05),
                       font.label = list(size = 15))

#####**species:sex**#####
ggarrange(combo1, combo2, 
          labels = c("A", "B"),
          vjust=2,
          ncol = 2, nrow = 1, widths = c(1, 1),
          font.label = list(size = 15))


####ANALYSES OF FACET COUNT - C, M & F1 HYBRIDS####

#full statistical model

hist(eye$whole_count_log) 

w1 <- lm(whole_count_log~type*sex+tibia_length_log, data=eye)
plot(w1, which=1)
plot(w1, which=2)
hist(resid(w1))

drop1(w1, test="Chisq")
w2 <- update(w1,.~.-type:sex) #type:sex interaction ns
drop1(w2, test="Chisq") #all three significant

plot(allEffects(w2))
plot(w2, which = 2) 
hist(resid(w2))

Anova(w2, test.statistic="Chisq") 
#type p>0.001
#sex p>0.001
#tibia p>0.001

Anova(w1, test.statistic="Chisq") 
#type:sex p=0.244


#posthoc comparisons
#by species type
library(multcomp)
posthoc <- glht(w2, linfct=mcp(type="Tukey"))
summary(posthoc,adjusted(type="bonferroni"))


#####plot of facet count based on MAM####

######**species**####
model_means_w2 <- emmeans(object = w2, specs = ~ type, type = "response") 

#p-value comparison table
pwpm(model_means_w2, adjust="bonferroni",  diffs = F)

#effect size
eff_size(model_means_w2, sigma = sigma(w2), edf = df.residual(w2))


#add letters to each mean to indicate difference
model_means_cld_w2 <- cld(object = model_means_w2,
                          reversed=T,
                          adjust="bonferroni",
                          Letters = letters,
                          alpha = 0.05)

#set marginal means at a dataframe
facet_w2 <- as.data.frame(summary(model_means_cld_w2))

#species

#color selection for consistenty with other papers
library(RColorBrewer)
display.brewer.pal(12, 'Paired')
brewer.pal(12, 'Paired')

facet_w2$type <- factor(facet_w2$type, c("CP", "CPxMP", "MPxCP", "MP"))

ggplot(facet_w2, aes(type,emmean))+
  geom_errorbar(data=facet_w2,mapping=aes(x=type, ymin=upper.CL, ymax=lower.CL), color="grey45",
                width=0, linewidth=1, position = position_dodge(0.4))+
  geom_point(data=facet_w2, aes(type,emmean, color=type), size=4,
             position = position_dodge(0.4))+
  geom_text(aes(label = .group, y = emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=4.5) +
  labs(x="",y="e.m.m. log10 (facet count)")+
  scale_y_continuous(breaks=seq(4.11, 4.17, 0.02), limits=c(4.11, 4.175))+
  scale_x_discrete(labels=c("C", "CxM" , "MxC" ,"M"))+
  scale_color_manual(values = c("#1F78B4", "#FF7F00" ,"grey40" ,"#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(plot.margin = margin(1,0.5,0.5,0.75, "cm"))+
  #theme_update(axis.text.x=element_text(vjust=0.5, angle=45))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(axis.text.x=element_text(face = "bold"))+
  theme(text = element_text(family = "Times New Roman"))+
  ggtitle("facet count") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


######**sex**####
model_means_w2.1 <- emmeans(object = w2, specs = ~ sex, type = "response") 

#p-value comparison table
pwpm(model_means_w2.1, adjust="bonferroni",  diffs = F)

#effect size
eff_size(model_means_w2.1, sigma = sigma(w2), edf = df.residual(w2))

#add letters to each mean to indicate difference
model_means_cld_w2.1 <- cld(object = model_means_w2.1,
                            reversed = T,
                            adjust="bonferroni",
                            Letters = letters,
                            alpha = 0.05)

#set marginal means at a dataframe
facet_w2.1 <- as.data.frame(summary(model_means_cld_w2.1))

#plot marginal means
facet_w2.1$sex <- factor(facet_w2.1$sex, c("male", "female"))

#sex
ggplot(facet_w2.1, aes(sex,emmean))+
  geom_errorbar(data=facet_w2.1,mapping=aes(x=sex, ymin=upper.CL, ymax=lower.CL), color="grey45",
                width=0, linewidth=1)+
  geom_point(data=facet_w2.1, aes(sex,emmean), size=4)+
  geom_text(aes(label = .group, y = emmean),
            vjust=0.3, hjust=-0.3,size=4.5) +
  labs(x="",y="")+
  scale_y_continuous(breaks=seq(4.11, 4.17, 0.02), limits=c(4.11, 4.17),
                     position="right")+
  scale_x_discrete(labels = c("male", "female"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(legend.position=c(0.88,0.89), legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(1,0.25,0.25,0.7, "cm"))


######**species:sex**####
model_means_w2.2 <- emmeans(object = w2, specs = ~ type:sex, type = "response") 

#p-value comparison table
pwpm(model_means_w2.2, adjust="bonferroni",  diffs = F)

#add letters to each mean to indicate difference
model_means_cld_w2.2  <- cld(object = model_means_w2.2 ,
                             adjust = "bonferroni",
                             reversed=T,
                             Letters = letters,
                             alpha = 0.05)

#set marginal means at a dataframe
facet_w2.2  <- as.data.frame(summary(model_means_cld_w2.2 ))


#plot marginal means
facet_w2.2 $sex <- factor(facet_w2.2 $sex, c("male", "female"))
facet_w2.2 $type <- factor(facet_w2.2 $type, c("CP", "CPxMP", "MPxCP", "MP"))

#species:sex plot
combo3 <- ggplot(facet_w2.2 , aes(type,emmean))+
  geom_errorbar(data=facet_w2.2 ,mapping=aes(x=type, ymin=upper.CL, ymax=lower.CL, color=sex),
                width=0, linewidth=1, position = position_dodge(width=0.5))+
  geom_point(data=facet_w2.2 , aes(type,emmean, color=sex), size=3,
             position = position_dodge2(width=0.5))+
  geom_text(aes(label = .group, y = emmean),
            vjust=0.3, hjust=-0.3,size=4.5, position = position_dodge2(width=0.5)) +
  labs(x="",y="e.m.m.  log10 (facet count)", color="")+
  scale_y_continuous(breaks=seq(4.09, 4.19, 0.02), limits=c(4.09, 4.19),
                     position="left")+
  scale_x_discrete(labels=c("C", "CxM", "MxC", "M"))+
  scale_color_manual(values = c("black","grey"))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(1,0.9,0.5,0.5, "cm"))+
  theme(legend.text=element_text(size=11))+
  theme(legend.key=element_blank())+ 
  theme(axis.text.x=element_text(face = "bold"))+
  theme(text = element_text(family = "Times New Roman"))



####ANALYSES OF CORNEAL AREA - C, M, & F1 HYBRIDS#### 
hist(eye$area_log)

k1 <- lm(area_log~type*sex+tibia_length_log, data=eye)
plot(k1,which=1)
plot(k1,which=2)
hist(resid(k1))

drop1(k1, test="Chisq") #type:sex n.s.
k2 <- update(k1,.~.-type:sex)
drop1(k2, test="Chisq") #all effects significant 

plot(k2,which=2) 
hist(resid(k2))

Anova(k2, test="Chisq")
#type p>0.001
#sex p>0.001
#tibia p>0.001

plot(allEffects(k2))

#species comparison
posthoc_size2<-glht(k2, linfct=mcp(type="Tukey"))
summary(posthoc_size2, adjusted(type="bonferroni"))

Anova(k1, test="Chisq")
#type:sex p=0.44


#####plot of corneal area based on MAM####

######**species**####

model_means_k2 <- emmeans(object = k2, specs = ~ type, type = "response") 

#p-value comparison table
pwpm(model_means_k2, adjust="bonferroni",  diffs = F)

#effect size
eff_size(model_means_k2, sigma = sigma(k2), edf = df.residual(k2))

#add letters to each mean to indicate difference
model_means_cld_k2 <- cld(object = model_means_k2,
                          reversed=T,
                          adjust="bonferroni",
                          Letters = letters,
                          alpha = 0.05)

#set marginal means at a dataframe
area_k2 <- as.data.frame(summary(model_means_cld_k2))

area_k2$type <- factor(area_k2$type, c("CP", "CPxMP", "MPxCP", "MP"))

#species
ggplot(area_k2, aes(type,emmean))+
  geom_errorbar(data=area_k2,mapping=aes(x=type, ymin=upper.CL, ymax=lower.CL), color="grey45",
                width=0, linewidth=1, position = position_dodge(0.4))+
  geom_point(data=area_k2, aes(type,emmean,color=type), size=4,
             position = position_dodge(0.4))+
  geom_text(aes(label = .group, y = emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=4.5) +
  labs(x="",y="e.m.m. log10 (corneal area)")+
  scale_y_continuous(breaks=seq(0.82, 0.92, 0.02), limits=c(0.84, 0.905))+
  scale_x_discrete(labels=c("C", "CxM", "MxC","M"))+
  scale_color_manual(values = c("#1F78B4","#FF7F00","grey40","#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(1,0.5,0.5,0.75, "cm"))+
  theme(axis.text.x=element_text(face = "bold"))+
  theme(text = element_text(family = "Times New Roman"))+
  ggtitle("species effect") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))


######**sex**####
model_means2_k2.1 <- emmeans(object = k2, specs = ~ sex, type = "response") 

#p-value comparison table
pwpm(model_means2_k2.1, adjust="bonferroni",  diffs = F)

#effect size
eff_size(model_means2_k2.1, sigma = sigma(k2), edf = df.residual(k2))

#add letters to each mean to indicate difference
model_means_cld2_k2.1 <- cld(object = model_means2_k2.1,
                             reversed = T,
                             adjust="bonferroni",
                             Letters = letters,
                             alpha = 0.05)

#set marginal means at a dataframe
area_k2.1 <- as.data.frame(summary(model_means_cld2_k2.1))

#plot marginal means
area_k2.1$sex <- factor(area_k2.1$sex, c("male", "female"))

#sex
ggplot(area_k2.1, aes(sex,emmean))+
  geom_errorbar(data=area_k2.1,mapping=aes(x=sex, ymin=upper.CL, ymax=lower.CL), color="grey45",
                width=0, linewidth=1)+
  geom_point(data=area_k2.1, aes(sex,emmean), size=4)+
  geom_text(aes(label = .group, y = emmean),
            vjust=0.3, hjust=-0.3,size=4.5) +
  labs(x="",y="")+
  scale_y_continuous(breaks=seq(0.84, 0.90, 0.02), limits=c(0.84, 0.90),
   position="right")+
  scale_x_discrete(labels = c("male", "female"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(legend.position=c(0.88,0.89), legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(1,0.25,0.25,0.7, "cm"))


######**species:sex**####
model_means3_k2.2 <- emmeans(object = k2, specs = ~ type:sex, type = "response") 

#p-value comparison table
pwpm(model_means3_k2.2, adjust="bonferroni",  diffs = F)

#add letters to each mean to indicate difference
model_means_cld3_k2.2 <- cld(object = model_means3_k2.2,
                             adjust = "bonferroni",
                             reversed = T,
                             Letters = letters,
                             alpha = 0.05)

#set marginal means at a dataframe
area_k2.2 <- as.data.frame(summary(model_means_cld3_k2.2))


#plot marginal means
area_k2.2$sex <- factor(area_k2.2$sex, c("male", "female"))

#species:sex plot
combo4<-ggplot(area_k2.2, aes(type,emmean))+
  geom_errorbar(data=area_k2.2,mapping=aes(x=type, ymin=upper.CL, ymax=lower.CL, color=sex),
                width=0, linewidth=1, position = position_dodge(width=0.5))+
  geom_point(data=area_k2.2, aes(type,emmean, color=sex), size=3,
             position = position_dodge2(width=0.5))+
  geom_text(aes(label = .group, y = emmean),
            vjust=0.3, hjust=-0.3,size=4.5, position = position_dodge2(width=0.5)) +
  labs(x="",y="e.m.m.  log10 (corneal area)", color="")+
  scale_y_continuous(breaks=seq(0.82, 0.92, 0.02), limits=c(0.82, 0.92))+
  scale_x_discrete(labels=c("C", "CxM", "MxC", "M"))+
  scale_color_manual(values=c("black","grey65"))+
  theme(legend.position="", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(1,0.5,0.5,0.5, "cm"))+
  theme(legend.text=element_text(size=11))+
  theme(legend.key=element_blank())+ 
  theme(axis.text.x=element_text(face = "bold"))+
  theme(text = element_text(family = "Times New Roman"))

######**combo species:sex figure - all species######
ggarrange(combo3, combo4, 
          labels = c("A", "B"),
          vjust=2,
          ncol = 2, nrow = 1, widths = c(1, 1),
          font.label = list(size = 15))


####DO SPECIES DIFFER IN BODY SIZE?####

#pure species
hist(eye2$tibia_length_log) 
size <- lm(tibia_length_log~type+sex, data=eye2)
plot(size, which=1)
plot(size, which=2)
hist(resid(size))

drop1(size, test="Chisq") #sex n.s.
size2 <- update(size,.~.-sex)
drop1(size2, test="Chisq")

Anova(size2, test="Chisq")
#species p<0.001
plot(allEffects(size2))

Anova(size, test="Chisq")
#sex p=0.788

#with F1 hybrids
hist(eye$tibia_length_log) 
size3 <- lm(tibia_length_log~type+sex, data=eye)
plot(size3, which=1)
plot(size3, which=2)
hist(resid(size3))

drop1(size3, test="Chisq") #sex n.s.
size4 <- update(size3,.~.-sex)
drop1(size4, test="Chisq")

Anova(size4, test="Chisq")
#type p>0.001
plot(allEffects(size4))

Anova(size3, test="Chisq")
#sex p=0.44

#species comparison
posthoc_size<-glht(size4, linfct=mcp(type="Tukey"))
summary(posthoc_size, adjusted(type="bonferroni"))


######**Body size emm plot**######
model_means_size <- emmeans(object = size4, specs = ~ type, type = "response") 

#p-value comparison table
pwpm(model_means_size, adjust="bonferroni",  diffs = F)

#effect size
eff_size(model_means_size, sigma = sigma(size4), edf = df.residual(size4))

#add letters to each mean to indicate difference
model_means_cld_size <- cld(object = model_means_size,
                            reversed=T,
                            adjust="bonferroni",
                            Letters = letters,
                            alpha = 0.05)

#set marginal means at a dataframe
sizes <- as.data.frame(summary(model_means_cld_size))

#size
sizes$type <- factor(sizes$type, c("CP", "CPxMP", "MPxCP", "MP"))

ggplot(sizes, aes(type,emmean))+
  geom_errorbar(data=sizes,mapping=aes(x=type, ymin=upper.CL, ymax=lower.CL), color="grey45",
                width=0, linewidth=1, position = position_dodge(0.4))+
  geom_point(data=sizes, aes(type,emmean,color=type), size=4,
             position = position_dodge(0.4))+
  geom_text(aes(label = .group, y = emmean), position = position_dodge(0.4),
            vjust=0.3, hjust=-0.2, size=4.5) +
  labs(x="",y="e.m.m. log10 (tibia length)")+
  scale_y_continuous(breaks=seq(0.68, 0.76, 0.02), limits=c(0.68, 0.76))+
  scale_x_discrete(labels=c("C", "CxM", "MxC", "M"))+
  scale_color_manual(values = c("#1F78B4","#FF7F00", "grey40", "#FDBF6F"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(1,0.5,0.5,0.75, "cm"))+
  theme(axis.text.x=element_text(face = "bold"))+
  theme(text = element_text(family = "Times New Roman"))+
  ggtitle("body size") + theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))



####FDR P-VALUE ADJUSTMENTS FROM ALL LMs####

pvalues <- c(0.008775,
             2.06E-06,
             2.77E-05,
             7.19E-01,
             0.004985,
             3.84E-05,
             8.61E-14,
             0.102394,
             0.0001924,
             3.81E-09,
             2.64E-07,
             0.2441719,
             4.18E-06,
             1.14E-07,
             2.20E-16,
             0.4434,
             1.15E-10,
             7.89E-01,
             6.54E-13,
             0.4432,
             2.27E-16,
             0.1203373)

adjusted.p <- p.adjust(pvalues, method = "fdr", n = length(pvalues))
adjusted.p


####COMPARISON TO PREVIOUS PAPER - SEYMOURE ET AL. (2015)####

#means per species & sex from this study
#pure species only
eye2 %>%
  group_by(type:sex) %>%
  summarise_at(vars(area), list(name = mean))

#import data extracted from Seymoure et al. (2015)
method_comp <- read.csv("~/Documents/LMU/Eye morphology/selection manuscript/method_comp.csv")

eye2$sex<- factor(eye2$sex, c("male", "female"))
method_comp$sex<- factor(method_comp$sex, c("male", "female"))

facet_comp <- ggplot(eye2, aes(type,whole_count))+ facet_grid(~sex)+
  geom_point(data=eye2, aes(type,whole_count), size=2.5,shape=1,color="grey40",
             position = position_jitter(0.085))+
  geom_point(data=method_comp, aes(type,whole_count), size=3, shape=17,
             position = position_jitter(0.085))+
  labs(x="",y="facet count")+
  scale_y_continuous(breaks=seq(11000, 17000, 1000), limits=c(11000, 17000))+
  scale_x_discrete(labels=c("H. cydno","H. melp"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(0.4,0.4,0.4,0.75, "cm"))+
  theme(axis.text.x=element_text(face = "bold"))+
  theme(strip.background = element_rect(fill = "transparent"))+ 
  theme(strip.text = element_text(size = 14, face = "bold"))+
  theme(text = element_text(family = "Times New Roman"))


area_comp <- ggplot(eye2, aes(type,area))+ facet_grid(~sex)+
  geom_point(data=eye2, aes(type,area), size=2.5,shape=1,color="grey40",
             position = position_jitter(0.085))+
  geom_point(data=method_comp, aes(type,area), size=3, shape=17,
             position = position_jitter(0.085))+
  labs(x="",y=bquote(bold('corneal area'~(mm^2))))+
  scale_y_continuous(breaks=seq(5.5, 9.5, 1.00), limits=c(5.4, 9.6))+
  scale_x_discrete(labels=c("H. cydno","H. melp"))+
  theme_update(axis.title.y=element_text(size=14,face="bold",hjust=0.5,vjust=2.5,angle=90))+
  theme(legend.position="none", legend.key=element_blank(),legend.title=element_blank(),
        legend.text=element_text(size=10))+
  theme(plot.margin = margin(0.4,0.4,0.4,1.08, "cm"))+
  theme(axis.text.x=element_text(face = "bold"))+
  theme(strip.background = element_rect(fill = "transparent"))+ 
  theme(strip.text = element_text(size = 14, face = "bold"))+
  theme(text = element_text(family = "Times New Roman"))


#combine plots
library(ggpubr)

ggarrange(facet_comp, area_comp, 
          labels = c("A", "B"),
          vjust=1.23, hjust=-0.25,
          ncol = 1, nrow = 2, widths = c(1, 1))


####SELECTION TESTS#####

library("Pstat")

#copy type (species) column and move to the first position
eye2.1 <- cbind(type_1 = eye2$type, eye2[, !(names(eye2) %in% "type")])

#####1) raw measurements####
raw_0.33 <- Pst(eye2.1, ci=1, csh=0.33, va=c("whole_count", "area"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

raw_0.33$level <- c(0.33,0.33)


raw_0.67 <- Pst(eye2.1, ci=1, csh=0.67, va=c("whole_count", "area"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

raw_0.67$level <- c(0.67,0.67)


raw_1 <- Pst(eye2.1, ci=1, csh=1.0, va=c("whole_count", "area"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

raw_1$level <- c(1,1)


raw_1.33 <- Pst(eye2.1, ci=1, csh=1.33, va=c("whole_count", "area"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

raw_1.33$level <- c(1.33,1.33)


raw_1.5 <- Pst(eye2.1, ci=1, csh=1.5, va=c("whole_count", "area"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)

raw_1.5$level <- c(1.5,1.5)


raw_2 <- Pst(eye2.1, ci=1, csh=2, va=c("whole_count", "area"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

raw_2$level <- c(2,2)


raw_3 <- Pst(eye2.1, ci=1, csh=3, va=c("whole_count", "area"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

raw_3$level <- c(3,3)


raw_4 <- Pst(eye2.1, ci=1, csh=4, va=c("whole_count", "area"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

raw_4$level <- c(4,4)


raw_values <- rbind(raw_0.33, raw_0.67, raw_1, raw_1.33, raw_1.5, raw_2, raw_3, raw_4)


#####2) log-transformed measurements####
log_0.33 <- Pst(eye2.1, ci=1, csh=0.33, va=c("whole_count_log", "area_log"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

log_0.33$level <- c(0.33,0.33)


log_0.67 <- Pst(eye2.1, ci=1, csh=0.67, va=c("whole_count_log", "area_log"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

log_0.67$level <- c(0.67,0.67)


log_1 <- Pst(eye2.1, ci=1, csh=1, va=c("whole_count_log", "area_log"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

log_1$level <- c(1,1)


log_1.33 <- Pst(eye2.1, ci=1, csh=1.33, va=c("whole_count_log", "area_log"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

log_1.33$level <- c(1.33,1.33)


log_1.5 <- Pst(eye2.1, ci=1, csh=1.5, va=c("whole_count_log", "area_log"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)

log_1.5$level <- c(1.5,1.5)


log_2 <- Pst(eye2.1, ci=1, csh=2, va=c("whole_count_log", "area_log"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

log_2$level <- c(2,2)


log_3 <- Pst(eye2.1, ci=1, csh=3, va=c("whole_count_log", "area_log"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

log_3$level <- c(3,3)


log_4 <- Pst(eye2.1, ci=1, csh=4, va=c("whole_count_log", "area_log"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

log_4$level <- c(4,4)


log_values <- rbind(log_0.33, log_0.67, log_1, log_1.33, log_1.5, log_2, log_3, log_4)



#####3) allometric-scaled####

library(allomr)

#maybe different per species & sex? need to do separate?

######CP males#####

#subset the data
eye2.1.cp.m <- eye2.1[which(eye2.1$type=="CP"),]
eye2.1.cp.m <- eye2.1.cp.m[which(eye2.1.cp.m$sex=="male"),]

#apply allometric scaling
count_al_cpm <- as.data.frame(allomr(eye2.1.cp.m$tibia_length, eye2.1.cp.m$whole_count))
area_al_cpm <- as.data.frame(allomr(eye2.1.cp.m$tibia_length, eye2.1.cp.m$area))

#add scaled values back to dataset
eye2.1.cp.m2 <- cbind(eye2.1.cp.m, count_al_cpm$Yx)
colnames(eye2.1.cp.m2)[colnames(eye2.1.cp.m2) == "count_al_cpm$Yx"] <- "whole_count_al"

eye2.1.cp.m2 <- cbind(eye2.1.cp.m2, area_al_cpm$Yx)
colnames(eye2.1.cp.m2)[colnames(eye2.1.cp.m2) == "area_al_cpm$Yx"] <- "area_al"


######MP males####

#subset the data
eye2.1.mp.m <- eye2.1[which(eye2.1$type=="MP"),]
eye2.1.mp.m <- eye2.1.mp.m[which(eye2.1.mp.m$sex=="male"),]

#apply allometric scaling
count_al_mpm <- as.data.frame(allomr(eye2.1.mp.m$tibia_length, eye2.1.mp.m$whole_count))
area_al_mpm <- as.data.frame(allomr(eye2.1.mp.m$tibia_length, eye2.1.mp.m$area))

#add scaled values back to dataset
eye2.1.mp.m2 <- cbind(eye2.1.mp.m, count_al_mpm$Yx)
colnames(eye2.1.mp.m2)[colnames(eye2.1.mp.m2) == "count_al_mpm$Yx"] <- "whole_count_al"

eye2.1.mp.m2 <- cbind(eye2.1.mp.m2, area_al_mpm$Yx)
colnames(eye2.1.mp.m2)[colnames(eye2.1.mp.m2) == "area_al_mpm$Yx"] <- "area_al"


######CP females####

#subset the data
eye2.1.cp.f <- eye2.1[which(eye2.1$type=="CP"),]
eye2.1.cp.f <- eye2.1.cp.f[which(eye2.1.cp.f$sex=="female"),]

#apply allometic scaling
count_al_cpf <- as.data.frame(allomr(eye2.1.cp.f$tibia_length, eye2.1.cp.f$whole_count))
area_al_cpf <- as.data.frame(allomr(eye2.1.cp.f$tibia_length, eye2.1.cp.f$area))

#add scaled values back to dataset
eye2.1.cp.f2 <- cbind(eye2.1.cp.f, count_al_cpf$Yx)
colnames(eye2.1.cp.f2)[colnames(eye2.1.cp.f2) == "count_al_cpf$Yx"] <- "whole_count_al"

eye2.1.cp.f2 <- cbind(eye2.1.cp.f2, area_al_cpf$Yx)
colnames(eye2.1.cp.f2)[colnames(eye2.1.cp.f2) == "area_al_cpf$Yx"] <- "area_al"


######MP females####

#subset the data
eye2.1.mp.f <- eye2.1[which(eye2.1$type=="MP"),]
eye2.1.mp.f <- eye2.1.mp.f[which(eye2.1.mp.f$sex=="female"),]

#apply allometric scaling
count_al_mpf <- as.data.frame(allomr(eye2.1.mp.f$tibia_length, eye2.1.mp.f$whole_count))
area_al_mpf <- as.data.frame(allomr(eye2.1.mp.f$tibia_length, eye2.1.mp.f$area))

#add scaled values back to dataset
eye2.1.mp.f2 <- cbind(eye2.1.mp.f, count_al_mpf$Yx)
colnames(eye2.1.mp.f2)[colnames(eye2.1.mp.f2) == "count_al_mpf$Yx"] <- "whole_count_al"

eye2.1.mp.f2 <- cbind(eye2.1.mp.f2, area_al_mpf$Yx)
colnames(eye2.1.mp.f2)[colnames(eye2.1.mp.f2) == "area_al_mpf$Yx"] <- "area_al"


#combine datasets for analyses

eye2.2 <- (rbind(eye2.1.cp.m2, eye2.1.mp.m2, eye2.1.cp.f2, eye2.1.mp.f2))


#and now with allomertic-scaled measurements
all_0.33 <- Pst(eye2.2, ci=1, csh=0.33, va=c("whole_count_al", "area_al"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)
all_0.33$level <- c(0.33,0.33)


all_0.67 <- Pst(eye2.2, ci=1, csh=0.67, va=c("whole_count_al", "area_al"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)
all_0.67$level <- c(0.67,0.67)


all_1 <- Pst(eye2.2, ci=1, csh=1, va=c("whole_count_al", "area_al"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)
all_1$level <- c(1.00,1.00)


all_1.33 <- Pst(eye2.2, ci=1, csh=1.33, va=c("whole_count_al", "area_al"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)
all_1.33$level <- c(1.33,1.33)


all_1.5 <- Pst(eye2.2, ci=1, csh=1.5, va=c("whole_count_al", "area_al"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)
all_1.5$level <- c(1.50,1.50)


all_2 <- Pst(eye2.2, ci=1, csh=2, va=c("whole_count_al", "area_al"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)
all_2$level <- c(2.00,2.00)


all_3 <- Pst(eye2.2, ci=1, csh=3, va=c("whole_count_al", "area_al"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)
all_3$level <- c(3.00,3.00)


all_4 <- Pst(eye2.2, ci=1, csh=4, va=c("whole_count_al", "area_al"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)
all_4$level <- c(4.00,4.00)


all_values <- rbind(all_0.33, all_0.67, all_1, all_1.33, all_1.5, all_2, all_3, all_4)


######combine pst values into one data frame####

pst_values <- rbind(raw_values, log_values, all_values)


####Fst values from Martin et al. 2013####

fst <- read.csv("~/Documents/LMU/Eye morphology/selection manuscript/fst_martin2013.csv")

ggplot(fst,aes(chi_rosFst)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "grey95") +
  geom_density(fill="grey85", alpha=0.6)+
  labs(x= "Fst H. cydno - H. melpomene", y="density")


######**Plotted with Pst values**####

#with only facet count
pst_values1 <- pst_values[which(pst_values$Quant_Varia=="whole_count_al"),]

quantile(fst$chi_rosFst, 0.95, na.rm=T)

pst_values1$Pst_Values

x.expression <- expression(bold(F[ST] ~ "H. cydno -" ~ "H. melpomene"))

count<- ggplot(fst,aes(chi_rosFst)) + 
  geom_histogram(aes(y = ..density..),
                 colour = "grey70", fill = "grey90") +
  geom_density(colour="grey70",fill="lightgrey", alpha=0.4)+
  geom_vline(aes(xintercept=0.8788755, color="0.33"), size=1.3, linetype="dashed",key_glyph = "rect")+
  geom_vline(aes(xintercept=0.9364345, color="0.67"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9564986, color="1.00"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9669353, color="1.33"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9705724, color="1.50"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9777657, color="2.00"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9850665, color="3.00"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9887579, color="4.00"), size=1.3, linetype="dashed")+
  scale_color_manual(
    name = expression(bold(c/h^2)), 
    values = c("0.33" = "yellow", "0.67" = "orange", "1.00" = "darkorange", "1.33" = "red", "1.50" = "red2",
               "2.00" = "red3", "3.00" = "red4", "4.00" = "black"))+
  geom_vline(xintercept=0.584, linetype="dotted", color="black", size=0.75)+
  labs(x= x.expression, y="density")+
  scale_x_continuous(limits = c(0,1.0))+
  scale_y_continuous(limits = c(0,3.4), expand = c(0, 0.02))+
  theme(legend.key=element_blank(), legend.title=element_text(size=12),
        legend.title.align=0.5,legend.text=element_text(size=12))+
  theme(plot.margin = margin(0.5,0,0.5,0.5, "cm"))+
  ggtitle("facet count") + theme(plot.title = element_text(hjust = 0.5, size=15,face="bold"))+
  theme(text = element_text(family = "Times New Roman"))


#with only area
pst_values2 <- pst_values[which(pst_values$Quant_Varia=="area_al"),]

quantile(fst$chi_rosFst, 0.95, na.rm=T)

pst_values2$Pst_Values

area <- ggplot(fst,aes(chi_rosFst)) + 
  geom_histogram(aes(y = ..density..),
                 colour = "grey70", fill = "grey90") +
  geom_density(colour="grey70",fill="lightgrey", alpha=0.4)+
  geom_vline(aes(xintercept=0.9608007, color="0.33"), size=1.3, linetype="dashed",key_glyph = "rect")+
  geom_vline(aes(xintercept=0.9803010, color="0.67"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9867153, color="1.00"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9899785, color="1.33"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9911042, color="1.50"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9933132, color="2.00"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9955322, color="3.00"), size=1.3, linetype="dashed")+
  geom_vline(aes(xintercept=0.9966454, color="4.00"), size=1.3, linetype="dashed")+
  scale_color_manual(
    name = expression(bold(c/h^2)), 
    values = c("0.33" = "yellow", "0.67" = "orange", "1.00" = "darkorange", "1.33" = "red", "1.50" = "red2",
               "2.00" = "red3", "3.00" = "red4", "4.00" = "black"))+
  geom_vline(xintercept=0.584, linetype="dotted", color="black", size=0.75)+
  labs(x= x.expression, y="density")+
  scale_x_continuous(limits = c(0,1.0))+
  scale_y_continuous(limits = c(0,3.4), expand = c(0, 0.02),position = "right")+
  theme(legend.key=element_blank(), legend.title=element_text(size=12),
        legend.title.align=0.5,legend.text=element_text(size=12))+
  theme(legend.position='none')+
  theme(plot.margin = margin(0.5,0.5,0.5,0.1, "cm"))+
  ggtitle("corneal area") + theme(plot.title = element_text(hjust = 0.5, size=15,face="bold"))+
  theme(text = element_text(family = "Times New Roman"))


#combine plots
library(ggpubr)

ggarrange(count, area, 
          #labels = c("A", "B"),
          #vjust=1.23, hjust=0.8,
          ncol = 2, nrow = 1, widths = c(1.2, 1))

######calculate p-values for each Pst####

#use only cydno vs. melp fst values
fst2 <- as.data.frame(fst$chi_rosFst)

#remove NAs from the Fst values
fst2.1 <- complete.cases(fst2)
fst3 <- fst2[fst2.1, ]


#calculate the proportion of the Fst distribution above each Pst value
proportions <- lapply(pst_values$Pst_Values, function(x) sum(fst3 > x) / length(fst3))

#combine the x_values and proportions into a data frame
p.values <- data.frame(x = pst_values$Pst_Values, p_value = unlist(proportions))

pst_values <- cbind(pst_values, p.values)

#double-check that Pst values are aligned in new dataset
pst_values$Pst_Values - pst_values$x #all zeros, everything is correct

#export to make a table
write.csv(pst_values,file='~/Documents/LMU/Eye morphology/selection manuscript/pst_table.csv')



####SELECTION TESTS PER SEX SEPARATELY#####

library("Pstat")


eye3m <- eye2.1[which(eye2.1$sex=="male"),]

#####1) raw measurements_male####
raw_0.33_m <- Pst(eye3m, ci=1, csh=0.33, va=c("whole_count", "area"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)

raw_0.33_m$level <- c(0.33,0.33)


raw_0.67_m <- Pst(eye3m, ci=1, csh=0.67, va=c("whole_count", "area"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

raw_0.67_m$level <- c(0.67,0.67)


raw_1_m <- Pst(eye3m, ci=1, csh=1.0, va=c("whole_count", "area"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

raw_1_m$level <- c(1,1)


raw_1.33_m <- Pst(eye3m, ci=1, csh=1.33, va=c("whole_count", "area"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

raw_1.33_m$level <- c(1.33,1.33)


raw_1.5_m <- Pst(eye3m, ci=1, csh=1.5, va=c("whole_count", "area"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)

raw_1.5_m$level <- c(1.5,1.5)


raw_2_m <- Pst(eye3m, ci=1, csh=2, va=c("whole_count", "area"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

raw_2_m$level <- c(2,2)


raw_3_m <- Pst(eye3m, ci=1, csh=3, va=c("whole_count", "area"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

raw_3_m$level <- c(3,3)


raw_4_m <- Pst(eye3m, ci=1, csh=4, va=c("whole_count", "area"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

raw_4_m$level <- c(4,4)


raw_values_m <- rbind(raw_0.33_m, raw_0.67_m, raw_1_m, 
                      raw_1.33_m, raw_1.5_m, raw_2_m, raw_3_m, raw_4_m)


#####2) log-transformed measurements_male####
log_0.33_m <- Pst(eye3m, ci=1, csh=0.33, va=c("whole_count_log", "area_log"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

log_0.33_m$level <- c(0.33,0.33)


log_0.67_m <- Pst(eye3m, ci=1, csh=0.67, va=c("whole_count_log", "area_log"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

log_0.67_m$level <- c(0.67,0.67)


log_1_m <- Pst(eye3m, ci=1, csh=1, va=c("whole_count_log", "area_log"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

log_1_m$level <- c(1,1)


log_1.33_m <- Pst(eye3m, ci=1, csh=1.33, va=c("whole_count_log", "area_log"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

log_1.33_m$level <- c(1.33,1.33)


log_1.5_m <- Pst(eye3m, ci=1, csh=1.5, va=c("whole_count_log", "area_log"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)

log_1.5_m$level <- c(1.5,1.5)


log_2_m <- Pst(eye3m, ci=1, csh=2, va=c("whole_count_log", "area_log"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

log_2_m$level <- c(2,2)


log_3_m <- Pst(eye3m, ci=1, csh=3, va=c("whole_count_log", "area_log"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

log_3_m$level <- c(3,3)


log_4_m <- Pst(eye3m, ci=1, csh=4, va=c("whole_count_log", "area_log"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

log_4_m$level <- c(4,4)


log_values_m <- rbind(log_0.33_m, log_0.67_m, log_1_m, log_1.33_m, 
                      log_1.5_m, log_2_m, log_3_m, log_4_m)


#####3) allometric-scaled_male####

eye3m_al <- (rbind(eye2.1.cp.m2, eye2.1.mp.m2))

#with allomertic-scaled measurements
all_0.33_m <- Pst(eye3m_al, ci=1, csh=0.33, va=c("whole_count_al", "area_al"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

all_0.33_m$level <- c(0.33,0.33)


all_0.67_m <- Pst(eye3m_al, ci=1, csh=0.67, va=c("whole_count_al", "area_al"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

all_0.67_m$level <- c(0.67,0.67)


all_1_m <- Pst(eye3m_al, ci=1, csh=1, va=c("whole_count_al", "area_al"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

all_1_m$level <- c(1,1)


all_1.33_m <- Pst(eye3m_al, ci=1, csh=1.33, va=c("whole_count_al", "area_al"), boot=1000,
                Pw = c("CP","MP"), pe = 0.95)

all_1.33_m$level <- c(1.33,1.33)


all_1.5_m <- Pst(eye3m_al, ci=1, csh=1.5, va=c("whole_count_al", "area_al"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)

all_1.5_m$level <- c(1.5,1.5)


all_2_m <- Pst(eye3m_al, ci=1, csh=2, va=c("whole_count_al", "area_al"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

all_2_m$level <- c(2,2)


all_3_m <- Pst(eye3m_al, ci=1, csh=3, va=c("whole_count_al", "area_al"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

all_3_m$level <- c(3,3)


all_4_m <- Pst(eye3m_al, ci=1, csh=4, va=c("whole_count_al", "area_al"), boot=1000,
             Pw = c("CP","MP"), pe = 0.95)

all_4_m$level <- c(4,4)


all_values_m <- rbind(all_0.33_m, all_0.67_m, all_1_m, 
                      all_1.33_m, all_1.5_m, all_2_m, all_3_m, all_4_m)


#####combine pst values into one male data frame with p-values####

pst_values_m <- rbind(raw_values_m, log_values_m, all_values_m)

#calculate p-values 

#use only cydno vs. melp fst values, fst3

#calculate the proportion of the Fst distribution above each Pst value
proportions_m <- lapply(pst_values_m$Pst_Values, function(x) sum(fst3 > x) / length(fst3))

# Combine the x_values and proportions into a data frame
p.values_m <- data.frame(x = pst_values_m$Pst_Values, p_value = unlist(proportions_m))

pst_values_m <- cbind(pst_values_m, p.values_m)

#double-check that Pst values are aligned in new dataset
pst_values_m$Pst_Values - pst_values_m$x #all zeros, everything is correct

#export to make a table
write.csv(pst_values_m,file='~/Documents/LMU/Eye morphology/selection manuscript/pst_table_male.csv')


#####1) raw measurements_female####

eye3f <- eye2.1[which(eye2.1$sex=="female"),]


raw_0.33_f <- Pst(eye3f, ci=1, csh=0.33, va=c("whole_count", "area"), boot=1000,
                  Pw = c("CP","MP"), pe = 0.95)
raw_0.33_f$level <- c(0.33,0.33)


raw_0.67_f <- Pst(eye3f, ci=1, csh=0.67, va=c("whole_count", "area"), boot=1000,
                  Pw = c("CP","MP"), pe = 0.95)
raw_0.67_f$level <- c(0.67,0.67)


raw_1.0_f <- Pst(eye3f, ci=1, csh=1.0, va=c("whole_count", "area"), boot=1000,
                  Pw = c("CP","MP"), pe = 0.95)
raw_1.0_f$level <- c(1.0,1.0)


raw_1.33_f <- Pst(eye3f, ci=1, csh=1.33, va=c("whole_count", "area"), boot=1000,
                 Pw = c("CP","MP"), pe = 0.95)
raw_1.33_f$level <- c(1.33,1.33)


raw_1.5_f <- Pst(eye3f, ci=1, csh=1.5, va=c("whole_count", "area"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)
raw_1.5_f$level <- c(1.5,1.5)


raw_2.0_f <- Pst(eye3f, ci=1, csh=2.0, va=c("whole_count", "area"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)
raw_2.0_f$level <- c(2.0,2.0)


raw_3.0_f <- Pst(eye3f, ci=1, csh=3.0, va=c("whole_count", "area"), boot=1000,
                 Pw = c("CP","MP"), pe = 0.95)
raw_3.0_f$level <- c(3.0,3.0)


raw_4.0_f <- Pst(eye3f, ci=1, csh=4.0, va=c("whole_count", "area"), boot=1000,
                 Pw = c("CP","MP"), pe = 0.95)
raw_4.0_f$level <- c(4.0,4.0)


raw_values_f <- rbind(raw_0.33_f, raw_0.67_f, raw_1.0_f, raw_1.33_f, raw_1.5_f,
                      raw_2.0_f, raw_3.0_f, raw_4.0_f)


#####2) log-transformed measurements_female####
log_0.33_f <- Pst(eye3f, ci=1, csh=0.33, va=c("whole_count_log", "area_log"), boot=1000,
                  Pw = c("CP","MP"), pe = 0.95)
log_0.33_f$level <- c(0.33,0.33)


log_0.67_f <- Pst(eye3f, ci=1, csh=0.67, va=c("whole_count_log", "area_log"), boot=1000,
                  Pw = c("CP","MP"), pe = 0.95)
log_0.67_f$level <- c(0.67,0.67)


log_1.0_f <- Pst(eye3f, ci=1, csh=1.0, va=c("whole_count_log", "area_log"), boot=1000,
                  Pw = c("CP","MP"), pe = 0.95)
log_1.0_f$level <- c(1.0,1.0)


log_1.33_f <- Pst(eye3f, ci=1, csh=1.33, va=c("whole_count_log", "area_log"), boot=1000,
                 Pw = c("CP","MP"), pe = 0.95)
log_1.33_f$level <- c(1.33,1.33)


log_1.5_f <- Pst(eye3f, ci=1, csh=1.5, va=c("whole_count_log", "area_log"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)
log_1.5_f$level <- c(1.5,1.5)


log_2.0_f <- Pst(eye3f, ci=1, csh=2, va=c("whole_count_log", "area_log"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)
log_2.0_f$level <- c(2.0,2.0)


log_3.0_f <- Pst(eye3f, ci=1, csh=3, va=c("whole_count_log", "area_log"), boot=1000,
                 Pw = c("CP","MP"), pe = 0.95)
log_3.0_f$level <- c(3.0,3.0)


log_4.0_f <- Pst(eye3f, ci=1, csh=4, va=c("whole_count_log", "area_log"), boot=1000,
                 Pw = c("CP","MP"), pe = 0.95)
log_4.0_f$level <- c(4.0,4.0)


log_values_f <- rbind(log_0.33_f, log_0.67_f, log_1.0_f, log_1.33_f, log_1.5_f,
                      log_2.0_f, log_3.0_f, log_4.0_f)


#####3) allometric-scaled_female####

eye3f_al <- (rbind(eye2.1.cp.f2, eye2.1.mp.f2))

#with allomertic-scaled measurements
all_0.33_f <- Pst(eye3f_al, ci=1, csh=0.33, va=c("whole_count_al", "area_al"), boot=1000,
                  Pw = c("CP","MP"), pe = 0.95)
all_0.33_f$level <- c(0.33,0.33)


all_0.67_f <- Pst(eye3f_al, ci=1, csh=0.67, va=c("whole_count_al", "area_al"), boot=1000,
                  Pw = c("CP","MP"), pe = 0.95)
all_0.67_f$level <- c(0.67,0.67)


all_1.0_f <- Pst(eye3f_al, ci=1, csh=1.0, va=c("whole_count_al", "area_al"), boot=1000,
                  Pw = c("CP","MP"), pe = 0.95)
all_1.0_f$level <- c(1.0,1.0)


all_1.33_f <- Pst(eye3f_al, ci=1, csh=1.33, va=c("whole_count_al", "area_al"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)
all_1.33_f$level <- c(1.33,1.33)


all_1.5_f <- Pst(eye3f_al, ci=1, csh=1.5, va=c("whole_count_al", "area_al"), boot=1000,
                  Pw = c("CP","MP"), pe = 0.95)
all_1.5_f$level <- c(1.5,1.5)


all_2.0_f <- Pst(eye3f_al, ci=1, csh=2, va=c("whole_count_al", "area_al"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)
all_2.0_f$level <- c(2.0,2.0)


all_3.0_f <- Pst(eye3f_al, ci=1, csh=3, va=c("whole_count_al", "area_al"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)
all_3.0_f$level <- c(3.0,3.0)


all_4.0_f <- Pst(eye3f_al, ci=1, csh=4, va=c("whole_count_al", "area_al"), boot=1000,
               Pw = c("CP","MP"), pe = 0.95)
all_4.0_f$level <- c(4.0,4.0)


all_values_f <- rbind(all_0.33_f, all_0.67_f, all_1.0_f, all_1.33_f, all_1.5_f,
                      all_2.0_f, all_3.0_f, all_4.0_f)

#####combine pst values into one female data frame with p-values####

pst_values_f <- rbind(raw_values_f, log_values_f, all_values_f)

#calculate p-values 
#use only cydno vs. melp fst values, fst3

#calculate the proportion of the Fst distribution above each Pst value
proportions_f <- lapply(pst_values_f$Pst_Values, function(x) sum(fst3 > x) / length(fst3))

# Combine the x_values and proportions into a data frame
p.values_f <- data.frame(x = pst_values_f$Pst_Values, p_value = unlist(proportions_f))

pst_values_f <- cbind(pst_values_f, p.values_f)

#double-check that Pst values are aligned in new dataset
pst_values_f$Pst_Values - pst_values_f$x #all zeros, everything is correct

#export to make a table
write.csv(pst_values_f,file='~/Documents/LMU/Eye morphology/selection manuscript/pst_table_female.csv')


####ALLOMETRY#####
library(smatr)

#####1) as an interactive term####

#smatr cannot handle 2 grouping terms, so make an interaction
eye2$inter <- interaction(eye2$type, eye2$sex)

eye2$inter<- factor(eye2$inter, c("CP.male", "CP.female", 
                                "MP.male", "MP.female"))


######testing differences in slope (beta)###### 


#facet count
a1 <- sma(whole_count_log~tibia_length_log*inter, robust=F, data=eye2, multcomp = T, multcompmethod = "adjust")
multcompmatrix(a1, sort=T)
summary(a1)
plot(a1)


########**plot facet count major axis regressions**######
dev.off()
par(mgp=c(3.5,0.6,0))
par(mar = c(5, 6, 1.25, 1.25))
par(lwd=1)

plot(a1, col=c("#1F78B4", "#1F78B4", "#FDBF6F", "#FDBF6F"),
     pch=c(17,2,16,1),
     lwd=1,
     cex=0.9,
     lty=c(1,2,1,2),
     xlim=c(0.635,0.80), ylim=c(4.03,4.25),
     las=1,cex.axis = 1.3,
     xlab = "",
     ylab = "",
     tck=-0.015,
     family = "Times New Roman")

#title(xlab = "log10 (hind tibia length)", line = 2,cex.lab = 1.2,font.lab=2,family = "Times New Roman")           
title(ylab = "log10 (facet count)", line = 3,cex.lab = 1.2,font.lab=2,family = "Times New Roman") 
title("facet count", adj = 0.5, line = 0.6, family = "Times New Roman", cex.main=1.35)


legend("topleft", c("C", "M"),
       col = c("#1F78B4", "#FDBF6F"), pch=c(17,16),
       bty="n")


#area
a2 <- sma(area_log~tibia_length_log*inter, data=eye2, multcomp = T, multcompmethod = "adjust")
multcompmatrix(a2, sort=T)
summary(a2)
plot(a2)

########**plot corneal area major axis regressions**######

plot(a2, col=c("#1F78B4", "#1F78B4", "#FDBF6F", "#FDBF6F"),
     pch=c(17,2,16,1),
     lwd=1,
     cex=0.9,
     lty=c(1,2,1,2),
     xlim=c(0.635,0.80), ylim=c(0.73,1.00),
     las=1,cex.axis = 1.3,
     xlab = "",
     ylab = "",
     tck=-0.015,
     family = "Times New Roman")

title(xlab = "log10 (hind tibia length)", line = 2,cex.lab = 1.2,font.lab=2,family = "Times New Roman")           
title(ylab = "log10 (corneal area)", line = 3,cex.lab = 1.2,font.lab=2,family = "Times New Roman") 
title("corneal area", adj = 0.5, line = 0.6, family = "Times New Roman", cex.main=1.35)

legend("topleft", c("C", "M"),
       col = c("#1F78B4", "#FDBF6F"), pch=c(17,16),
       bty="n")




######testing shift in elevation  (grade shift - a)#####

#facet count
a3 <- sma(whole_count_log~tibia_length_log+inter, type="elevation", 
          data=eye2, multcomp = T, multcompmethod = "adjust")
multcompmatrix(a3,  sort=T)
summary(a3)
plot(a3)

#area
a4 <- sma(area_log~tibia_length_log+inter, type="elevation", 
          data=eye2, multcomp = T, multcompmethod = "adjust" )
multcompmatrix(a4,  sort=T)
summary(a4)
plot(a4)


######testing shift along common axis#####

#facet count
a5 <- sma(whole_count_log~tibia_length_log+inter, type="shift", 
          data=eye2, multcomp = T, multcompmethod = "adjust")
multcompmatrix(a5,  sort=T)
summary(a5)
plot(a5)


#area
a6 <- sma(area_log~tibia_length_log+inter, type="shift", 
          data=eye2, multcomp = T, multcompmethod = "adjust")
multcompmatrix(a6,  sort=T)
summary(a6)
plot(a6)

#####2) as a single term, separated by sex####

#divide datasets by sex
eye2.m <- eye2[which(eye2$sex=="male"),]
eye2.f <- eye2[which(eye2$sex=="female"),]

eye.m <- eye[which(eye$sex=="male"),]
eye.f <- eye[which(eye$sex=="female"),]


######testing differences in slope (beta)###### 

#facet count - male
b1 <- sma(whole_count_log~tibia_length_log*type, data=eye2.m)
summary(b1)
plot(b1)

#facet count - female
b1.1 <- sma(whole_count_log~tibia_length_log*type, data=eye2.f)
summary(b1.1)
plot(b1.1)

par(mfrow = c(1, 2))
plot(b1,main = "males")
plot(b1.1,main = "females")
par(mfrow = c(1, 1))

#area - male
b2 <- sma(area_log~tibia_length_log*type, data=eye2.m)
summary(b2)
plot(b2)

#area - female
b2.1 <- sma(area_log~tibia_length_log*type, data=eye2.f)
summary(b2.1)
plot(b2.1)


######testing shift in elevation  (grade shift - a)#####

#facet count - male
b3 <- sma(whole_count_log~tibia_length_log+type, type="elevation", data=eye2.m)
summary(b3)
plot(b3)

#facet count - female
b3.1 <- sma(whole_count_log~tibia_length_log+type, type="elevation", data=eye2.f)
summary(b3.1)
plot(b3.1)

plot(eye2.f$whole_count_log~eye2.f$tibia_length_log,
     col = eye2.f$type,
     pch = 16)


#area - male
b4 <- sma(area_log~tibia_length_log+type, type="elevation", data=eye2.m)
summary(b4)
plot(b4)

#area - female
b4.1 <- sma(area_log~tibia_length_log+type, type="elevation", data=eye2.f)
summary(b4.1)
plot(b4.1)


######testing shift along common axis#####

#facet count - male
b5 <- sma(whole_count_log~tibia_length_log+type, type="shift", data=eye2.m)
summary(b5)
plot(b5)

#facet count - female
b5.1 <- sma(whole_count_log~tibia_length_log+type, type="shift", data=eye2.f)
summary(b5.1)
plot(b5.1)


#area - male
b6 <- sma(area_log~tibia_length_log+type, type="shift", data=eye2.m)
summary(b6)
plot(b6)

#area - female
b6.1 <- sma(area_log~tibia_length_log+type, type="shift", data=eye2.f)
summary(b6.1)
plot(b6.1)


#####FDR P-VALUE CORRECTION #####
pvalues2 <- c(0.11696,
              0.0065491,
              0.50544,
              0.455,
              0.67507,
              0.7595,
              0.26785,
              1.6192e-07,
              1.482e-07,
              4.9142e-09)

adjusted.p2 <- p.adjust(pvalues2, method = "fdr", n = length(pvalues2))
adjusted.p2



#####3) now including F1 hybrids####
eye$inter <- interaction(eye$type, eye$sex)

#facet count
c1 <- sma(whole_count_log~tibia_length_log*inter, data=eye, multcomp = T, multcompmethod = "adjust")
multcompmatrix(c1, sort=T)
summary(c1)
plot(c1)

eye.m <- eye[which(eye$sex=="male"),]

c1m <- sma(whole_count_log~tibia_length_log*inter, data=eye.m, multcomp = T, multcompmethod = "adjust")
multcompmatrix(c1m, sort=T)
summary(c1m)
plot(c1m)

eye.f <- eye[which(eye$sex=="female"),]

c1f <- sma(whole_count_log~tibia_length_log*inter, data=eye.f, multcomp = T, multcompmethod = "adjust")
multcompmatrix(c1f, sort=T)
summary(c1f)
plot(c1f)


########**plot facet count major axis regressions: C, M, F1**######
dev.off()
par(mgp=c(3.5,0.6,0))
par(mar = c(5, 6, 2, 1.25))
par(lwd=1)

eye$inter<- factor(eye$inter, c("CP.male", "CP.female", 
                                "CPxMP.male", "CPxMP.female",
                                "MPxCP.male", "MPxCP.female",
                                "MP.male", "MP.female"))

#males
plot(c1, col=c("#1F78B4","#1F78B4", "#FF7F00", "#FF7F00", "grey40","grey40", "#FDBF6F","#FDBF6F"),
     #pch=c(17,2, 15,0, 18,5,  19,1),
     pch="",
     lwd=3,
     cex=0.9,
     lty=c(1,0,1,0,1,0,1,0),
     xlab = "",
     ylab = "log10 (facet count)",
     cex.axis = 1.65,cex.lab = 1.7,font.lab = 2, las=1,
     xlim=c(0.635,0.80), ylim=c(4.03,4.25),
     family = "Times New Roman", tck=-0.015)

legend("topleft", c("C","CxM", "MxC", "M"),
       col = c("#1F78B4","#FF7F00", "grey40", "#FDBF6F"), lty= 1, lwd=3,
       bty="n")

plot(c1m, col=c("#1F78B4", "#FF7F00", "grey40", "#FDBF6F"),
     #pch=c(17,15,18,19),
     pch=19,
     lwd=3,
     cex=0.25,
     lty=c(1,1,1,1),
     xlim=c(0.635,0.80), ylim=c(4.03,4.25),
     las=1,cex.axis = 1.3,
     xlab = "",
     ylab = "",
     tck=-0.015,
     family = "Times New Roman")

title(xlab = "log10 (hind tibia length)", line = 2,cex.lab = 1.2,font.lab=2,family = "Times New Roman")           
title(ylab = "log10 (facet count)", line = 3,cex.lab = 1.2,font.lab=2,family = "Times New Roman") 
title("males", adj = 0.5, line = 0.6, family = "Times New Roman", cex.main=1.35)

legend("topleft", c("C","CxM", "MxC", "M"),
       col = c("#1F78B4","#FF7F00", "grey40", "#FDBF6F"), lty= 1, lwd=3,
       bty="n")

#females
plot(c1, col=c("#1F78B4","#1F78B4", "#FF7F00", "#FF7F00", "grey40","grey40", "#FDBF6F","#FDBF6F"),
     #pch=c(17,2, 15,0, 18,5,  19,1),
     pch="",
     lwd=3,
     cex=0.9,
     lty=c(0,1,0,1,0,1,0,1),
     xlab = "",
     ylab = "log10 (facet count)",
     cex.axis = 1.65,cex.lab = 1.7,font.lab = 2, las=1,
     xlim=c(0.635,0.80), ylim=c(4.03,4.25),
     family = "Times New Roman", tck=-0.015)

legend("topleft", c("C","CxM", "MxC", "M"),
       col = c("#1F78B4","#FF7F00", "grey40", "#FDBF6F"), lty= 1, lwd=3,
       bty="n")


plot(c1f, col=c("#1F78B4", "#FF7F00", "grey40", "#FDBF6F"),
     #pch=c(17,15,18,19),
     pch=19,
     lwd=3,
     cex=0.25,
     lty=c(1,1,1,1),
     xlim=c(0.635,0.80), ylim=c(4.03,4.25),
     las=1,cex.axis = 1.3,
     xlab = "",
     ylab = "",
     tck=-0.015,
     family = "Times New Roman")

title(xlab = "log10 (hind tibia length)", line = 2,cex.lab = 1.2,font.lab=2,family = "Times New Roman")           
title(ylab = "log10 (facet count)", line = 3,cex.lab = 1.2,font.lab=2,family = "Times New Roman") 
title("females", adj = 0.5, line = 0.6, family = "Times New Roman", cex.main=1.35)

legend("topleft", c("C","CxM", "MxC", "M"),
       col = c("#1F78B4","#FF7F00", "grey40", "#FDBF6F"), lty= 1, lwd=3,
       bty="n")


#area
d1 <- sma(area_log~tibia_length_log*inter, data=eye, multcomp = T, multcompmethod = "adjust")
multcompmatrix(d1, sort=T)
summary(d1)
plot(d1)

eye.m <- eye[which(eye$sex=="male"),]

d1m <- sma(area_log~tibia_length_log*inter, data=eye.m, multcomp = T, multcompmethod = "adjust")
multcompmatrix(d1m, sort=T)
summary(d1m)
plot(d1m)

eye.f <- eye[which(eye$sex=="female"),]

d1f <- sma(area_log~tibia_length_log*inter, data=eye.f, multcomp = T, multcompmethod = "adjust")
multcompmatrix(d1f, sort=T)
summary(d1f)
plot(d1f)


########**plot corneal area major axis regressions: C, M, F1**######
dev.off()
par(mgp=c(3.5,0.6,0))
par(mar = c(5, 6, 2, 1.25))
par(lwd=1)

eye$inter<- factor(eye$inter, c("CP.male", "CP.female", 
                                "CPxMP.male", "CPxMP.female",
                                "MPxCP.male", "MPxCP.female",
                                "MP.male", "MP.female"))

#males
plot(d1m, col=c("#1F78B4", "#FF7F00", "grey40", "#FDBF6F"),
     #pch=c(17,15,18,19),
     pch=19,
     lwd=3,
     cex=0.25,
     lty=c(1,1,1,1),
     xlim=c(0.635,0.80), ylim=c(0.73,1.00),
     las=1,cex.axis = 1.3,
     xlab = "",
     ylab = "",
     tck=-0.015,
     family = "Times New Roman")

title(xlab = "log10 (hind tibia length)", line = 2,cex.lab = 1.2,font.lab=2,family = "Times New Roman")           
title(ylab = "log10 (corneal area)", line = 3,cex.lab = 1.2,font.lab=2,family = "Times New Roman") 
title("males", adj = 0.5, line = 0.6, family = "Times New Roman", cex.main=1.35)

legend("topleft", c("C","CxM", "MxC", "M"),
       col = c("#1F78B4","#FF7F00", "grey40", "#FDBF6F"), lty= 1, lwd=3,
       bty="n")

#females
plot(d1f, col=c("#1F78B4", "#FF7F00", "grey40", "#FDBF6F"),
     #pch=c(17,15,18,19),
     pch=19,
     lwd=3,
     cex=0.25,
     lty=c(1,1,1,1),
     xlim=c(0.635,0.80), ylim=c(0.73,1.00),
     las=1,cex.axis = 1.3,
     xlab = "",
     ylab = "",
     tck=-0.015,
     family = "Times New Roman")

title(xlab = "log10 (hind tibia length)", line = 2,cex.lab = 1.2,font.lab=2,family = "Times New Roman")           
title(ylab = "log10 (corneal area)", line = 3,cex.lab = 1.2,font.lab=2,family = "Times New Roman") 
title("females", adj = 0.5, line = 0.6, family = "Times New Roman", cex.main=1.35)

legend("topleft", c("C","CxM", "MxC", "M"),
       col = c("#1F78B4","#FF7F00", "grey40", "#FDBF6F"), lty= 1, lwd=3,
       bty="n")





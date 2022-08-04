getwd()
setwd("C:/Users/...")

# Loading packages
########################
library(lmerTest)
library(LMERConvenienceFunctions)
library(lme4)
library(languageR)
library(lsmeans)
library(optimx)
library(ggplot2)
library(lattice)
library(finalfit)
library(forcats)
library(plyr)
library(reshape2)

########################
# Data preparation
########################
Data <- read.delim("preTest_raw_data_plausibility.txt", encoding = "UTF-8")
Data$Participant <- as.factor(Data$Participant)
Data$Neologism <- as.factor(Data$Neologism)
Data$Condition <- as.factor(Data$Condition)
Data$Figure <- as.factor(Data$Figure)
Data$Regularity <- as.factor(Data$Regularity)
Data$Pattern <- as.factor(Data$Pattern)
Data$Score <- as.numeric(Data$Score)
Data$Unknown <- as.factor(Data$Unknown)
Data_no_NA <- droplevels(subset(Data, Score != "NA"))

#############################
# Verification of unknown words
#############################
Data_Unk <- subset(Data, Unknown == "Unknown")
ggplot(data = Data_Unk, aes(x=Neologism)) + geom_bar(stat = "count", position=position_dodge(width=0.8), width=0.7) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Unkown words")

# Now we drop the neologism that have more than 2 occurrences of "Unknown" 'unkown' in data 
# > absorptivite = 23Q05 (3), alambic = 31Q09 (6), bruine = 21Q02 (4), compressibilite = 23Q04 (5), croupier = 11Q02 (3), enclume = 31Q07 (3), 
# glanage = 41Q03 (9)
Data_Sco <- subset(Data, QN != "23Q05" & QN != "31Q09" & QN != "21Q02" & QN != "23Q04" & QN != "11Q02" & QN != "31Q07" & QN != "41Q03")
Data_Sco$QN = droplevels(Data_Sco$QN)
summary(Data_Sco)

#############################
# Data visualisation
#############################
tgc_Score_neo <- summarySE(Data_no_NA, measurevar="Score", groupvars=c("Neologism", "Regularity","Figure", "Pattern"))
tgc_Score_cdt <- summarySE(Data_no_NA, measurevar="Score", groupvars=c("condition", "Regularity", "Figure"))

ggplot(tgc_Score_neo, aes(x=Neologism, y=Score, fill = Regularity)) + facet_grid(.~Pattern, scales = "free_x") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") +
  stat_summary(fun = "mean", geom="bar", position = position_dodge(1)) + 
  geom_errorbar(aes(ymin=Score-ci, ymax=Score+ci), size=.3, width=.2, position=position_dodge(.8)) + ylim(0, 7) +
  theme_minimal() + theme(legend.position="bottom", legend.direction = "vertical", axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
                          strip.background = element_blank(), strip.text.x = element_blank()) +
  guides(fill=guide_legend(ncol=2))

ggplot(tgc_Score_cdt, aes(x=condition, y=Score, fill = condition)) +
  geom_bar(position=position_dodge(width=0.8), stat="identity", width=0.7) + expand_limits(y = c(0, 1)) +
  geom_text(aes(label=round(Score, 2)), vjust=-2, color="black", position = position_dodge(0.8), size=4.5) +
  scale_fill_manual(name = "Experimental conditions", 
                    labels = c("Highly regular metaphor", "Weakly regular metaphor", "Highly regular metonymy", "Weakly regular metonymy"),
                    values = c("#3D52A1", "#B4DDF7", "#ED875E", "#FFE3AA")) +
  geom_errorbar(aes(ymin=Score-ci, ymax=Score+ci), size=.3, width=.2, position=position_dodge(.8)) + ylim(0, 7) +
  theme_minimal() + theme(legend.position="bottom", legend.direction = "vertical", legend.text=element_text(size=12), axis.title.x=element_blank(), axis.text.x = element_blank(), 
                          axis.ticks.x=element_blank(), axis.text=element_text(size=12)) +
  guides(fill=guide_legend(ncol=2))

DataP11 <- subset(Data, Pattern == "P11")
DataP11$Pattern = droplevels(DataP11$Pattern)
DataP12 <- subset(Data, Pattern == "P12")
DataP12$Pattern = droplevels(DataP12$Pattern)
DataP13 <- subset(Data, Pattern == "P13")
DataP13$Pattern = droplevels(DataP13$Pattern)
DataP21 <- subset(Data, Pattern == "P21")
DataP21$Pattern = droplevels(DataP21$Pattern)
DataP22 <- subset(Data, Pattern == "P22")
DataP22$Pattern = droplevels(DataP22$Pattern)
DataP23 <- subset(Data, Pattern == "P23")
DataP23$Pattern = droplevels(DataP23$Pattern)
DataP31 <- subset(Data, Pattern == "P31")
DataP31$Pattern = droplevels(DataP31$Pattern)
DataP32 <- subset(Data, Pattern == "P32")
DataP32$Pattern = droplevels(DataP32$Pattern)
DataP33 <- subset(Data, Pattern == "P33")
DataP33$Pattern = droplevels(DataP33$Pattern)
DataP41 <- subset(Data, Pattern == "P41")
DataP41$Pattern = droplevels(DataP41$Pattern)
DataP42 <- subset(Data, Pattern == "P42")
DataP42$Pattern = droplevels(DataP42$Pattern)
DataP43 <- subset(Data, Pattern == "P43")
DataP43$Pattern = droplevels(DataP43$Pattern)
summary(DataP11)
summary(DataP12)

ggplot(DataP11, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) + labs(title = "Median Score for Pattern 11") +
  scale_x_discrete(labels = c("mentor", "croupier", "vigile", "barman", "pompier", "majordome", "acrobate", "arbitre", "infirmier", "troubadour"))
ggplot(DataP12, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) + labs(title = "Median Score for Pattern 12") +
  scale_x_discrete(labels = c("pistache","orange", "cacahuete", "abricot", "olive", "cerise", "haricot", "datte", "pasteque", "escalope"))
ggplot(DataP13, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) + labs(title = "Median Score for Pattern 13") +
  scale_x_discrete(labels = c("chrysalide", "sapin", "colline", "grotte", "cratere", "dune", "falaise", "eboulis", "banquise", "fougere"))
ggplot(DataP21, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) + labs(title = "Median Score for Pattern 21") +
  scale_x_discrete(labels = c("brise", "bruine", "typhon", "bourrasque", "cyclone", "seisme", "orage", "crachin", "mistral", "blizzard"))
ggplot(DataP22, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) + labs(title = "Median Score for Pattern 22") +
  scale_x_discrete(labels = c("estomac", "moustache", "trachee", "tempe", "genou", "thorax", "orteil", "cheveu", "epaule", "paupiere"))
ggplot(DataP23, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) + labs(title = "Median Score for Pattern 23") +
  scale_x_discrete(labels = c("insalubrite", "obesite", "etancheite", "compressibilite", "absorbtivite", "fluorescence", "corpulence", "sveltesse", "humidite", "adherence"))
ggplot(DataP31, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) + labs(title = "Median Score for Pattern 31") +
  scale_x_discrete(labels = c("chevalet", "babyfoot", "chaudiere", "parcmetre", "piano", "microondes", "enclume", "photocopieuse", "alambic", "congelateur"))
ggplot(DataP32, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) + labs(title = "Median Score for Pattern 32") +
  scale_x_discrete(labels = c("jambe", "dos", "doigt", "cuisse", "poignet", "poing", "abdos", "uterus", "gorge", "sourcil"))
ggplot(DataP33, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) + labs(title = "Median Score for Pattern 33") +
  scale_x_discrete(labels = c("savoir", "laicite","elegance", "generosite", "spiritualite", "pauvrete", "ivresse", "humour", "paranoia", "fougue"))
ggplot(DataP41, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) + labs(title = "Median Score for Pattern 41") +
  scale_x_discrete(labels = c("ablation", "reparation", "glanage", "racket", "bouclage", "perquisition", "liquidation", "broyage", "largage", "verification"))
ggplot(DataP42, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") + labs(title = "Median Score for Pattern 42") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) +
  scale_x_discrete(labels = c("paternalisme", "serviabilite", "sincerite", "clemence", "bienveillance", "neutralite", "sexisme", "arrogance", "docilite", "immaturite"))
ggplot(DataP43, aes(x=QN, y=Score)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + xlab(label = "Neologism") + labs(title = "Median Score for Pattern 43") +
  stat_summary(fun.y = "median", geom="bar", position = position_dodge(1)) + 
  stat_summary(fun.ymin = "min", fun.ymax = "max", geom = "errorbar", color = "grey40",position = position_dodge(1), width=.2) +
  scale_x_discrete(labels = c("safran", "vinaigre", "levure", "moutarde", "vanille", "lait", "ketchup", "sel", "beurre", "huile"))


#############################
# Correlation between mean plausibility and mean regularity by pattern
#############################
pattern <- c("P11", "P12", "P13", "P14", "P15", "P16", "P21", "P22", "P23", "P24", "P25", "P26")
plausibility <- c(4.028708, 3.398058, 4.173913, 5.487685, 4.382775, 3.935644, 3.159204, 4.419048, 3.588517, 5.423077, 4.214286, 4.649746)
regularity <- c(2.481481, 2.592593, 2.666667, 4.407407, 4.592593, 5.518519, 3.000000, 3.148148, 3.555556, 4.333333, 4.629630, 5.481481)
cor(plausibility, regularity, method = "spearman")

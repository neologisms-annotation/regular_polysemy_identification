getwd()
setwd("C:/Users/...")

########################
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
  # absorptivite = 23Q05 (3), alambic = 31Q09 (6), bruine = 21Q02 (4), compressibilite = 23Q04 (5), 
  # croupier = 11Q02 (3), enclume = 31Q07 (3), glanage = 41Q03 (9)
Data_Sco <- subset(Data, QN != "23Q05" & QN != "31Q09" & QN != "21Q02" & QN != "23Q04" & QN != "11Q02" & QN != "31Q07" & QN != "41Q03")
Data_Sco$QN = droplevels(Data_Sco$QN)
summary(Data_Sco)

#############################
# Plausibility scores description
#############################
# Data of the 7 words per pattern that have the highest score
Data <- droplevels(subset(Data, Neologism != "croupier" & Neologism != "arbitre" & Neologism != "infirmier" & 
                            Neologism != "abricot" & Neologism != "haricot" & Neologism != "escalope" &
                            Neologism != "sapin" & Neologism != "grotte" & Neologism != "cratere" & 
                            Neologism != "brise" & Neologism != "bruine" & Neologism != "crachin" &
                            Neologism != "estomac" & Neologism != "thorax" & Neologism != "orteil" & 
                            Neologism != "compressibilite" & Neologism != "absorptivite" & Neologism != "adherence" &
                            Neologism != "piano" & Neologism != "enclume" & Neologism != "alambic" & 
                            Neologism != "dos" & Neologism != "cuisse" & Neologism != "abdos" &
                            Neologism != "generosite" & Neologism != "spiritualite" & Neologism != "humour" & 
                            Neologism != "reparation" & Neologism != "glanage" & Neologism != "bouclage" &
                            Neologism != "sincerite" & Neologism != "neutralite" & Neologism != "docilite" & 
                            Neologism != "safran" & Neologism != "moutarde" & Neologism != "huile"))
NeoMean <- summarySE(na.omit(Data), measurevar="Score", groupvars = c("Neologism", "Pattern"))
summarySE(NeoMean, measurevar="Score", groupvars = c("Pattern"))
summarySE(na.omit(Data), measurevar="Score", groupvars = c("Pattern"))

#############################
# Correlation between plausibility and regularity
#############################
DataCOR <- read_excel("preTest_results.xlsx")
# Coding variables as nominal variables
DataCOR$Pattern <- as.factor(DataCOR$Pattern)
DataCOR$Senses <- as.factor(DataCOR$Senses)
summary(DataCOR)
cor.test(DataCOR$Reg_num, DataCOR$Plausibility_w, method = "pearson") # mean plausibility per target word

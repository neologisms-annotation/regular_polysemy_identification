################################################################################################
# Script for the data analysis
################################################################################################
getwd()
setwd("C:/Users/...")


#############################
# Loading packages
#############################
library(lmerTest)
library(LMERConvenienceFunctions)
library(lme4)
library(languageR)
library(lsmeans)
library(optimx)
library(lattice)
library(finalfit)
library(effects)
library(plyr)
library(forcats)
library(ordinal)
library(MASS)
library(readxl)
library(ggplot2)
library(gridExtra)
library(car)
library(jtools)
library(kableExtra)
library(partR2)
library(MuMIn)

# Preparing data
#############################
# Preparing data
#############################
Data_part <- read.delim("exp_data_part.txt", encoding = "ANSI")
Data_part$Gender <- as.factor(Data_part$Gender)
Data_part$Version <- as.factor(Data_part$Version)
Data_part <- droplevels(subset(Data_part, Participant != "P006" & Participant != "P013" &  Participant != "P072" &
                             Participant != "P087" & Participant != "P124" & Participant != "P133"))
summary(Data_part)

Datac0 <- read.delim("exp_data.txt", encoding = "ANSI")
Datac0 <- droplevels(subset(Datac0, Novelty != "Distractor" & Novelty != "Practice"))
Datac0$Participant <- as.factor(Datac0$Participant)
Datac0$ID <- as.factor(Datac0$ID)
Datac0$Word <- as.factor(Datac0$Word)
Datac0$Identification <- as.factor(Datac0$Identification)
Datac0$Pattern <- as.factor(Datac0$Pattern)
Datac0$C_S1 <- as.factor(Datac0$C_S1)
Datac0$C_S2 <- as.factor(Datac0$C_S2)
Datac0$C_change <- as.factor(Datac0$C_change)
Datac0$Plausibility <- as.numeric(Datac0$Score)
Datac0$Response <- as.factor(Datac0$Response)
Datac0$Version <- as.factor(Datac0$Version)
Datac0$Response_num <- as.factor(Datac0$Response)
levels(Datac0$Response_num)[levels(Datac0$Response_num)=="N"] <- "0"
levels(Datac0$Response_num)[levels(Datac0$Response_num)=="Y"] <- "1"
Datac0$Response_num <- as.numeric(Datac0$Response_num)
Datac0$Identification <- as.factor(Datac0$Identification)
Datac0$Regularity <- as.numeric(Datac0$Reg_num)
Datac0$Figure <- as.factor(Datac0$Figure)
Datac0$Novelty <- as.factor(Datac0$Novelty)

# All data from a given participant will be eliminated if positive answer rate on neutral stimuli is greater than 2.5 SD above the average rate
Datac1 <- droplevels(subset(Datac0, Novelty == "Source"))
Datac1$mn_iden_part <- ave(Datac1$Response_num, Datac1$Participant, FUN=mean)
mn_original <- mean(Datac1$Response_num[Datac1$Novelty == "Source"])
sd_original <- sd(Datac1$Response_num[Datac1$Novelty == "Source"])
max_iden <- mn_original+(2.5*sd_original)
min_iden <- mn_original-(2.5*sd_original)

Datac2 <- droplevels(subset(Datac1, mn_iden_part < max_iden & mn_iden_part > min_iden))
levels(Datac1$Participant)
levels(Datac2$Participant) # Participant P006, P013,  P072, P087, P124 need to be removed
# We also remove P133 and since their mean is also really high

# Answers shorter than 1000ms will be eliminated from the analysis data (3862 obs. vs. 3858 obs.) > 4 trials, i.e. 0.1014% of Datac0
Data00 <- droplevels(subset(Datac0, RT > 1))
Data0 <- droplevels(subset(Datac0, RT > 1, Participant != "P006" & Participant != "P013" &  Participant != "P072" &
                             Participant != "P087" & Participant != "P124" & Participant != "P133"))
Data <- droplevels(subset(Datac0, RT > 1 & Participant != "P006" & Participant != "P013" &  Participant != "P072" &
                            Participant != "P087" & Participant != "P124" & Participant != "P133" & Novelty == "Novel"))
summary(Data)

##################################
# Main analyses
##################################
#############################
# Effect of novelty
#############################
# Choice of the best model for the random intercepts
glmer_rnd0 = glmer(Identification ~ 1 + (1 | Word),data=Data0, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_rnd1 = glmer(Identification ~ 1 + (1 | Participant),data=Data0, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_rnd2 = glmer(Identification ~ 1 + (1 | Pattern),data=Data0, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_rnd3 = glmer(Identification ~ 1 + (1 | Participant) + (1 | Word), data=Data0, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_rnd4 = glmer(Identification ~ 1 + (1 | Participant) + (1 | Pattern), data=Data0, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_rnd5 = glmer(Identification ~ 1 + (1 | Pattern) + (1 | Word), data=Data0, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_rnd6 = glmer(Identification ~ 1 + (1 | Participant) + (1 | Pattern) + (1 | Word), data=Data0, family="binomial", control=glmerControl(optimizer="bobyqa"))
# Choose model with lowest AIC
AIC(glmer_rnd0, glmer_rnd1, glmer_rnd2, glmer_rnd3, glmer_rnd4, glmer_rnd5, glmer_rnd6)

glmer_nov = glmer(Identification ~ Novelty + (1 | Participant) + (1 | Pattern) + (1 | Word), data=Data0, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmer_nov, type = "III")
plot(predictorEffects(glmer_nov, ~ Novelty), lines=list(multiline=TRUE), main = NULL,
     axes=list(grid=TRUE, y=list(lab="Rate of identification", type ="response", lim=c(0, 1))))
glmer_nov = glmer(Identification ~ Novelty + (1 | Participant) + (1 | Pattern) + (1 | Word), data=Data0, family="binomial", control=glmerControl(optimizer="bobyqa"))
contrasts(Data0$Novelty) <-contr.Sum(levels(Data0$Novelty))
summary(glmer_nov)

# We test the same effect without removing any participant
glmer_nov2 = glmer(Identification ~ Novelty + (1 | Participant) + (1 | Pattern) + (1 | Word), data=Data00, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmer_nov2, type = "III")
plot(predictorEffects(glmer_nov2, ~ Novelty))
contrasts(Data00$Novelty) <-contr.Sum(levels(Data00$Novelty))
summary(glmer_nov2)

#############################
# Effect of figure, regularity and frequency
#############################
# Choice of the best model for the random intercepts
glmer_ran0 = glmer(Identification ~ 1 + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_ran1 = glmer(Identification ~ 1 + (1 | Participant),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_ran2 = glmer(Identification ~ 1 + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
# Choose model with lowest AIC
AIC(glmer_ran0, glmer_ran1, glmer_ran2)

# Choice of the best model for the fixed effects
glmer_id1 = glmer(Identification ~ Reg_num  + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id2 = glmer(Identification ~ Figure + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id3 = glmer(Identification ~ Freq + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id4 = glmer(Identification ~ Reg_num + Figure + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id5 = glmer(Identification ~ Freq + Reg_num  + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id6 = glmer(Identification ~ Freq + Figure + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id7 = glmer(Identification ~ Freq + Reg_num + Figure + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
AIC(glmer_id1, glmer_id2,glmer_id3,glmer_id4,glmer_id5,glmer_id6,glmer_id7)

# check interaction
glmer_id4i = glmer(Identification ~ Reg_num * Figure  + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
AIC(glmer_id4i,glmer_id4)
Anova(glmer_id4i, type = "III")

# check random slopes
glmer_idR1 = glmer(Identification ~ Reg_num + Figure  + (1 + Reg_num | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_idR2 = glmer(Identification ~ Reg_num + Figure  + (1 + Figure | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_idR3 = glmer(Identification ~ Reg_num + Figure  + (1 + Reg_num + Figure | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
AIC(glmer_id2,glmer_idR1, glmer_idR2, glmer_idR3)
#glmer_idenR1 do not converge

Anova(glmer_idR2)
glmer_id_R2 = glmer(Identification ~ Regularity + Figure  + (1 + Figure | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
plot(predictorEffects(glmer_id_R2, ~ Regularity + Figure), lines=list(multiline=TRUE),
     axes=list(grid=TRUE, y=list(lab="Prob(identification)", type ="response", lim=c(0, 1))))
contrasts(Data$Figure) <-contr.Sum(levels(Data$Figure))
glmer_id_R2 = glmer(Identification ~ Regularity + Figure  + (1 + Figure | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
summary(glmer_id_R2)
r.squaredGLMM(glmer_id_R2)

##################################
# Post-hoc analyses
##################################
#############################
# Plausibility
#############################
# Mean plausibility per pattern
glmerPL = glmer(Identification ~ Figure + Plausibility + (1 + Figure | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmerPL, type = "III")
r.squaredGLMM(glmerPL)
plot(predictorEffects(glmerPL, ~ Figure + Plausibility), rows = 1, cols = 2)
contrasts(Data$Figure) <-contr.Sum(levels(Data$Figure))
glmerPL = glmer(Identification ~ Figure + Plausibility + (1 + Figure | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
summary(glmerPL)

#############################
# Concreteness
#############################
glmerC1 = glmer(Identification ~ Reg_num + Figure + C_S1 + (1 + Figure | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerC2 = glmer(Identification ~ Reg_num + Figure + C_S2 + (1 + Figure | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerC3 = glmer(Identification ~ Reg_num + Figure + C_S1 * C_S2 + (1 + Figure | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmerC3, type = "III")
r.squaredGLMM(glmerC3)
summary(glmerC3)

#############################
# Test of the variables as only predictor
#############################
# Simple single predictors
glmerREG = glmer(Identification ~ Reg_num + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerFIG = glmer(Identification ~ Figure  + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerCS1 = glmer(Identification ~ C_S1 + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerCS2 = glmer(Identification ~ C_S2 + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerPLA = glmer(Identification ~ Plausibility + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerFRE = glmer(Identification ~ Freq + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmerREG, type = "III")
r.squaredGLMM(glmerREG)




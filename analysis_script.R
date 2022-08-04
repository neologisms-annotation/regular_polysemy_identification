################################################################################################
# Script for the data analysis
################################################################################################
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
library(partR2) # for adjusted R square

# Preparing data
########################
Data_part <- read.delim("exp_data_part.txt", encoding = "ANSI")
Data_part$Gender <- as.factor(Data_part$Gender)
Data_part$Version <- as.factor(Data_part$Version)
Data0 <- droplevels(subset(Datac0, RT > 1, Participant != "P006" & Participant != "P013" &  Participant != "P072" &
                             Participant != "P087" & Participant != "P124" & Participant != "P133"))
summary(Data_part)

Datac0 <- read.delim("exp_data.txt", encoding = "ANSI")
# Coding variables as nominal variables
Datac0 <- droplevels(subset(Datac0, Novelty != "Distractor" & Novelty != "Practice"))
Datac0$Participant <- as.factor(Datac0$Participant)
Datac0$ID <- as.factor(Datac0$ID)
Datac0$Word <- as.factor(Datac0$Word)
Datac0$Identification <- as.factor(Datac0$Identification)
Datac0$Pattern <- as.factor(Datac0$Pattern)
Datac0$C_S1 <- as.factor(Datac0$C_S1)
Datac0$C_S2 <- as.factor(Datac0$C_S2)
Datac0$C_change <- as.factor(Datac0$C_change)
Datac0$Response <- as.factor(Datac0$Response)
Datac0$Reg_num_fac <- as.factor(Datac0$Reg_num)
Datac0$Version <- as.factor(Datac0$Version)
Datac0$Response_num <- as.factor(Datac0$Response)
levels(Datac0$Response_num)[levels(Datac0$Response_num)=="N"] <- "0"
levels(Datac0$Response_num)[levels(Datac0$Response_num)=="Y"] <- "1"
Datac0$Response_num <- as.numeric(Datac0$Response_num)
Datac0$Identification <- as.factor(Datac0$Identification)
Datac0$Reg_degree <- as.factor(Datac0$Reg_degree)
Datac0$Regularity <- as.numeric(Datac0$Reg_num)
Datac0$Figure <- as.factor(Datac0$Figure)
Datac0$Novelty <- as.factor(Datac0$Novelty)
# We create a variable "condition" that is only used to generate figures
Datac0$Condition <- with(Datac0, interaction(Reg_degree, Figure), drop = TRUE )
levels(Datac0$Condition)[levels(Datac0$Condition)=="High.Metaphor"] <- "1"
levels(Datac0$Condition)[levels(Datac0$Condition)=="Low.Metaphor"] <- "2"
levels(Datac0$Condition)[levels(Datac0$Condition)=="High.Metonymy"] <- "3"
levels(Datac0$Condition)[levels(Datac0$Condition)=="Low.Metonymy"] <- "4"
summary(Datac0)

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
Data0 <- droplevels(subset(Datac0, RT > 1, Participant != "P006" & Participant != "P013" &  Participant != "P072" &
                             Participant != "P087" & Participant != "P124" & Participant != "P133"))
Data <- droplevels(subset(Datac0, RT > 1 & Participant != "P006" & Participant != "P013" &  Participant != "P072" &
                            Participant != "P087" & Participant != "P124" & Participant != "P133" & Novelty == "Novel"))
summary(Data)


################################################################################
# Data visulisation
################################################################################
#############################
# Participants results on source item identification
#############################
Data_source <- droplevels(subset(Datac1, Novelty == "Source"))
tgc_source <- ddply(Data_source, c("Participant"), summarise,
                      N    = length(Identification),
                      prob = length(Identification[Identification == "Y"])/N,
                      wald_ci = 1.96*(sqrt(prob*(1-prob)/N))) # Wald 95% confidence interval
ggplot(tgc_source, aes(x=Participant, y=prob)) + expand_limits(y = c(0, 1)) +
  geom_bar(position=position_dodge(width=0.8), stat="identity", width=0.7) +      # Thinner lines
  geom_errorbar(aes(ymin=prob-wald_ci, ymax=prob+wald_ci), size=.3,    # Thinner lines + with 95% confidence interval
                width=.2, position=position_dodge(.8)) +
  ylab("Rate of source word identification") + xlab("Participant") +
  theme_minimal() + theme(axis.text.x = element_text(angle=90, vjust= 0.3), legend.position="bottom", legend.direction = "vertical",
                          axis.title.x = element_blank(), strip.text.x = element_blank()) +
  guides(fill=guide_legend(ncol=2))

#############################
# Identification in general
#############################
# Run the functions length, mean, and sd on the value of "iden" for each group, 
# broken down by reg_degree + figure + identifier
tgc_iden_cdt <- ddply(Data, c("Condition", "Novelty"), summarise,
                      N    = length(Identification),
                      prob = length(Identification[Identification == "Y"])/N,
                      wald_ci = 1.96*(sqrt(prob*(1-prob)/N))) # Wald 95% confidence interval

ggplot(tgc_iden_cdt, aes(x=Novelty, y=prob, fill = Condition)) + 
  facet_grid(. ~ Condition, scales = "free_x") + expand_limits(y = c(0, 1)) +
  geom_bar(position=position_dodge(width=0.8), stat="identity", width=0.7) +      # Thinner lines
  geom_errorbar(aes(ymin=prob-wald_ci, ymax=prob+wald_ci), size=.3,    # Thinner lines + with 95% confidence interval
                width=.2, position=position_dodge(.8)) +
  ylab("Empirical probability of identification") + xlab("Target word") +
  scale_fill_manual(name = "Experimental conditions", 
                    labels = c("Highly regular metaphor", "Weakly regular metaphor", "Highly regular metonymy", "Weakly regular metonymy"),
                    values = c("#3D52A1", "#B4DDF7", "#ED875E", "#FFE3AA")) +
  theme_minimal() + theme(axis.text.x = element_text(angle=0, vjust= 0.3), legend.position="bottom", legend.direction = "vertical",
                          axis.title.x = element_blank(), strip.text.x = element_blank()) +
  guides(fill=guide_legend(ncol=2))

################################################################################################
# Identification analysis
################################################################################################
##################################
# Effect of novelty
##################################
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
# df      AIC
# glmer_rnd0  2 15165.75
# glmer_rnd1  2 15476.86
# glmer_rnd2  2 15271.31
# glmer_rnd3  3 15001.31
# glmer_rnd4  3 15114.34
# glmer_rnd5  3 15133.56
# glmer_rnd6  4 14969.12    best model

glmer_nov = glmer(Identification ~ Novelty + (1 | Participant) + (1 | Pattern) + (1 | Word), data=Data0, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmer_nov, type = "III")
plot(predictorEffects(glmer_nov, ~ Novelty))

contrasts(Data0$Novelty) <-contr.Sum(levels(Data0$Novelty))
summary(glmer_nov)

# R2s and partial R2s
#(R2_nov <- partR2(glmer_nov,  partvars = c("Novelty"), R2_type = "marginal", nboot = 100, CI = 0.95))
#forestplot(R2_nov, type = "R2", line_size = 0.7, text_size = 14, point_size = 3)

##################################
# Effect of figure, regularity and frequency
##################################
# Choice of the best model for the random intercepts
glmer_ran0 = glmer(Identification ~ 1 + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_ran1 = glmer(Identification ~ 1 + (1 | Participant),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_ran2 = glmer(Identification ~ 1 + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
# Choose model with lowest AIC
AIC(glmer_ran0, glmer_ran1, glmer_ran2)
# df      AIC
# glmer_ran0  2 7536.784
# glmer_ran1  2 8031.095
# glmer_ran2  3 7002.374    best model


# Choice of the best model for the fixed effects
# main effects only:
glmer_id1 = glmer(Identification ~ Reg_num  + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id2 = glmer(Identification ~ Figure + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id3 = glmer(Identification ~ Freq + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id4 = glmer(Identification ~ Reg_num + Figure + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id5 = glmer(Identification ~ Freq + Reg_num  + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id6 = glmer(Identification ~ Freq + Figure + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id7 = glmer(Identification ~ Freq + Reg_num + Figure + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
AIC(glmer_id1, glmer_id2,glmer_id3,glmer_id4,glmer_id5,glmer_id6,glmer_id7)
# df      AIC
# glmer_id1  4 6998.462
# glmer_id2  4 6962.870
# glmer_id3  4 7004.327
# glmer_id4  5 6960.787
# glmer_id5  5 7000.461
# glmer_id6  5 6964.654
# glmer_id7  6 6962.417
Anova(glmer_id2)
Anova(glmer_id4)
Anova(glmer_id7)

glmer_id1 = glmer(Identification ~ Reg_num  + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id01 = glmer(Identification ~ Reg_degree  + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
AIC(glmer_id1, glmer_id01)
# df      AIC
# glmer_id1   4 6998.462
# glmer_id01  4 6996.338    degrees of regularity is slightly better

# check interaction
glmer_id4i = glmer(Identification ~ Reg_num * Figure  + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
AIC(glmer_id4i,glmer_id4)
# df      AIC
# glmer_id4i  6 6960.356
# glmer_id4   5 6960.787
Anova(glmer_id4i, type = "III")

# check random slopes
glmer_idR1 = glmer(Identification ~ Reg_num + Figure  + (1 + Reg_num | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_idR2 = glmer(Identification ~ Reg_num + Figure  + (1 + Figure | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_idR3 = glmer(Identification ~ Reg_num + Figure  + (1 + Reg_num + Figure | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
AIC(glmer_id2,glmer_idR1, glmer_idR2, glmer_idR3)
#glmer_idenR1 do not converge
# df      AIC
# glmer_id2   4 6962.870
# glmer_idR1  7 6963.355
# glmer_idR2  7 6949.197    best model
# glmer_idR3 10 6951.557

Anova(glmer_idR2)
# Chisq Df Pr(>Chisq)    
# Reg_num  4.1974  1    0.04049 *  
# Figure  48.2206  1  3.809e-12 ***

plot(predictorEffects(glmer_idR2, ~ Regularity + Figure), lines=list(multiline=TRUE),
     axes=list(grid=TRUE, y=list(lab="Prob(identification)", type ="response", lim=c(0, 1))))

contrasts(Data$Figure) <-contr.Sum(levels(Data$Figure))
glmer_idR2 = glmer(Identification ~ Regularity + Figure  + (1 + Figure | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
summary(glmer_idR2)

# Maria say that the R-squared works well for ANOVA but it doesn't quite well for more complex models :o
library(MuMIn)  # most accepted and used one
r.squaredGLMM(glmer_idR2)
# R2m       R2c
# theoretical 0.10325511 0.4075210
# delta       0.09070357 0.3579833
# The marginal R squared values are those associated with your fixed effects,
 # the conditional ones are those of your fixed effects plus the random effects.
 # Usually we will be interested in the marginal effects.

glmer_id_0 = glmer(Identification ~ Reg_num  + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_id_1 = glmer(Identification ~ Figure  + (1 | Participant) + (1 | Word),data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmer_id_1, type = "III")
r.squaredGLMM(glmer_id_1)


##### We check the same model without P32 and P12
Data_10P <- droplevels(subset(Data, Pattern != "P12" & Pattern != "P32"))

glmer_id10p1 = glmer(Identification ~ Freq + Reg_num + Figure + (1 + Figure | Participant) + (1 | Word),data=Data_10P, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmer_id10p1, type = "III")
glmer_id10p2 = glmer(Identification ~ Figure + (1 + Figure | Participant) + (1 | Word),data=Data_10P, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmer_id10p2, type = "III")
#   Chisq               Df Pr(>Chisq)    
#  (Intercept) 39.074  1  4.081e-10 ***
#  Figure      59.254  1  1.386e-14 ***

##################################
# Is there other predictors ?
##################################
# Here, we check the effect of three supplementary variables
  # C_S1 : the concreteness of the sense 1 (concrete vs. abstract)
  # C_S2 : the concreteness of the sense 2 (concrete vs. abstract)
  # C_Change : indicate if there is a change frome concrete to abstract or from abstract to concrete between S1 and S2 (y/n)
  # Score : is the by-pattern mean plausibility score given to stimuli by participant during the preliminary experiment (numerical)

# We have the following hypothesese
  # H1 : C_S1 : concrete S1 will be more detected than abstract ones
  # H2 : C_S2 : concrete S2 will be more detected than abstract ones
  # H3 : C_change : there may be an interaction with one or the two other variables
        # probably with C_S2 : 
          # abstract with no change will be less detected than concrete with no change
          # abstract and concrete with change would have a detection rate in between abstract and concrete with no change

# Simple single predictors
glmerC01 = glmer(Identification ~ C_S1 + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerC02 = glmer(Identification ~ C_S2 + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerC03 = glmer(Identification ~ C_change + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerC04 = glmer(Identification ~ Score + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerC05 = glmer(Identification ~ Freq + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmerC05, type = "III")

# Interaction
glmerC1 = glmer(Identification ~ C_S1*C_change + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerC2 = glmer(Identification ~ C_S2*C_change + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmerC3 = glmer(Identification ~ C_S1*C_change + C_S2*C_change + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmerC3, type = "III")
r.squaredGLMM(glmerC04)

# With our principal predictors
glmer_C31 = glmer(Identification ~ Reg_num + Figure + C_S1 * C_change + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_C32 = glmer(Identification ~ Reg_num + Figure + C_S2 * C_change + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_C33 = glmer(Identification ~ Reg_num + Figure + Score + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
glmer_C34 = glmer(Identification ~ Figure + Score + (1 | Participant) + (1 | Word), data=Data, family="binomial", control=glmerControl(optimizer="bobyqa"))
Anova(glmer_C34, type = "III")
r.squaredGLMM(glmer_C34)
summary(glmer_C34)

plot(predictorEffects(glmer_C31, ~ Reg_num + Figure + C_S1 * C_change))
plot(predictorEffects(glmer_C32, ~ Reg_num + Figure + C_S2 * C_change))
plot(predictorEffects(glmer_C33, ~ Reg_num + Figure + Score), rows = 1, cols = 3)
plot(predictorEffects(glmer_C34, ~ Figure + Score))


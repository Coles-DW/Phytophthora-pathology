#data
EL <- read.csv(file = 'Stats_EL_3.csv')

#stats
library(visreg)
library(car)
library(emmeans)
library(multcompView)
library(multcomp)
install.packages("ggeffects")
library(ggeffects)
#test for normality of data (If p-value >0.05 means data not sig. diff. from normal distr. Therefore assume normality)

shapiro.test(EL$TRT)
shapiro.test(EL$dsiRNA2)
shapiro.test(EL$dsiRNAsc)

str(EL)
EL$HPI<-as.factor(EL$HPI)
EL$TRT<-as.factor(EL$TRT)


#create model:
dsi1 <- lm(REL ~ TRT*HPI, data=EL)

#look @ marginal effects

mydf <- ggpredict(dsi1, terms = c("HPI", "TRT"))
mydf

plot(mydf)




dsi2 <- lm(dsiRNA2 ~ HPI, data=EL)
dsisc <- lm(dsiRNAsc ~ HPI, data=EL)

summary(dsi1)
summary(dsi2)
summary(dsisc)

#Anova:
anova(dsi1)
anova(dsi2)
anova(dsisc)

#residualplot
residualPlot(sa)
residualPlot(ja)
residualPlot(aba)
residualPlot(saja)

#qqplot
qqPlot(sa)
qqPlot(ja)
qqPlot(aba)
qqPlot(saja)

#Statistical difference in means between treatments
sa.emm <- emmeans(sa, ~Treatment)
cld(sa.emm)
ja.emm <- emmeans(ja, ~Treatment)
cld(ja.emm)
aba.emm <- emmeans(aba, ~Treatment)
cld(aba.emm)
saja.emm <- emmeans(saja, ~Treatment)
cld(saja.emm)


#read data
EL2 <- read.csv("12_EL.csv")
head(EL2)

library(ggplot2)
library(plyr)
library(tidyverse)


str(EL2)
EL2$Treatment<-as.factor(EL2$Treatment)


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

EL3 <- data_summary(EL2, varname="REL",
                    groupnames= "Treatment")
str(EL3)

head(EL3)

ggplot(EL3, aes(x=Treatment, y=REL)) +
  geom_bar(stat="identity", position=position_dodge())+
  labs(x="\n Treatment", y = "% Electrolyte leakage \n")+
  scale_y_continuous(limits=c(0,100), expand = c(0,0))+
  theme_classic()+
  scale_fill_grey()+
  theme(text=element_text(family="Times New Roman" , size=18))+
  geom_errorbar(aes(ymin=REL-sd, ymax=REL+sd), width=.2,
                position=position_dodge(.9))

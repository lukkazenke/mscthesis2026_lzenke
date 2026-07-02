#########MSc Thesis Statistics- PART 1 FACS - Lukka Zenke - last modified: 09 Sept 2025########
###############################################################################################

#FILES#
Plantgrowth <- read.delim("~/Desktop/BioOcean/Doing Science/[Bioc-ws-23-24] Doing Science block course/Plantgrowth.txt")

# PREPARATIONS #
str(Plantgrowth) #always check, here nutrients needs to be transformed to categorical
Plantgrowth$Nutrients<-as.factor(Plantgrowth$Nutrients)
str(Plantgrowth)
table(Plantgrowth$Nutrients) #is it a balanced design? how many observations per level

#PLOTTING#
dotchart(Plantgrowth$Biomass, groups=factor(Plantgrowth$Nutrients), xlab="Biomass")
#checking normality:
hist(Plantgrowth$Biomass, breaks=14)
qqnorm(Plantgrowth$Biomass)
qqline(Plantgrowth$Biomass)
shapiro.test(Plantgrowth$Biomass) #suprisingly is significant although data is normal, don't use for bigger numbers!
#Boxplot
boxplot(Plantgrowth$Biomass~Plantgrowth$Nutrients,
        xlab="Nutrients", ylab="Biomass (g/m2)")
#Box-whisker Plot with ggplot2
library(ggplot2)
library(Hmisc)
Figure1<-ggplot(Plantgrowth, aes(Nutrients, Biomass)) #Layer 1: defining the axes
Figure1+
  geom_boxplot() #Layer2: plotting boxplots
Figure1+
  geom_boxplot()+
  labs(x="Nutrients addition in %", y="Biomass (g/m2") #Layer3: labels for axes
Figure1+
  geom_boxplot()+
  labs(x="Nutrients addition in %", y="Biomass (g/m2")+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(colour="black"),
        axis.title=element_text(size=16, colour="black"),
        axis.text=element_text(size=12))
#Mean-error plot with ggplot2
Figure2<-ggplot(Plantgrowth, aes(Nutrients, Biomass))
Figure2+
  stat_summary(fun=mean, geom="bar",
               colour="white", linewidth=1,
               fill="grey60")+
  stat_summary(fun.data=mean_cl_normal,
               geom="errorbar", position=position_dodge(width=0.9),
               width=0.0, colour="black", linetype=1,
               size=0.8)+ #adds Cls as error bars
  theme(panel.background = element_blank(),
        axis.line=element_line(colour="black"),
        axis.title=element_text(size=16, colour="black"),
        axis.text=element_text(size=12)) #check visually if means are different

#MODELLING#
#information about treatment levels
tapply(Plantgrowth$Biomass, Plantgrowth$Nutrients, mean) #look at mean of each group
tapply(Plantgrowth$Biomass, Plantgrowth$Nutrients, sd) #check variances
Model1<-aov(Plantgrowth$Biomass~Plantgrowth$Nutrients)

#DIAGNOSTICS INTERLUDE#
#homogeneity of variances
plot(resid(Model1)~fitted(Model1)) #in bands because of categorical independent variable
abline(h=0, lwd=2, lty=2, col="black") #here fitted values are group mean for each group
fligner.test(Plantgrowth$Biomass~Plantgrowth$Nutrients) #should be insignificant
#normality of errors
hist(resid(Model1))
shapiro.test(resid(Model1))
qqnorm(resid(Model1))
qqline(resid(Model1))
#influential data points
plot(cooks.distance(Model1), type="h")
#diagnostics say: All is Well

#CALL FOR OUTPUT#
summary(Model1) #shows ANOVA table = highly significant effect
summary.lm(Model1) #regression output of ANOVA -> can get R2 here
#Posthoc test
TukeyHSD(Model1) #Tukey-Test, compares all groups with each other

#RETROSPECTIVE POWER ANALYSIS#
power.anova.test(groups=6,n=NULL,
                 between.var=2371534, #read from ANOVA table Mean Sq
                 within.var=50537, #Mean Sq Residuals
                 sig.level=0.05,
                 power=1) #search for n where power is 1 or 0.95

#OTHER ANOVAS#
#if data is not normal, use Kruskal-Wallis ANOVA:
kruskal.test(Plantgrowth$Biomass~Plantgrowth$Nutrients) #has less power because it works with ranks
library(pgirmess)
kruskalmc(Plantgrowth$Biomass,Plantgrowth$Nutrients) #post-hoc test for Kruskal-Wallis
#if variances are not homogenous use welch-adjustes ANOVA:
oneway.test(Plantgrowth$Biomass~Plantgrowth$Nutrients) #affects df because of correction

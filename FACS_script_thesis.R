####### FACS analysis - Lukka Zenke thesis - last updated 16 March 2026#####
rm(list=ls(all=TRUE))
###LOAD DATA ###
FACS_uncut<-read.csv("MScthesis_FACS_data.csv", header = TRUE)
FACS <- FACS_uncut[-c(37,38,39,40, 41:45),] #cut DSM507
#put the variable in the desired order bc R otherwise sorts alphabetically
FACS$vibrio_strain <- factor(FACS$vibrio_strain, c("Hal054","Hal281", "M60_M31", "M60_M70", "SA48"))
FACS$ankyrin_present <- as.factor(FACS$ankyrin_present)
FACS$self_vs_foreign <- as.factor(FACS$self_vs_foreign)
FACS$time_min <- as.factor(FACS$time_min)
FACS$sample_id <- as.factor(FACS$sample_id)
str(FACS)


###CHECK NORMALITY###
#checking for outliers in dependent variable
dotchart(FACS$phago_avg_DAPI_adjusted)
#checking for balanced design
table(FACS$vibrio_strain, FACS$time_min) #here unbalanced due to losses in DMZ507
#checking for normality
hist(FACS$phago_avg_DAPI_adjusted, breaks=20, col="pink") #a little weird, but still normal
qqnorm(FACS$phago_avg_DAPI_adjusted)
qqline(FACS$phago_avg_DAPI_adjusted) #fine by me
shapiro.test(FACS$phago_avg_DAPI_adjusted) #useless as usual

###CHECK HOMOGENIETY OF VARIANCES
library(car)
leveneTest(FACS$phago_avg_DAPI_adjusted~FACS$vibrio_strain, data = FACS) #variances are OK!

#MEANS AND SDs
aggregate(FACS$phago_avg_DAPI_adjusted ~ FACS$vibrio_strain + FACS$time_min, FUN = mean)
aggregate(FACS$phago_avg_DAPI_adjusted ~ FACS$vibrio_strain + FACS$time_min, FUN = sd)

###PLOTTING###
library(ggplot2)
boxplot(FACS$phago_avg_DAPI_adjusted~FACS$vibrio_strain) #probably nothing significant

#self vs non-self plot
#self_plot <- ggplot(FACS, aes(x = vibrio_strain, y = phago_avg_DAPI_adjusted, colour = self_vs_foreign)) +
  #geom_jitter(width = 0.2, alpha = 0.6, color = "black")+
 # geom_boxplot(outlier.shape = NA) +  # hide default outliers
  #scale_colour_manual(name= "Origin", values=c("#CC6666", "#9999CC"),
                      labels= c("foreign", "H. panicea"))+
 # labs(title = " ",
     #  x = "Vibrio strain",
     #  y = "% Phagocytosis") +
 # theme_minimal()
  
#print(self_plot)

#plot both time points as separate bars under one species
ank_plot_time <- ggplot(FACS, aes(
  x = vibrio_strain,
  y = phago_avg_DAPI_adjusted,
  colour = time_min
)) +
  geom_boxplot(
    aes(fill = factor(ankyrin_present)),
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    aes(fill = factor(ankyrin_present)),
    alpha = 0.6,
    size = 1.5,
    position = position_jitterdodge(
      dodge.width = 0.8,
      jitter.width = 0.15
    ),
    color = "black"
  ) +
  scale_fill_manual(
    name = "Ankyrin",
    values = c("positive" = "#CC6666", "negative" = "#9999CC"),
    labels = c("absent", "present")
  ) +
  scale_colour_manual(
    name = "Timepoint",
    values = c("#3b3b3b", "#000000"),
    labels = c("30 min", "60min")
  ) +
  labs(
    x = "Vibrio strain",
    y = "% Phagocytosis",
    title = ""
  ) +
  theme_minimal()

print(ank_plot_time)


###same as above but different time point is striped
library(ggpattern)
ank_plot_time2 <- ggplot(
  FACS,
  aes(
    x = vibrio_strain,
    y = phago_avg_DAPI_adjusted,
    #alpha = self_vs_foreign,
    group = interaction(vibrio_strain, time_min)
  )
) +
  geom_boxplot_pattern(
    aes(
      fill    = factor(ankyrin_present),
      pattern = time_min
    ),
    outlier.shape = NA,
    position = position_dodge(width = 0.8),
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.03
  ) +
  geom_jitter(
    aes(fill = factor(ankyrin_present)),
    size = 1,
    width = 0.15,
    #position = position_jitterdodge(
      #dodge.width = 0.8,
      #jitter.width = 0.15
    #),
    color = "black"
  ) +
  scale_fill_manual(
    name = "Ankyrin",
    values = c("positive" = "#CC6666", "negative" = "#9999CC"),
    labels = c("absent", "present")
  ) +
  scale_pattern_manual(
    name = "Timepoint",
    values = c("30" = "none", "60" = "stripe"),
    labels = c("30 min", "60 min")
  ) +
  #scale_alpha_manual(
    #name = "Origin",
    #values = c("self" = 0.4, "foreign" = 0.9)
  #) +
  labs(
    x = "Vibrio strain",
    y = "% Phagocytosis"
  ) +
  theme_minimal()+
  theme(axis.text=element_text(size=10.5),
        axis.title=element_text(size=14),
        axis.line.y = element_line(color = "grey", linewidth = 0.2))

print(ank_plot_time2)



###MODELLING - MULTIWAY ANOVA###
library(car)
#use Type III sums of squares because of inbalanced design
options(contrasts=c("contr.sum","contr.poly"))
#define Model
  #vibrio strains
Model1<-aov(FACS$phago_avg_DAPI_adjusted~FACS$vibrio_strain*FACS$time_min)
Model1.5<-aov(FACS$phago_avg_DAPI_adjusted~FACS$vibrio_strain)
  #ankyrin present + time
Model2<-aov(FACS$phago_avg_DAPI_adjusted~FACS$ankyrin_present*FACS$time_min)
  #self vs foreign + time
Model3<-aov(FACS$phago_avg_DAPI_adjusted~FACS$self_vs_foreign*FACS$time_min)
  #self vs foreign + ankyrin present + time
Model4<-aov(FACS$phago_avg_DAPI_adjusted~FACS$self_vs_foreign*FACS$time_min*FACS$ankyrin_present)
  #including background as a predictor
Model5<-aov(FACS$phago_avg_DAPI_adjusted~FACS$self_vs_foreign*FACS$time_min*FACS$ankyrin_present+FACS$background)
  #including initial Vibrio concentration as a predictor
Model6<-aov(FACS$phago_avg_DAPI_adjusted~FACS$self_vs_foreign*FACS$time_min*FACS$ankyrin_present+FACS$t0_vibrio)


#OUTPUT
Anova(Model1, type="III") #strain & time significant, effect of interaction also significant
Anova(Model1.5, type = "III") #strain significant
Anova(Model2, type="III") #time significant, ankyrin and interaction not
Anova(Model3, type="III") #self vs foreign signifcant & time, interaction not
Anova(Model4, type="III") #self vs foreign signifcant & time, some interactions
Anova(Model5, type="III") #background is not  significant
Anova(Model6, type="III") #initial concentration is not  significant
AIC(Model3,Model4) #model 4 which includes origin is better at describing the situation
AIC(Model5,Model4) #can't be done

#EFFECT SIZE
library(rstatix)
res.aov <- FACS %>% anova_test(phago_avg_DAPI_adjusted ~self_vs_foreign*ankyrin_present*time_min)
res.aov #decimal numbers can be translated to effect size in %

#EMMS to see how much % is the difference
library(emmeans)
emm <- emmeans(Model4, ~ self_vs_foreign)
emm_cond <- emmeans(Model4, ~ self_vs_foreign | time_min * ankyrin_present)
pairs(emm_cond)
emmeans(Model4, ~ time_min)
pairs(emmeans(Model4, ~ time_min), type = "response")

# 1. Get condition-specific means
emm_cond_time <- emmeans(Model4, ~ time_min | self_vs_foreign * ankyrin_present)
# 2. Convert to dataframe
emm_df <- as.data.frame(emm_cond_time)
# 3. Calculate % change
percent_by_condition <- emm_df %>%
  pivot_wider(
    id_cols = c(self_vs_foreign, ankyrin_present),
    names_from = time_min,
    values_from = emmean
  ) %>%
  mutate(percent_change = (`60` - `30`) / `30` * 100)

percent_by_condition
#phagocytosis change per strain (grrr SA48)
strain_changes <- FACS %>%
  group_by(vibrio_strain, time_min) %>%
  summarise(mean_phago = mean(phago_avg_DAPI_adjusted, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = time_min, values_from = mean_phago) %>%
  mutate(percent_change = (`60` - `30`) / `30` * 100)

strain_changes

#DIAGNOSTICS INTERLUDE#
#homogeneity of variances
plot(resid(Model4)~fitted(Model4))
abline(h=0, lwd=2, lty=2, col="black") #oof I don't know, its' supposed to look like a starry sky
fligner.test(FACS$phago_avg_DAPI_adjusted~FACS$self_vs_foreign) #insignificant - all is well
fligner.test(FACS$phago_avg_DAPI_adjusted~FACS$ankyrin_present) #insignificant - all is well
fligner.test(FACS$phago_avg_DAPI_adjusted~FACS$time_min) #insignificant - all is well
fligner.test(FACS$phago_avg_DAPI_adjusted~interaction(FACS$self_vs_foreign,FACS$ankyrin_present)) #insignificant - all is well
fligner.test(FACS$phago_avg_DAPI_adjusted~interaction(FACS$self_vs_foreign,FACS$time_min)) #insignificant - all is well
fligner.test(FACS$phago_avg_DAPI_adjusted~interaction(FACS$ankyrin_present,FACS$time_min)) #significant, so variances differ, but does it matter?
fligner.test(FACS$phago_avg_DAPI_adjusted~interaction(FACS$ankyrin_present,FACS$time_min, FACS$self_vs_foreign)) #insignificant - all is well
#normality of errors
hist(resid(Model4), breaks=20) #little wild but ok
qqnorm(resid(Model4))
qqline(resid(Model4)) #looks good
shapiro.test(resid(Model4)) #not significant, all is well
#influential data points
plot(cooks.distance(Model4), type="h") #looks good

#Tukex's test to check exavtly which pairs are different
tukey <- TukeyHSD(Model4)
print(tukey)

tukey_species <- TukeyHSD(Model1.5)
print(tukey_species)

tukey_species_time <- TukeyHSD(Model1)
print(tukey_species_time)

#compact letter display
library(multcompView)
cld <- multcompLetters4(Model1,tukey_species_time)
print(cld) #there are nice letters between different plots at the very end

cld_species <- multcompLetters4(Model1.5,tukey_species)
print(cld_species)

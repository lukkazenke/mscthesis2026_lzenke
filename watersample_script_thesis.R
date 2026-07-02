####### Water sample analysis - Lukka Zenke thesis - last updated 16 March 2026 #####
rm(list=ls(all=TRUE))
###LOAD DATA ###
water_uncut<-read.csv("MSc_water_rawdata.csv", header = TRUE)
water <- water_uncut[-c(78),] #cut rows with NA and also 

#loading controls separately
water_control <-read.csv("thesis_watersample_controls.csv", header = TRUE)
water_control$replicate <- as.factor(water_control$replicate)
#put the variable in the desired order bc R otherwise sorts alphabetically
water$vibrio_strain <- as.factor(water$vibrio_strain)
water$strain_replicate <- as.factor(water$strain_replicate)
water$ankyrin_present <- as.factor(water$ankyrin_present)
water$self_vs_foreign <- as.factor(water$self_vs_foreign)
water$time_min <- as.numeric(levels(water$time_min))[water$time_min]
water$sample_id <- as.factor(water$sample_id)
water$final_vibrio.mL <- as.integer(water$final_vibrio.mL)
water$norm_vibrio.mL <- as.integer(water$norm_vibrio.mL)
str(water)
unique(water$time_min)

### TIME MEAN ###
library(dplyr)

#normalized mean per replicate perr strain for exponential calculations
water_mean_replicate <- water %>%
  group_by(strain_replicate,time_min) %>%
  summarise(
    vibrio_mean = mean(norm_vibrio.mL, na.rm = TRUE),
    .groups = "drop"
  )

#normalized mean per strain for plotting
water_mean_vibrio <- water %>%
  group_by(vibrio_strain,time_min,ankyrin_present) %>%
  summarise(
    vibrio_mean = mean(norm_vibrio.mL, na.rm = TRUE),
    .groups = "drop"
  )


### PLOTTING ###
library(ggplot2)

#plotting controls
ggplot(water_control, aes(x = time_min, y = final_vibrio.mL, colour = replicate)) +
  geom_jitter(size = 1, width = 0.2) +
  #do we want a linear smooth, non-linear smooth or even a line?
  geom_smooth(
    method = "lm", span = 1.2, se = FALSE, linewidth = 0.5) +
  #scale_color_manual(
  #name = "Vibrio Strain",                     # legend title
  #values = c("Temp_ºC" = "#CC6666",
  # "PSU"    = "#9999CC",
  #"pH" = "black"),
  #labels = c("pH", "Salinity [PSU]", "Temperature [ºC]" )   # legend labels
  #) +
  labs(x = "Time [min]", y = "Vibrio/mL") +
  scale_x_continuous(breaks = c(0,10,20,30, 40, 50, 60)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme_minimal()

ggplot(water, aes(x = time_min, y = final_vibrio.mL, colour = strain_replicate)) +
  geom_jitter(size = 1, width = 0.2) +
  #do we want a linear smooth, non-linear smooth or even a line?
  geom_smooth(
    method = "loess", span = 1.2, se = FALSE, linewidth = 0.5) +
  #scale_color_manual(
  #name = "Vibrio Strain",                     # legend title
  #values = c("Temp_ºC" = "#CC6666",
  # "PSU"    = "#9999CC",
  #"pH" = "black"),
  #labels = c("pH", "Salinity [PSU]", "Temperature [ºC]" )   # legend labels
  #) +
  labs(x = "Time [min]", y = "Vibrio/mL") +
  scale_x_continuous(breaks = c(0,10,20,30, 40, 50, 60)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme_minimal()

##-----#plotting Vibrio strains#------##
#use semi-log scale like Riisgard et al (Denmark)
ggplot(water, aes(x = time_min, y = norm_vibrio.mL,
                  colour = vibrio_strain, linetype = ankyrin_present )) +
  #geom_jitter(size = 1, width = 0.2) +
  geom_point(
    data = water_mean_vibrio,
    aes(time_min, vibrio_mean),
    size = 1
  )+
  #do we want a linear smooth (lm), non-linear smooth (loess)?
  geom_smooth(
   data = water_mean_vibrio,
    aes(time_min, vibrio_mean),
    method = "lm", span = 1.2, se = FALSE, linewidth = 0.5) +
  #no smoothing just connecting
  #geom_line(
   # data = water_mean_vibrio,
    #aes(time_min, vibrio_mean),
  #)+
  labs(x = "Time [min]", y = "Vibrio/mL") +
  scale_x_continuous(breaks = c(0,10,20,30, 40, 50, 60)) +
  scale_y_continuous(trans='log10',
                     limits = c(400000,1500000),
                     labels = function(x) format(x, scientific = TRUE)) +
  #colourblind friendly palette
  scale_color_manual(
   name = "Vibrio strain",
    values = c("#0e288e","#164c64","#a53606", "#b32db5","#881a58"),
    labels = c("Hal054", "Hal281","M60_M31", "M60_M70","SA48"))+
  scale_linetype_discrete(name = "Ankyrin",
                          labels = c("absent", "present"))+
  theme_minimal()

### CALCULATING EXPONENTIAL CURVES FOR UPTAKE RATES####
#exponential fit for each vibrio_strain based on the vibrio cell mean per each time point
exp_stats <- water %>%
  group_by(strain_replicate) %>%
  group_modify(~{
    
    df <- .x
    
    fit <- nls(
      final_vibrio.mL ~ a * exp(b * time_min),
      data = df,
      start = list(
        a = min(df$final_vibrio.mL),
        b = 0.01
      )
    )
    
    y <- df$final_vibrio.mL
    y_hat <- predict(fit)
    
    r2 <- 1 - sum((y - y_hat)^2) / sum((y - mean(y))^2)
    
    tibble(
      a = coef(fit)["a"],
      b = coef(fit)["b"],
      R2 = r2
    )
  })

exp_stats #view R2

# Create a column with the actual exponential function as a string
exp_stats_func <- exp_stats %>%
  mutate(
    function_str = paste0("y = ", round(a, 3), " * exp(", round(b, 5), " * x)"),
    R2 = round(R2, 3)
  ) %>%
  select(strain_replicate, function_str, R2)

exp_stats_func #view table

#fit linear functions to controls
lin_stats <- water_control %>%
  group_by(replicate) %>%
  group_modify(~{
    
    df <- .x
    
    fit <- lm(final_vibrio.mL ~ time_min, data = df)
    
    y <- df$final_vibrio.mL
    y_hat <- predict(fit)
    
    r2 <- 1 - sum((y - y_hat)^2) / sum((y - mean(y))^2)
    
    tibble(
      a = coef(fit)[1],   # intercept
      b = coef(fit)[2],   # slope
      R2 = r2
    )
  })
#display as table
lin_stats_func <- lin_stats %>%
  mutate(
    function_str = paste0(
      "y = ",
      round(a, 3),
      " + ",
      round(b, 5),
      " * x"
    ),
    R2 = round(R2, 3)
  ) %>%
  select(replicate, function_str, R2)

lin_stats_func



##----------STATISTICS---------------##
w_stats<-read.csv("water_stats_finaldata.csv", header = TRUE)
w_stats <- w_stats[-c(1:6,13,14,19),] #cut controls and rows with negative filtration

w_stats$vibrio_strain <- as.factor(w_stats$vibrio_strain)
w_stats$ankyrin_present <- as.factor(w_stats$ankyrin_present)
w_stats$self_vs_foreign <- as.factor(w_stats$self_vs_foreign)
w_stats$time_min <- as.factor(w_stats$duraction_min)
w_stats$perc_change <- as.numeric(w_stats$perc_change)
str(w_stats)

#mean F/V per strain
aggregate(w_stats$filtration_mL_min.1_mL.1~w_stats$vibrio_strain, FUN = mean)
aggregate(w_stats$filtration_mL_min.1_mL.1~w_stats$vibrio_strain, FUN = sd)
#mean F/V overall
mean(w_stats$filtration_mL_min.1_mL.1)
sd(w_stats$filtration_mL_min.1_mL.1)
#mean intial Vibrio per strain 
aggregate(w_stats$Vibrio.mL_start~w_stats$vibrio_strain, FUN = mean)
aggregate(w_stats$Vibrio.mL_start~w_stats$vibrio_strain, FUN = sd)

##--PLOTTING--##
library(ggplot2)

#barplots for mL water filtration rate
filt_plot <- ggplot(w_stats, aes(x = vibrio_strain, y = filtration_mL_min.1_mL.1,
                                   fill= ankyrin_present))+
  geom_boxplot()+
  scale_fill_manual(
    name = "Ankyrin",
    values = c("yes" = "#CC6666", "no"="#9999CC"),
    labels = c("absent", "present")
  ) +
  geom_jitter(
    aes(fill = factor(self_vs_foreign)),
    width = 0.05, size = 1.5,
    colour = "black")+
  labs(x = "Vibrio strain", y = "Filtration rate [mL*min-1*mL-1]") +
  scale_y_continuous(limits = c(0,5))+
  theme_minimal()

print(filt_plot)


##--STATISTICS--##
library(car)
#1) check normality & outliers
dotchart(w_stats$filtration_mL_min.1_mL.1)
hist(w_stats$filtration_mL_min.1_mL.1) #looks ok!
qqnorm(w_stats$filtration_mL_min.1_mL.1) #doesn't look good
qqline(w_stats$filtration_mL_min.1_mL.1)
shapiro.test(w_stats$filtration_mL_min.1_mL.1) #insignificant, all is well
#2)check for balanced design
table(w_stats$vibrio_strain) #unbalance because some had to be excluded
#3)check for homogeniety of variances
leveneTest(w_stats$filtration_mL_min.1_mL.1~w_stats$vibrio_strain, data = w_stats) 
#variances are OK!


##ANOVA MODELLING##
options(contrasts=c("contr.sum","contr.poly"))
#vibrio strains
Model0<- aov(w_stats$filtration_mL_min.1_mL.1~w_stats$vibrio_strain)
#only-ankyrin
Model1<- aov(w_stats$filtration_mL_min.1_mL.1~w_stats$ankyrin_present)
#ankyrin + starting concentration
Model1.5<- aov(w_stats$filtration_mL_min.1_mL.1~w_stats$ankyrin_present*w_stats$Vibrio.mL_start)
#also checking starting concentration if that is messing something up
Model2 <- aov(w_stats$filtration_mL_min.1_mL.1~w_stats$Vibrio.mL_start)
#ankyrin + self vs foreign
Model3 <- aov(w_stats$filtration_mL_min.1_mL.1~w_stats$ankyrin_present*w_stats$self_vs_foreign)
#ankyrin + self vs foreign + background + starting concentration + size
Model4 <- aov(w_stats$filtration_mL_min.1_mL.1~w_stats$ankyrin_present*w_stats$self_vs_foreign
              + w_stats$background + w_stats$sponge_weight_g + w_stats$Vibrio.mL_start)

#OUTPUT
#use Type III sums of squares because of inbalanced design
Anova(Model0, type="III") #vibrio strain overall not significant
Anova(Model1, type="III") #ankyrin is significant
Anova(Model1.5, type="III") #ankyrin not significant anymore, starting also not significant
Anova(Model2, type="III") #starting concentration isn't
Anova(Model3, type="III") #ankyrin remains significant, origin is not
Anova(Model4, type="III") #initial concentration and size significant!!

AIC(Model1, Model3) #ank-only model is best
AIC(Model4, Model3) #model without background is better
AIC(Model1, Model4) #model with size and initial concentration is best

#EFFECT SIZE
library(rstatix)
res.aov.water <- w_stats %>% anova_test(filtration_mL_min.1_mL.1 ~ ankyrin_present*self_vs_foreign
                                        + background + sponge_weight_g + Vibrio.mL_start)
res.aov.water #ges = generalized eta squared =  effect size


##DIAGNOSTICS INTERLUDE#
#homogeneity of variances
plot(resid(Model4)~fitted(Model4))
abline(h=0, lwd=2, lty=2, col="black") #good!
#normality of errors
hist(resid(Model4), breaks=20) #little wild but ok
qqnorm(resid(Model4))
qqline(resid(Model4)) #looks good
shapiro.test(resid(Model4)) #not significant, all is well
#influential data points
plot(cooks.distance(Model4), type="h") #looks good

#how much does it change?
model_log <- lm(log(filtration_mL_min.1_mL.1) ~ ankyrin_present*self_vs_foreign +
                  background + sponge_weight_g + Vibrio.mL_start,
                data = w_stats)
coef(model_log)["sponge_weight_g"] #change decimal into percent

##plot relationship of weight and filtration rate
library(ggplot2)
library(ggpmisc) # for showing equation and R2 in plot (https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph)

ggplot(w_stats, aes(x = sponge_vol_mL, y = filtration_mL_min.1_mL.1)) +
  geom_point( size = 2, aes (colour = factor(ankyrin_present)))+
  scale_colour_manual(name = "Ankyrin",
     values = c("yes" = "#CC6666", "no"="#9999CC"),
    labels = c("absent", "present"))+
 geom_smooth(method = "lm", colour = "black", se = FALSE)+
   #stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  labs(x = "sponge volume [mL]", y = "Filtration rate [mL*min-1*mL-1]") +
  scale_y_continuous(limits = c(0,4))+
  scale_x_continuous(limits = c(0,10))+
  theme_minimal()

#plot that compares linear model to Riisgard 2022
model_power <- lm(log(filtration_mL_min.1_mL.1) ~ log(sponge_vol_mL),
                  data = w_stats)

a <- exp(coef(model_power)[1])
b <- coef(model_power)[2]

r2_power <- summary(model_power)$r.squared
ggplot(w_stats, aes(x = sponge_vol_mL, y = filtration_mL_min.1_mL.1)) +
  geom_point(size = 2, aes(colour = factor(ankyrin_present))) +
  scale_colour_manual(name = "Ankyrin",
                      values = c("yes" = "#CC6666", "no"="#9999CC"),
                      labels = c("absent", "present")) +
  # Linear model
  stat_poly_line(colour = "black", se = FALSE) +
  #stat_poly_eq(use_label(c("eq", "R2"))) +
  # Power model curve based on Riisgard
  geom_function(fun = function(x) a * x^b,
                colour = "blue", linewidth = 1,
                xlim = range(w_stats$sponge_weight_g)) +
  # Power model label
  #annotate("text",
           #x = 7, y = 3.5,
           #label = paste0("y = ", round(a, 2), "x^", round(b, 2),
                          #"\nR² = ", round(r2_power, 2)),
          # hjust = 0, colour = "blue") +
  #layout
  labs(x = "sponge volume [mL]",
       y = "F/V [mL*min-1*mL-1]") +
  scale_y_continuous(limits = c(0,5)) +
  scale_x_continuous(limits = c(2,10)) +
  theme_minimal()

##relationship of intial vibrio concentration and filtration
plot(w_stats$Vibrio.mL_start, w_stats$filtration_mL_min.1_mL.1)
lines(lowess(w_stats$Vibrio.mL_start, w_stats$filtration_mL_min.1_mL.1), col="red")
  #shows: not linear but rather messy, highly sponge dependent

#########################################
####percentage change of Vibrios########
#MEANS AND SDs
aggregate(w_stats$perc_change ~ w_stats$vibrio_strain, FUN = mean)
aggregate(w_stats$perc_change ~ w_stats$vibrio_strain, FUN = sd)

##--STATISTICS--##
library(car)
#1) check normality & outliers
dotchart(w_stats$perc_change) #no outliers
hist(w_stats$perc_change) #opposite of normal
qqnorm(w_stats$perc_change) #doesn't look good
qqline(w_stats$perc_change)
shapiro.test(w_stats$perc_change) #significant, not normal
#2)check for balanced design
table(w_stats$vibrio_strain) #unbalance because some had to be excluded
#3)check for homogeniety of variances
leveneTest(w_stats$perc_change~w_stats$vibrio_strain, data = w_stats) 
#test not significant = variances are OK!

#MODELLING
#we are starting with linear model
model <- lm(perc_change ~ ankyrin_present * self_vs_foreign, data = w_stats)
plot(model)
model2 <- lm(perc_change ~ vibrio_strain, data = w_stats)
summary(model2)
#using permutation based ANOVA to account for deviations from normality
library(permuco)
perm_model1 <- aovperm(perc_change ~ ankyrin_present*self_vs_foreign, data = w_stats, np = 5000)
summary(perm_model1) #nothing significant
#check vibrio strains individually
perm_model2 <- aovperm(perc_change ~ vibrio_strain, data = w_stats, np = 5000)
summary(perm_model2) #not significant
#check background
perm_model3 <- aovperm(perc_change ~ ankyrin_present*self_vs_foreign + background, data = w_stats, np = 5000)
summary(perm_model3)
#with time_min
perm_model4 <- aovperm(perc_change ~ ankyrin_present*self_vs_foreign*time_min, data = w_stats, np = 5000)
summary(perm_model4) #significant to noones surprise

#################################################
### CHECK UNSTAINED PORTION OF THE WATER
###LOAD DATA ###
water_unstained<-read.csv("watersample_data_unstained.csv", header = TRUE)
water_unstained$vibrio_strain <-as.factor(water_unstained$vibrio_strain )
water_unstained$date <-as.factor(water_unstained$date)
water_unstained$Total_calc <-as.integer(water_unstained$Total_calc)
water_unstained$unstained_calc <-as.integer(water_unstained$unstained_calc)
water_unstained$ankyrin_present <- as.factor(water_unstained$ankyrin_present)
str(water_unstained)

unstained <- water_unstained[-c(210:251),] 
unstained$vibrio_strain <-as.factor(unstained$vibrio_strain)
str(unstained)
##boxplot##
library(ggplot2)
box_plot <- ggplot(unstained, aes(x = vibrio_strain, y = unstained_calc,
                                 ))+
  geom_boxplot()+
  theme_minimal()
print(box_plot)

#statistics#
library(car)
#1) check normality & outliers
dotchart(water_unstained$unstained_calc)
hist(water_unstained$unstained_calc) #doesn't look good
qqnorm(water_unstained$unstained_calc) #doesn't look good
qqline(water_unstained$unstained_calc)
shapiro.test(water_unstained$unstained_calc) #significant, not normal!
#2)check for homogeniety of variances
leveneTest(water_unstained$unstained_calc~water_unstained$date, data = water_unstained) 
#variances are NOT OK!

#not normal + wild variances = Kruskal Wallis Test
kruskal.test(water_unstained$unstained_calc~water_unstained$date)
#significant differences! probably because of M31 (19.08.)
#checking pair-wise
pairwise.wilcox.test(water_unstained$unstained_calc,water_unstained$date,
                     p.adjust = "bonferroni")

#check for vibrio strains instead of dates
kruskal.test(unstained$unstained_calc~unstained$vibrio_strain)
  #highly significant
#checking pair-wise
pairwise.wilcox.test(unstained$unstained_calc,unstained$vibrio_strain,
                     p.adjust = "bonferroni")
#all vibrio strains are significantly different to each other, except Hal281 & Hal054

#check between ank vs non-ank
kruskal.test(water_unstained$unstained_calc~water_unstained$ankyrin_present)
#significant differences!

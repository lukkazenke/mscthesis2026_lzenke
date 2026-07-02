####### Phagocytic cells and Vibrio counts - Lukka Zenke thesis - last updated 24.02.2026 #########
rm(list=ls(all=TRUE))
###LOAD DATA ###
vib_counts_uncut <- read.csv("vibrio_count_data.csv", header = TRUE)
#vib_counts <- vib_counts_uncut[-c(541:900),] 
vib_counts<-subset(vib_counts_uncut,vibrio_no!="NA") #cut rows with empty cells

#make sure everything has the right format
vib_counts$vibrio_strain <- as.factor(vib_counts$vibrio_strain)
vib_counts$ankyrin_present <- as.factor(vib_counts$ankyrin_present)
vib_counts$self_vs_foreign <- trimws(vib_counts$self_vs_foreign)
vib_counts$self_vs_foreign <- as.factor(vib_counts$self_vs_foreign)
vib_counts$time_min <- as.factor(vib_counts$time_min)
vib_counts$cell_type <- as.factor(vib_counts$cell_type)
vib_counts$vibrio_no <- as.numeric(vib_counts$vibrio_no )
str(vib_counts)

###-----VIBRIO CELLS PER CELL-----####
dotchart(vib_counts$vibrio_no)
hist(vib_counts$vibrio_no) #count data --> poisson
boxplot(vib_counts$vibrio_no~vib_counts$ankyrin_present)

###GLM###
#check underdispersion
mean(vib_counts$vibrio_no)
var(vib_counts$vibrio_no) #variance is much smaller than mean --> underdispersion!!

#using COM-Poisson for
library(glmmTMB)
model_cmp <- glmmTMB(vibrio_no ~ ankyrin_present * time_min * cell_type,
                     family = compois(),
                     data = vib_counts)
summary(model_cmp) 
#keep in mind that summary uses a level of a group the intercept 
#and then compares the other levels to that group
car::Anova(model_cmp, type = "III")
#three-way interaction is not significant, let's scrub it
model_cmp_no3 <- glmmTMB(vibrio_no ~ ankyrin_present * time_min +
                        ankyrin_present * cell_type +
                       time_min * cell_type,
                     family = compois(),
                     data = vib_counts)
summary(model_cmp_no3) 
#keep in mind that summary uses a level of a group the intercept 
#and then compares the other levels to that group
car::Anova(model_cmp_no3, type = "III")

###no interaction model
model_final <- glmmTMB(vibrio_no ~ ankyrin_present + time_min + cell_type + self_vs_foreign,
                   family = compois(),
                   data = vib_counts)
summary(model_final) 
car::Anova(model_final, type = "III")

AIC(model_cmp_no3, model_final) #no interaction is better!! 
 
#there is no difference between vibrios per cell type between treatments

#checking if SW background is significant
model_bckgr <- glmmTMB(vibrio_no ~ ankyrin_present + time_min + cell_type + self_vs_foreign + background,
                       family = compois(),
                       data = vib_counts)
summary(model_bckgr) 
car::Anova(model_bckgr, type = "III")

AIC(model_bckgr, model_final)

#compare regular poisson just to be sure
model_pois <- glmmTMB(vibrio_no ~ ankyrin_present * time_min,
                      family = poisson,
                      data = vib_counts)

AIC(model_pois, model_final) #cmp is WAY better

#Diagnostics interlude
library(DHARMa)
sim_res <- simulateResiduals(model_final)
plot(sim_res)
testDispersion(sim_res) #disperion parameter close to 1 = good!, p>0.05 = good!
model_cmp$sdr$pdHess #makes TRUE -> good!
sigma(model_final)

###effect sizes###
#1) Incidence Rate Ratios (IRRs)with confidence intervals
exp(confint(model_final, parm = "beta_"))

#2) Pairwise contrasts
library(emmeans)
emm <- emmeans(model_final, ~ ankyrin_present) #ankyrin
pairs(emm, type = "response")
pairs(emmeans(model_final, ~ cell_type), type = "response") #cell type
#L is different than the others, but not the others among each other
pairs(emmeans(model_final, ~ time_min), type = "response")


####---- PROPORTIONS OF CELL TYPES ---#####

#look at my data -> note that foreign/ankyrin/60min is missing!
ftable(vib_counts$self_vs_foreign,
      vib_counts$time_min,
      vib_counts$ankyrin_present)
ftable(vib_counts$vibrio_strain,
       vib_counts$time_min)

##PLOTTING##
library(ggplot2)
ggplot(vib_counts, aes(x =vibrio_strain, fill = cell_type))+
  geom_bar (position = "fill")+
  labs(
    x = "Vibrio strain",
    y = "Percentage of phagocytic cells",
    fill = "Cell types"
  ) +
  scale_fill_manual(
    name = "Cell type",
values = c("F+" = "#F4F4F4", "L" = "#A4A4A4", "M" = "#555555", "S" = "#050505"),
labels = c("flagellated", "large", "medium", "small")
) +
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()


##STATISTICS##
#response variable is categorical! => Multinomical logistic regression 
#just using glm -> no three way interaction because of missing data points
  #vibrio strain
model_cells_glm_strain <- glm(cell_type ~ vibrio_strain
                               + time_min,
                       family = "binomial", data = vib_counts)

summary(model_cells_glm_strain)
car::Anova(model_cells_glm_strain, type = "III")

#all the other factors
model_cells_glm <- glm(cell_type ~ ankyrin_present * time_min +
                         ankyrin_present * self_vs_foreign +
                         time_min * self_vs_foreign,
                       family = "binomial", data = vib_counts)

summary(model_cells_glm)
car::Anova(model_cells_glm, type = "III")
  #glm without interactions
model_cells_glm2 <- glm(cell_type ~ ankyrin_present + time_min + self_vs_foreign,
                       family = "binomial", data = vib_counts)

summary(model_cells_glm2)

AIC(model_cells_glm, model_cells_glm2) #with interactions is better

#DIAGNOSTICS INTERLUDE - GLM this time#
#homogeneity of variances
plot(resid(model_cells_glm)~fitted(model_cells_glm))
abline(h=0, lwd=2, lty=2, col="black") #performs better than gaussian model, no trend
#normality of errors
hist(resid(model_cells_glm), breaks=25) #two peaks that don't look nice
qqnorm(resid(model_cells_glm))
qqline(resid(model_cells_glm)) #qqplot slightly worse than Gaussian, but still ok
shapiro.test(resid(model_cells_glm)) #non significant - all is well
#influential data points
plot(cooks.distance(model_cells_glm), type="h")



####Vibrio cells per cell type######

#####PLOTTING####
#1)make two subsets for with and without ankyrin
library(dplyr)
df_ank <- vib_counts %>% filter(ankyrin_present == "yes")
df_noank <- vib_counts %>% filter(ankyrin_present == "no")
#2)make vibrio numbers a factor for plotting, sort cell types
df_ank$cell_type <- factor(df_ank$cell_type, c("F+","S","M","L"))
df_ank$vibrio_no <- as.factor(df_ank$vibrio_no)
df_noank$cell_type <- factor(df_noank$cell_type, c("F+","S","M","L"))
df_noank$vibrio_no <- as.factor(df_noank$vibrio_no)
#3) calculate total percentage for each cell type
df_prop_ank <- df_ank %>%
  count(cell_type, vibrio_no) %>%  # no time_min
  group_by(cell_type) %>%
  mutate(prop = n / sum(n))

df_prop_ank_time <- df_ank %>%
  count(cell_type, vibrio_no, time_min) %>% 
  group_by(time_min) %>% 
  mutate(prop = n / sum(n))

df_prop_no <- df_noank %>%
  count(cell_type, vibrio_no) %>%
  ungroup() %>%  # important: remove grouping
  mutate(prop = n / sum(n))  # now prop is fraction of all cells

df_prop_no_time <- df_noank %>%
  count(cell_type, vibrio_no,time_min) %>%  # with time_min
  group_by(time_min) %>%
  mutate(prop = n / sum(n))

#4a)plot ANK
ggplot(df_prop_ank_time,
       aes(x = cell_type,
           y = prop,
           fill = vibrio_no)) +
  geom_col(position = position_stack(),
           width = 0.8) +
  facet_wrap(~ time_min,
             labeller = labeller(time_min = function(x) paste0(x, "min"))) +
  labs(
    x = "Cell type",
    y = "Percentage of phagocytic cells",
    fill = "No. of Vibrio"
  ) +
  scale_fill_manual(values = c("#F2DADA","#E6A8A8","#CC6666","#A94444"
  )) +
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()

#4a)plot NO-ANK
ggplot(df_prop_no_time,
       aes(x = cell_type,
           y = prop,
           fill = vibrio_no)) +
  geom_col(position = position_stack(),
           width = 0.8) +
  facet_wrap(~ time_min,
             labeller = labeller(time_min = function(x) paste0(x, "min"))) +
  labs(
    x = "Cell type",
    y = "Percentage of phagocytic cells",
    fill = "No. of Vibrio"
  ) +
  scale_fill_manual(values = c("#E5E5F2","#C2C2E6","#9999CC","#6F6FB3")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()

#############
#plot that is seperated by ANK presence and not by time point
##################
df_prop3 <- vib_counts %>%
  count(ankyrin_present, cell_type, vibrio_no) %>%
  group_by(ankyrin_present) %>%
  mutate(prop = n / sum(n))

df_prop3$vibrio_no <- as.factor(df_prop3$vibrio_no)
df_prop3$cell_type <- factor(df_prop3$cell_type, c("F+","S","M","L"))

ggplot(df_prop3,
       aes(x = cell_type,
           y = prop,
           fill = vibrio_no)) +
  geom_col(position = position_stack(),
           width = 0.8) +
  facet_wrap(~ ankyrin_present)+
  labs(
    x = "Cell type",
    y = "Percentage of phagocytic cells",
    fill = expression(italic("Vibrio"))
  ) +
  scale_fill_manual(values = c("#E5E5F2","#C2C2E6","#9999CC","#6F6FB3"
  )) +
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()


## INDIVIDUAL PLOTS FOR ALL STRAINS WITH TIME###
  
#HAL281
#1)make separate dfs
df_Hal281 <- subset(vib_counts, vibrio_strain == "Hal281")
df_Hal281_prop <- df_Hal281 %>%
  count(cell_type, vibrio_no, time_min) %>%
  group_by(time_min) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(vibrio_no = factor(vibrio_no)) #make sure Vib count is a factor
df_Hal281_prop$cell_type <- factor(df_Hal281_prop$cell_type, c("F+","S","M","L"))

#2)plot, colours based on ANK or not
ggplot(df_Hal281_prop,
       aes(x = cell_type,
           y = prop,
           fill = vibrio_no)) +
  geom_col(position = position_stack(),
           width = 0.8) +
  facet_wrap(~ time_min,
             labeller = labeller(time_min = function(x) paste0(x, "min"))) +
  labs(
    x = "Cell type",
    y = "% of phagocytic cells",
    fill = "No. of Vibrio"
  ) +
  scale_fill_manual(values = c("#E6D5FF", "#A855F7", "#6A00B8", "#240046")) +
  #scale_colour_manual(values = c("#4DA6FF", "#006DCC", "#003F7F", "#002855")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+
  theme(axis.text=element_text(size=10.5),
        axis.title=element_text(size=14),
        axis.line.y = element_line(color = "grey", linewidth = 0.2))

#HAL054
#1)make separate dfs
df_Hal054 <- subset(vib_counts, vibrio_strain == "Hal054")
df_Hal054_prop <- df_Hal054 %>%
  count(cell_type, vibrio_no, time_min) %>%
  group_by(time_min) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(vibrio_no = factor(vibrio_no)) #make sure Vib count is a factor
df_Hal054_prop$cell_type <- factor(df_Hal054_prop$cell_type, c("F+","S","M","L"))

#2)plot, colours based on ANK or not
ggplot(df_Hal054_prop,
       aes(x = cell_type,
           y = prop,
           fill = vibrio_no)) +
  geom_col(position = position_stack(),
           width = 0.8) +
  facet_wrap(~ time_min,
             labeller = labeller(time_min = function(x) paste0(x, "min"))) +
  labs(
    x = "Cell type",
    y = "% of phagocytic cells",
    fill = "No. of Vibrio"
  ) +
  scale_fill_manual(values = c("#FFD6D6", "#FF5C5C", "#B30000", "#330000")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+
  theme(axis.text=element_text(size=10.5),
        axis.title=element_text(size=14),
        axis.line.y = element_line(color = "grey", linewidth = 0.2))

#SA48
#1)make separate dfs
df_SA48 <- subset(vib_counts, vibrio_strain == "SA48")
df_SA48_prop <- df_SA48 %>%
  count(cell_type, vibrio_no, time_min) %>%
  group_by(time_min) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(vibrio_no = factor(vibrio_no)) #make sure Vib count is a factor
df_SA48_prop$cell_type <- factor(df_SA48_prop$cell_type, c("F+","S","M","L"))

#2)plot, colours based on ANK or not
ggplot(df_SA48_prop,
       aes(x = cell_type,
           y = prop,
           fill = vibrio_no)) +
  geom_col(position = position_stack(),
           width = 0.8) +
  facet_wrap(~ time_min,
             labeller = labeller(time_min = function(x) paste0(x, "min"))) +
  labs(
    x = "Cell type",
    y = "% of phagocytic cells",
    fill = "No. of Vibrio"
  ) +
  scale_fill_manual(values = c("#E6D5FF", "#A855F7", "#6A00B8", "#240046")) +
  #scale_colour_manual(values = c("#4DA6FF", "#006DCC", "#003F7F", "#002855")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+
  theme(axis.text=element_text(size=10.5),
        axis.title=element_text(size=14),
        axis.line.y = element_line(color = "grey", linewidth = 0.2))

#M31
#1)make separate dfs
df_M31 <- subset(vib_counts, vibrio_strain == "M60_M31")
df_M31_prop <- df_M31 %>%
  count(cell_type, vibrio_no, time_min) %>%
  group_by(time_min) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(vibrio_no = factor(vibrio_no)) #make sure Vib count is a factor
df_M31_prop$cell_type <- factor(df_M31_prop$cell_type, c("F+","S","M","L"))

#2)plot, colours based on ANK or not
ggplot(df_M31_prop,
       aes(x = cell_type,
           y = prop,
           fill = vibrio_no)) +
  geom_col(position = position_stack(),
           width = 0.8) +
  facet_wrap(~ time_min,
             labeller = labeller(time_min = function(x) paste0(x, "min"))) +
  labs(
    x = "Cell type",
    y = "% of phagocytic cells",
    fill = "No. of Vibrio"
  ) +
  scale_fill_manual(values = c("#E6D5FF", "#A855F7", "#6A00B8", "#240046")) +
  #scale_colour_manual(values = c("#4DA6FF", "#006DCC", "#003F7F", "#002855")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+
  theme(axis.text=element_text(size=10.5),
        axis.title=element_text(size=14),
        axis.line.y = element_line(color = "grey", linewidth = 0.2))

#M70
#1)make separate dfs
df_M70 <- subset(vib_counts, vibrio_strain == "M60_M70")
df_M70_prop <- df_M70 %>%
  count(cell_type, vibrio_no, time_min) %>%
  group_by(time_min) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(vibrio_no = factor(vibrio_no)) #make sure Vib count is a factor
df_M70_prop$cell_type <- factor(df_M70_prop$cell_type, c("F+","S","M","L"))

#2)plot, colours based on ANK or not
ggplot(df_M70_prop,
       aes(x = cell_type,
           y = prop,
           fill = vibrio_no)) +
  geom_col(position = position_stack(),
           width = 0.8) +
  facet_wrap(~ time_min,
             labeller = labeller(time_min = function(x) paste0(x, "min"))) +
  labs(
    x = "Cell type",
    y = "% of phagocytic cells",
    fill = "No. of Vibrio"
  ) +
  scale_fill_manual(values = c("#FFD6D6", "#FF5C5C", "#B30000", "#330000")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+
  theme(axis.text=element_text(size=10.5),
        axis.title=element_text(size=14),
        axis.line.y = element_line(color = "grey", linewidth = 0.2))


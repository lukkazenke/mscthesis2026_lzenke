######## MSc thesis - environmental data analysis - 14. Nov 2025 #############
rm(list=ls(all=TRUE))
#import data
df1 <- read.csv("MSc_environmental_data.csv", header = TRUE)
ENV <- df1[-8,]

ENV$Strain<-as.factor(ENV$Strain)
ENV$Ankyrin.<-as.factor(ENV$Ankyrin.)
ENV$Date<-as.Date(ENV$Date, format="%Y-%m-%d")
str(ENV)

###AVERAGES###
mean(ENV$Temp_ºC) #18.34
sd(ENV$Temp_ºC) #0.70
mean(ENV$PSU) #14.56
sd(ENV$PSU) #1.35
mean(ENV$pH) #8.21
sd(ENV$pH) #0.15

###PLOTTING####
library(ggplot2)
library(tidyr)
library(dplyr)

df_long <- ENV %>% #fit all three parameters per date in one column
  pivot_longer(cols = c(Temp_ºC, PSU, pH),
               names_to = "variable",
               values_to = "value")
# Plot
ggplot(df_long, aes(x = Date, y = value, color = variable)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1)+
  scale_color_manual(
    name = "Parameter",                     # legend title
    values = c("Temp_ºC" = "#CC6666",
               "PSU"    = "#9999CC",
               "pH" = "black"),
    labels = c("pH", "Salinity [PSU]", "Temperature [ºC]" )   # legend labels
  ) +
  scale_y_continuous(limits = c(0, 25)) + 
  labs(x = " ", y = " ") +
  theme_minimal()


####STATISTICAL ANALYSIS####
#check if there were significant differences in abiotic factors between ankyrin positive and negative days
#statistical differences between all groups not possible because we only have 1 time poiint per group

hist(ENV$Temp_ºC, breaks= 6) #not normal
shapiro.test(ENV$Temp_ºC) #normal
hist(ENV$PSU) #normal
shapiro.test(ENV$PSU) #normal
hist(ENV$pH) #not normal
shapiro.test(ENV$pH) #not normal

t.test(ENV$Temp_ºC~ENV$Ankyrin) #no difference
kruskal.test(ENV$Temp_ºC~ENV$Ankyrin)#no difference

t.test(ENV$PSU~ENV$Ankyrin) #no difference
kruskal.test(ENV$PSU~ENV$Ankyrin)#no difference

t.test(ENV$pH~ENV$Ankyrin) #no difference
kruskal.test(ENV$pH~ENV$Ankyrin)#no difference

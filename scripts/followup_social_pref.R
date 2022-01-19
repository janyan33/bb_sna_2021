## setwd("C:/Users/janya/Desktop/R/bedbugs")

library(tidyverse)
library(ggplot2); theme_set(theme_classic())
library(lme4)
library(glmmTMB)
library(emmeans)

## Import and clean dataset
soc_pref <- read.csv("data/bb_soc_pref.csv", stringsAsFactors = TRUE) %>% 
            rename(sex = focal)

soc_pref$replicate <- as.factor(soc_pref$replicate)

## Removing NAs
soc_pref[soc_pref == "N/A"] <- NA
nrow(soc_pref[is.na(soc_pref$final_decision_side),])
soc_pref[is.na(soc_pref$final_decision_side),]

soc_pref <- soc_pref[!(is.na(soc_pref$final_decision_scent)),]

## Creating separate data sets for the 4 treatments (2 controls combined)
data_con <- soc_pref %>% 
            filter(treatment == "a" | treatment == "b")

data_con$final_decision_scent <- factor(data_con$final_decision_scent, levels = c("Mated F", "Control"))

data_male_choice <- soc_pref %>% 
                    filter(treatment == "c")
data_male_choice <- droplevels(data_male_choice)

data_female_choice <- soc_pref %>% 
                      filter(treatment == "d")
data_female_choice <- droplevels(data_female_choice)

data_virgins <- soc_pref %>% 
                filter(treatment == "e")
data_virgins <- droplevels(data_virgins)

## Plotting results
sex_labs <- c("Female focals", "Male focals")
names(sex_labs) <- c("F", "M")

## SOCIAL PREFERENCE BAR GRAPHS
## Treatments A and B: Controls
ggplot(data = data_con, aes(x = final_decision_scent, fill = sex)) + geom_bar(position = "dodge") + 
       facet_grid(~sex, labeller = labeller(sex = sex_labs)) + 
       theme(legend.position = "none", text = element_text(size = 18)) + 
       labs(y = "Count", x = "") + ylim(0, 30) + scale_fill_manual(values = c("sandybrown", "skyblue3"))

## Treatment C: Male focals, mated f vs. mated m scent
ggplot(data = data_male_choice, aes(x = final_decision_scent)) + geom_bar(position = "dodge", fill = "skyblue3") + 
       facet_grid(~sex, labeller = labeller(sex = sex_labs)) + 
       theme(legend.position = "none", text = element_text(size = 18)) + 
       labs(y = "Count", x = "") + ylim(0, 30)

## Treatment D: Female focals, mated f vs. mated m scent
ggplot(data = data_female_choice, aes(x = final_decision_scent)) + geom_bar(position = "dodge", fill = "sandybrown") + 
       facet_grid(~sex, labeller = labeller(sex = sex_labs)) + 
       theme(legend.position = "none", text = element_text(size = 18)) + 
       labs(y = "Count", x = "") + ylim(0, 30)

## Treatment E: Male focals, mated f vs. virgin f scent
ggplot(data = data_virgins, aes(x = final_decision_scent)) + geom_bar(position = "dodge", fill = "skyblue3") + 
      facet_grid(~sex, labeller = labeller(sex = sex_labs)) + 
      theme(legend.position = "none", text = element_text(size = 18)) + 
      labs(y = "Count", x = "") + ylim(0, 30)

## Models
con_model <- glmer(final_decision_scent ~ sex + (1|replicate), family = "binomial", 
                 data = data_con)
summary(con_model)

# TREATMENTS A AND B STATS
tab_con <- table(data_con$sex, data_con$final_decision_scent)
chisq.test(tab_con) #need to fix this one I think

con_model <- glmer(final_decision_scent ~ sex + (1|replicate), family = "binomial", 
                   data = data_con)
summary(con_model)

# TREATMENT C STATS
(tab_c <- table(data_male_choice$final_decision_scent))
chisq.test(tab_c)

male_model <- glmer(final_decision_scent ~ (1|replicate), family = "binomial", data = data_male_choice)
summary(male_model)

# TREATMENT D STATS
(tab_d <- table(data_female_choice$final_decision_scent))
chisq.test(tab_d)

female_model <- glmer(final_decision_scent ~ (1|replicate), family = "binomial", data = data_female_choice)
summary(female_model)

# TREATMENT E STATS
(tab_e <- table(data_virgins$final_decision_scent))
chisq.test(tab_e)

virgin_model <- glmer(final_decision_scent ~ (1|replicate), family = "binomial", data = data_virgins)
summary(virgin_model)

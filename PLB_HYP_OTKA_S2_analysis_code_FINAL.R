

#######################################
#  PLB-HYP OTKA S2 analysis code plan #
#######################################

#######################################
#              Packages               #
#######################################

library(lme4) # for lmer()
library(tidyverse) # for tidy code
library(gsheet)
library(lattice) # for qqmath()
library(lmeresampler) # for bootstrap.merMod() and confint.lmeresamp()

#######################################
#         Custom functions            #
#######################################

# convenience functions may be added here

FPQ_recoder <- function(varname){
  as.numeric(recode(varname,
                    "Extreme" = "5", 
                    "Very Much" = "4",
                    "A Fair Amount" = "3",
                    "A Little" = "2",
                    "Not At All" = "1",
                    .default = NA_character_))
}




#######################################
#            Data management          #
#######################################

data_pre_pre = read.csv("https://raw.githubusercontent.com/kekecsz/PLB-HYP-OTKA-S2/main/data/PLB-HYP-OTKA-S2-data-raw.csv")


# eligibility check

# general study eligibility:

## Those are eligible to participate in this study who:
## Are at least 18 years old
## Have not participated in this study before
## Have not tried hypnosis before
## Have not visited a university course on hypnosis before
## Do not have a mental disability
## Have not been diagnosed with schizophrenia or other forms of psychosis and who do not have a mental illness with symptoms of delusions or paranoia
## Have not been diagnosed with any form of epilepsy
## Have no current injury, skin disease, or frostbite on the arms
## Have no cardiovascular disease or Raynaud's phenomenon
## Have no history of fainting easily
## Are not suffering of chronic pain
## Are not suffering from substance use disorder
## Are not under the influence of pain medication or any other substance that affects state of consciousness or pain


# not test sessions
data_nontest = data_pre_pre %>% 
  filter(live.test..1.0. == 1)

# Some columns contain "," as a decimal instead of ".", we replace "," with "." to fix this.
cols_to_numerize = c("X1_paintolerance", "paintolerance_BL", "X2_paintolerance", "X1_painintensity", "X2_painintensity", "painintensity_BL", "watertemp_BL", "X1_watertemp", "X2_watertemp")
data_nontest[,cols_to_numerize] = sapply(data_nontest[,cols_to_numerize], function(x) gsub(",", ".", x))
data_nontest[,cols_to_numerize] = sapply(data_nontest[,cols_to_numerize],as.numeric)

# Recode FPQ choices with numerical values
data_nontest = data_nontest %>% 
  mutate(FPQ_1_num = FPQ_recoder(FPQ_1),
         FPQ_2_num = FPQ_recoder(FPQ_2),
         FPQ_3_num = FPQ_recoder(FPQ_3),
         FPQ_4_num = FPQ_recoder(FPQ_4),
         FPQ_5_num = FPQ_recoder(FPQ_5),
         FPQ_6_num = FPQ_recoder(FPQ_6),
         FPQ_7_num = FPQ_recoder(FPQ_7),
         FPQ_8_num = FPQ_recoder(FPQ_8),
         FPQ_9_num = FPQ_recoder(FPQ_9),
         FPQ_10_num = FPQ_recoder(FPQ_10),
         FPQ_11_num = FPQ_recoder(FPQ_11),
         FPQ_12_num = FPQ_recoder(FPQ_12),
         FPQ_13_num = FPQ_recoder(FPQ_13),
         FPQ_14_num = FPQ_recoder(FPQ_14),
         FPQ_15_num = FPQ_recoder(FPQ_15),
         FPQ_16_num = FPQ_recoder(FPQ_16),
         FPQ_17_num = FPQ_recoder(FPQ_17),
         FPQ_18_num = FPQ_recoder(FPQ_18),
         FPQ_19_num = FPQ_recoder(FPQ_19),
         FPQ_20_num = FPQ_recoder(FPQ_20),
         FPQ_21_num = FPQ_recoder(FPQ_21),
         FPQ_22_num = FPQ_recoder(FPQ_22),
         FPQ_23_num = FPQ_recoder(FPQ_23),
         FPQ_24_num = FPQ_recoder(FPQ_24),
         FPQ_25_num = FPQ_recoder(FPQ_25),
         FPQ_26_num = FPQ_recoder(FPQ_26),
         FPQ_27_num = FPQ_recoder(FPQ_27),
         FPQ_28_num = FPQ_recoder(FPQ_28),
         FPQ_29_num = FPQ_recoder(FPQ_29),
         FPQ_30_num = FPQ_recoder(FPQ_30))

# Compute FPQ_total score and pain tolerance change compared to baseline (positive values mean larger values in the experimental trial compared to baseline)
data_nontest = data_nontest %>% 
  rowwise() %>%
  mutate(FPQ_total = mean(c(FPQ_1_num, FPQ_2_num, FPQ_3_num, FPQ_4_num, FPQ_5_num, 
                            FPQ_6_num, FPQ_7_num, FPQ_8_num, FPQ_9_num, FPQ_10_num, 
                            FPQ_11_num, FPQ_12_num, FPQ_13_num, FPQ_14_num, FPQ_15_num, 
                            FPQ_16_num, FPQ_17_num, FPQ_18_num, FPQ_19_num, FPQ_20_num, 
                            FPQ_21_num, FPQ_22_num, FPQ_23_num, FPQ_24_num, FPQ_25_num, 
                            FPQ_26_num, FPQ_27_num, FPQ_28_num, FPQ_29_num, FPQ_30_num)),
         X1_paintolerance_change = X1_paintolerance - paintolerance_BL,
         X2_paintolerance_change = X2_paintolerance - paintolerance_BL,
         X1_painintensity_change = X1_painintensity - painintensity_BL,
         X2_painintensity_change = X2_painintensity - painintensity_BL,
         mean_water_temp = mean(c(watertemp_BL, X1_watertemp, X2_watertemp))
  ) %>% 
  ungroup()

# Recode age range, height range, which_real_hypnosis_ENG, and prefered_procedure_ENG
data_nontest = data_nontest %>%
  mutate(age = as.numeric(recode(age_range,
                                 "under 18" = "17", 
                                 "18 - 24" = "21",
                                 "25 - 34" = "30",
                                 "35 - 44" = "40",
                                 "45 - 54" = "50",
                                 "55 - 64" = "60",
                                 "65 - 74" = "70",
                                 "75 - 84" = "80",
                                 "85 or older" = "90",
                                 .default = NA_character_)),
         underaged = as.numeric(recode(age_range,
                                       "under 18" = "1",
                                       .default = "0")),
         height = as.numeric(recode(height_range,
                                    "smaller than 150 cm" = "145", 
                                    "150-155 cm" = "152",
                                    "156-160 cm" = "158",
                                    "161-165 cm" = "163",
                                    "166-170 cm" = "168",
                                    "171-175 cm" = "173",
                                    "176-180 cm" = "178",
                                    "181-185 cm" = "183",
                                    "186-190 cm" = "188",
                                    "191-195 cm" = "192",
                                    "196-200 cm" = "198",
                                    "201-205 cm" = "203",
                                    "205-210 cm" = "207",
                                    "higher than 210 cm" = "212",
                                    .default = NA_character_)),
         which_real_hypnosis_ENG = case_when((which_real_fieldname == "${e://Field/hypnosis_type_1}" & hypnosis_type_1 == "white noise hypnosis") ~ "whitenoise",
                                             (which_real_fieldname == "${e://Field/hypnosis_type_2}" & hypnosis_type_2 == "white noise hypnosis") ~ "whitenoise",
                                             (which_real_fieldname == "${e://Field/hypnosis_type_1}" & hypnosis_type_1 == "hypnotic relaxation") ~ "relaxation",
                                             (which_real_fieldname == "${e://Field/hypnosis_type_2}" & hypnosis_type_2 == "hypnotic relaxation") ~ "relaxation"),
        prefered_procedure_ENG = recode(prefered_procedure,
                                                          "fehér zaj hipnózis" = "white noise hypnosis",
                                                          "hipnotikus relaxáció" = "hypnotic relaxation")
         
  )

# Create deceived_by_sham variable
# only those get "Yes" on this variable who suspected deception AND were successful at guessing which condition was the real hypnosis

data_nontest = data_nontest %>%
  mutate(deceived_by_sham = case_when(suspected_control == "No" ~ "Yes",
                                      (suspected_control == "Yes" & which_real_hypnosis_ENG == "whitenoise") ~ "Yes",
                                      (suspected_control == "Yes" & which_real_hypnosis_ENG == "relaxation") ~ "No"))

### Long format

# convert data to long format to be able to use linear mixed model

data_long = data_nontest %>% 
  gather(key = "condition_paintolerance", value = "paintolerance_change", c(X1_paintolerance_change, X2_paintolerance_change)) %>%
  gather(key = "condition_painintensity", value = "painintensity_change", c(X1_painintensity_change, X2_painintensity_change)) %>%
  filter(substring(condition_paintolerance, 1, 2) == substring(condition_painintensity, 1, 2)) %>% 
  gather(key = "condition_exp_painred", value = "expectancy_pain_reduction", c(X1_exppainred_block1_1, X2_exppainred_block2_1)) %>%
  filter(substring(condition_paintolerance, 1, 2) == substring(condition_exp_painred, 1, 2)) %>% 
  gather(key = "condition_exp_hypdepth", value = "expectancy_hypnosis_depth", c(X1_exphypdepth, X2_exphypdepth)) %>%
  filter(substring(condition_paintolerance, 1, 2) == substring(condition_exp_hypdepth, 1, 2)) %>%
  gather(key = "condition_hypdepth", value = "hypnosis_depth", c(X1_hypdepth, X2_hypdepth)) %>%
  filter(substring(condition_paintolerance, 1, 2) == substring(condition_hypdepth, 1, 2)) %>%
  arrange(uniqueID)


# Creating variable block_number, which takes value 1 for the first experimental block and 2 for the second experimental bloxk
data_long = data_long %>%
  mutate(block_number = as.numeric(recode(condition_paintolerance,
                                          "X1_paintolerance_change" = "1",
                                          "X2_paintolerance_change" = "2")))

# Create condition variable
data_long = data_long %>%
  mutate(condition = case_when((group == "relaxation - whitenoise" & block_number == 1) ~ "relaxation",
                               (group == "whitenoise - relaxation" & block_number == 2) ~ "relaxation",
                               (group == "relaxation - whitenoise" & block_number == 2) ~ "whitenoise",
                               (group == "whitenoise - relaxation" & block_number == 1) ~ "whitenoise"))




#######################################
#         Dataset versioins           #
#######################################

# Only include participants eligible for inclusion in the study,
# that is, exclude people who tried hypnosis before or are under 18.
data_eligible_except_for_pain = data_nontest %>% 
  filter(triedhypnosisbefore == "No",
         underaged == 0)

nrow(data_eligible_except_for_pain)

# Only include participants eligible for inclusion in the study and pain analysis

# number excluded for high water temperature
sum((data_eligible_except_for_pain$watertemp_BL > 3) | (data_eligible_except_for_pain$X1_watertemp > 3) | (data_eligible_except_for_pain$X2_watertemp > 3))
#number excluded for not being deceived by sham
table(data_eligible_except_for_pain$deceived_by_sham)[1]

data_eligible_for_pain_intensity = data_eligible_except_for_pain %>% 
  filter(watertemp_BL < 3.01, ## water temperature not higher than 3 Celsius 
  X1_watertemp < 3.01, ## water temperature not higher than 3 Celsius
  X2_watertemp < 3.01, ## water temperature not higher than 3 Celsius
  deceived_by_sham == "Yes") ## who were deceived by the sham hypnosis

nrow(data_eligible_for_pain_intensity)

# Only include participants eligible for the confirmatory analysis



data_conf = data_eligible_for_pain_intensity %>% 
  filter(paintolerance_BL < 110.01) %>%  ## who have a pain tolerance not higher than 110 seconds at baseline
  filter((X1_paintolerance <110.01) | (X2_paintolerance <110.01)) ## the data of those will be excluded from the confirmatory analysis who have a pain tolerance higher than 110 seconds on both experimental blocks

nrow(data_conf)

sum(data_eligible_for_pain_intensity$paintolerance_BL < 110.01)
sum(data_eligible_for_pain_intensity$paintolerance_BL < 110.01)


### The same dataset versions for the long format data

# Only include participants eligible for inclusion in the study

data_eligible_except_for_pain_long = data_long %>% 
  filter(triedhypnosisbefore == "No",
         underaged == 0)

nrow(data_eligible_except_for_pain_long)

# Only include participants eligible for inclusion in the study and pain analysis

data_eligible_for_pain_intensity_long = data_eligible_except_for_pain_long %>% 
  filter(watertemp_BL < 3.01, ## water temperature not higher than 3 Celsius 
         X1_watertemp < 3.01, ## water temperature not higher than 3 Celsius
         X2_watertemp < 3.01, ## water temperature not higher than 3 Celsius
         deceived_by_sham == "Yes") ## who were deceived by the sham hypnosis

nrow(data_eligible_for_pain_intensity_long)

# Only include participants eligible for the confirmatory analysis

# confirmatory analysis eligibility:
## Only those participants’ data will be included in the confirmatory analysis, who 
## fit the general study eligibility criteria mentioned above, 
## who have a pain tolerance not higher than 110 seconds at baseline, 
## of whom the water temperature in the water tank during the experiment was not higher than 3 Celsius in any of the measurement points, 
## and who were deceived by the sham hypnosis. Being “deceived by sham” is determined by looking at two variables: 1) those we do not report having been suspicious that one of the procedures was not evidence based hypnosis induction will be treated as “deceived by sham”, also, 2) those who report having been suspicious that one of the procedures was not evidence based hypnosis induction, but did not correctly identify which is the real evidence based induction will also be treated as “deceived by sham”. 
## Furthermore, the data of those will be excluded from the confirmatory analysis who have a pain tolerance higher than 110 seconds on both experimental blocks (data of those who have a pain tolerance higher than 110 seconds on only one of the experimental blocks will be retained).


data_conf_long = data_eligible_for_pain_intensity_long %>% 
  filter(paintolerance_BL < 110.01) %>%  ## who have a pain tolerance not higher than 110 seconds at baseline
  filter(X1_paintolerance <110.01 | X2_paintolerance <110.01) ## the data of those will be excluded from the confirmatory analysis who have a pain tolerance higher than 110 seconds on both experimental blocks

nrow(data_conf_long)




### Datasets for sensitivity analysis, including participants who were not deceived by sham


data_eligible_for_pain_intensity_sensitivity = data_eligible_except_for_pain %>% 
  filter(watertemp_BL < 3.01, ## water temperature not higher than 3 Celsius 
         X1_watertemp < 3.01, ## water temperature not higher than 3 Celsius
         X2_watertemp < 3.01 ## water temperature not higher than 3 Celsius
         ) ## who were deceived by the sham hypnosis

nrow(data_eligible_for_pain_intensity_sensitivity)


data_conf_sensitivity = data_eligible_for_pain_intensity_sensitivity %>% 
  filter(paintolerance_BL < 110.01) %>%  ## who have a pain tolerance not higher than 110 seconds at baseline
  filter((X1_paintolerance <110.01) | (X2_paintolerance <110.01)) ## the data of those will be excluded from the confirmatory analysis who have a pain tolerance higher than 110 seconds on both experimental blocks

nrow(data_conf_sensitivity)



data_eligible_for_pain_intensity_long_sensitivity = data_eligible_except_for_pain_long %>% 
  filter(watertemp_BL < 3.01, ## water temperature not higher than 3 Celsius 
         X1_watertemp < 3.01, ## water temperature not higher than 3 Celsius
         X2_watertemp < 3.01 ## water temperature not higher than 3 Celsius
         ) ## who were deceived by the sham hypnosis

nrow(data_eligible_for_pain_intensity_long_sensitivity)

data_conf_long_sensitivity = data_eligible_for_pain_intensity_long_sensitivity %>% 
  filter(paintolerance_BL < 110.01) %>%  ## who have a pain tolerance not higher than 110 seconds at baseline
  filter(X1_paintolerance <110.01 | X2_paintolerance <110.01) ## the data of those will be excluded from the confirmatory analysis who have a pain tolerance higher than 110 seconds on both experimental blocks

nrow(data_conf_long_sensitivity)


#######################################
#            Data exploration         #
#######################################

# code for preliminary exploration of the data

### Deceived by sham table
table(data_eligible_except_for_pain$suspected_control, data_eligible_except_for_pain$which_real_hypnosis_ENG)

### Demographics

# eligible participants

data_eligible_except_for_pain%>%
  group_by(participant_gender)%>%
  summarize(n=n())%>%
  mutate(prop=paste0(round(n/sum(n),2)*100, "%"))

data_eligible_except_for_pain%>%
  group_by(age_range)%>%
  summarize(n=n())%>%
  mutate(prop=paste0(round(n/sum(n),2)*100, "%"))

data_eligible_except_for_pain%>%
  summarize(mean = round(mean(FPQ_total),2), sd = round(sd(FPQ_total),2))

data_eligible_except_for_pain%>%
  summarize(mean = round(mean(mean_water_temp),2), sd = round(sd(mean_water_temp),2))

# pain analysis sample

data_eligible_for_pain_intensity%>%
  group_by(participant_gender)%>%
  summarize(n=n())%>%
  mutate(prop=paste0(round(n/sum(n),2)*100, "%"))

data_eligible_for_pain_intensity%>%
  group_by(age_range)%>%
  summarize(n=n())%>%
  mutate(prop=paste0(round(n/sum(n),2)*100, "%"))

data_eligible_for_pain_intensity%>%
  summarize(mean = round(mean(FPQ_total),2), sd = round(sd(FPQ_total),2))

data_eligible_for_pain_intensity%>%
  summarize(mean = round(mean(mean_water_temp),2), sd = round(sd(mean_water_temp),2))


# confirmatory analysis sample

data_conf%>%
  group_by(participant_gender)%>%
  summarize(n=n())%>%
  mutate(prop=paste0(round(n/sum(n),2)*100, "%"))

data_conf%>%
  group_by(age_range)%>%
  summarize(n=n())%>%
  mutate(prop=paste0(round(n/sum(n),2)*100, "%"))

data_conf%>%
  summarize(mean = round(mean(FPQ_total),2), sd = round(sd(FPQ_total),2))

data_conf%>%
  summarize(mean = round(mean(mean_water_temp),2), sd = round(sd(mean_water_temp),2))


#######################################
#       Confirmatory analysis         #
#######################################

## H1: We hypothesize that sham hypnosis will evoke comparable expectations to conventional hypnosis in the following measures.

### H1.1 in expected pain reduction


data_conf_long %>% 
  group_by(condition) %>% 
  summarize(mean = mean(expectancy_pain_reduction),
            sd = sd(expectancy_pain_reduction))

data_conf_long %>% 
  ggplot() +
  aes(x = condition,y = expectancy_pain_reduction) +
  geom_violin(aes(fill = condition)) +
  geom_boxplot(width = 0.1) +
  geom_jitter(width = 0.1) +
  theme(text = element_text(size = 20))

SESOI_H1.1 = -10

mod_H1.1 = lmer(expectancy_pain_reduction ~ condition + (1|uniqueID), data = data_conf_long)
CI_H1.1 = confint(mod_H1.1, level = 0.9)["conditionwhitenoise",]
CI_lb_H1.1 = CI_H1.1[1]
CI_ub_H1.1 = CI_H1.1[2]

mod_H1.1
CI_H1.1

decision_H1.1 = if(SESOI_H1.1 < CI_lb_H1.1){"support_for_null"} else if (0 > CI_ub_H1.1){"support_for_alternative"} else {"inconclusive"} 
decision_H1.1

# Normality check

qqmath(mod_H1.1, id=0.05)

random_effects_mod_H1.1 = data.frame(ranef(mod_H1.1)[[1]])	
names(random_effects_mod_H1.1) = c("intercept")

qqmath(ranef(mod_H1.1))	

random_effects_mod_H1.1 %>% 	
  ggplot() +	
  aes(sample = intercept) +	
  stat_qq() +	
  stat_qq_line()	

# sensitivity analysis for H1.1

mod_H1.1_sensitivity = lmer(expectancy_pain_reduction ~ condition + (1|uniqueID), data = data_conf_long_sensitivity)
CI_H1.1_sensitivity = confint(mod_H1.1_sensitivity, level = 0.9)["conditionwhitenoise",]
CI_lb_H1.1_sensitivity = CI_H1.1_sensitivity[1]
CI_ub_H1.1_sensitivity = CI_H1.1_sensitivity[2]

mod_H1.1_sensitivity
CI_H1.1_sensitivity

decision_H1.1_sensitivity = if(SESOI_H1.1 < CI_lb_H1.1_sensitivity){"support_for_null"} else if (0 > CI_ub_H1.1_sensitivity){"support_for_alternative"} else {"inconclusive"} 
decision_H1.1_sensitivity


### H1.2 in expected hypnosis depth
data_conf_long %>% 
  group_by(condition) %>% 
  summarize(mean = mean(expectancy_hypnosis_depth),
            sd = sd(expectancy_hypnosis_depth))

data_conf_long %>% 
  ggplot() +
  aes(x = condition,y = expectancy_hypnosis_depth) +
  geom_violin(aes(fill = condition)) +
  geom_boxplot(width = 0.1) +
  geom_jitter(width = 0.1) +
  theme(text = element_text(size = 20))

SESOI_H1.2 = -1

mod_H1.2 = lmer(expectancy_hypnosis_depth ~ condition + (1|uniqueID), data = data_conf_long)
CI_H1.2 = confint(mod_H1.2, level = 0.9)["conditionwhitenoise",]
CI_lb_H1.2 = CI_H1.2[1]
CI_ub_H1.2 = CI_H1.2[2]

mod_H1.2
CI_H1.2

decision_H1.2 = if(SESOI_H1.2 < CI_lb_H1.2){"support_for_null"} else if (0 > CI_ub_H1.2){"support_for_alternative"} else {"inconclusive"} 
decision_H1.2

# Normality check

qqmath(mod_H1.2, id=0.05)

random_effects_mod_H1.2 = data.frame(ranef(mod_H1.2)[[1]])	
names(random_effects_mod_H1.2) = c("intercept")

qqmath(ranef(mod_H1.2))	

random_effects_mod_H1.2 %>% 	
  ggplot() +	
  aes(sample = intercept) +	
  stat_qq() +	
  stat_qq_line()	

# Sensitivity analysis for H1.2

mod_H1.2_sensitivity = lmer(expectancy_hypnosis_depth ~ condition + (1|uniqueID), data = data_conf_long_sensitivity)
CI_H1.2_sensitivity = confint(mod_H1.2_sensitivity, level = 0.9)["conditionwhitenoise",]
CI_lb_H1.2_sensitivity = CI_H1.2_sensitivity[1]
CI_ub_H1.2_sensitivity = CI_H1.2_sensitivity[2]

mod_H1.2_sensitivity
CI_H1.2_sensitivity

decision_H1.2_sensitivity = if(SESOI_H1.2 < CI_lb_H1.2_sensitivity){"support_for_null"} else if (0 > CI_ub_H1.2_sensitivity){"support_for_alternative"} else {"inconclusive"} 
decision_H1.2_sensitivity

## H2: We will test the hypothesis that conventional hypnosis will evoke greater hypnotic effects than sham hypnosis in the following measures.

### H2.1 Improved pain tolerance compared to baseline (improvement in the conventional hypnosis condition vs. the sham condition)

data_conf_long %>% 
  group_by(condition) %>% 
  summarize(mean(paintolerance_change),
            sd(paintolerance_change))

fig1 = data_conf_long %>% 
  ggplot() +
  aes(x = condition, y = paintolerance_change) +
  geom_violin(aes(fill = condition)) +
  geom_boxplot(width = 0.1) +
  geom_jitter(width = 0.1) +
  theme(text = element_text(size = 20)) +
  ylab("pain tolerance change (sec)") +
  scale_x_discrete(labels=c('relaxation', 'white noise')) +
  scale_fill_discrete(labels=c('relaxation', 'white noise'))
fig1

## save figure
# jpeg("fig1.jpg", units="in", width=8, height=6, res=300)
# fig1
# dev.off()


SESOI_H2.1 = -10

mod_H2.1 = lmer(paintolerance_change ~ age + participant_gender + height + FPQ_total + condition + mean_water_temp + (1|uniqueID), data = data_conf_long)
CI_H2.1 = confint(mod_H2.1, level = 0.9)["conditionwhitenoise",]
CI_lb_H2.1 = CI_H2.1[1]
CI_ub_H2.1 = CI_H2.1[2]

mod_H2.1
CI_H2.1

decision_H2.1 = if(SESOI_H2.1 < CI_lb_H2.1){"support_for_null"} else if (0 > CI_ub_H2.1){"support_for_alternative"} else {"inconclusive"} 
decision_H2.1


# Normality check

qqmath(mod_H2.1, id=0.05)

random_effects_mod_H2.1 = data.frame(ranef(mod_H2.1)[[1]])	
names(random_effects_mod_H2.1) = c("intercept")

qqmath(ranef(mod_H2.1))	

random_effects_mod_H2.1 %>% 	
  ggplot() +	
  aes(sample = intercept) +	
  stat_qq() +	
  stat_qq_line()	

### Normality assumption seems to be violated, running bootstrap procedure for CIs 

mySumm <- function(.) {
  s <- getME(., "sigma")
  c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta")))
}

boo1 <- bootstrap(model = mod_H2.1, .f = mySumm, type = "parametric", B = 20)

# bootstrapped confidence interval for condition
confint(boo1, level = 0.9)[6,]


# sensitivity analysis for H2.1

mod_H2.1_sensitivity = lmer(paintolerance_change ~ age + participant_gender + height + FPQ_total + condition + mean_water_temp + (1|uniqueID), data = data_conf_long_sensitivity)
CI_H2.1_sensitivity = confint(mod_H2.1_sensitivity, level = 0.9)["conditionwhitenoise",]
CI_lb_H2.1_sensitivity = CI_H2.1_sensitivity[1]
CI_ub_H2.1_sensitivity = CI_H2.1_sensitivity[2]

mod_H2.1_sensitivity
CI_H2.1_sensitivity

decision_H2.1_sensitivity = if(SESOI_H2.1 < CI_lb_H2.1_sensitivity){"support_for_null"} else if (0 > CI_ub_H2.1_sensitivity){"support_for_alternative"} else {"inconclusive"} 
decision_H2.1_sensitivity

### H2.2 in Hypnosis depth

data_conf_long %>% 
  group_by(condition) %>% 
  summarize(mean(hypnosis_depth),
            sd(hypnosis_depth))

fig2 = data_conf_long %>% 
  ggplot() +
  aes(x = condition, y = hypnosis_depth) +
  geom_violin(aes(fill = condition)) +
  geom_boxplot(width = 0.1) +
  geom_jitter(width = 0.1) +
  theme(text = element_text(size = 20)) +
  ylab("self reported hypnosis depth (0-10)") +
  scale_x_discrete(labels=c('relaxation', 'white noise')) +
  scale_fill_discrete(labels=c('relaxation', 'white noise'))
fig2

## save figure
# jpeg("fig2.jpg", units="in", width=8, height=6, res=300)
# fig2
# dev.off()

SESOI_H2.2 = -1

mod_H2.2 = lmer(hypnosis_depth ~ condition + (1|uniqueID), data = data_conf_long)
CI_H2.2 = confint(mod_H2.2, level = 0.9)["conditionwhitenoise",]
CI_lb_H2.2 = CI_H2.2[1]
CI_ub_H2.2 = CI_H2.2[2]

mod_H2.2
CI_H2.2

decision_H2.2 = if(SESOI_H2.2 < CI_lb_H2.2){"support_for_null"} else if (0 > CI_ub_H2.2){"support_for_alternative"} else {"inconclusive"} 
decision_H2.2


# Normality check

qqmath(mod_H2.2, id=0.05)

random_effects_mod_H2.2 = data.frame(ranef(mod_H2.2)[[1]])	
names(random_effects_mod_H2.2) = c("intercept")

qqmath(ranef(mod_H2.2))	

random_effects_mod_H2.2 %>% 	
  ggplot() +	
  aes(sample = intercept) +	
  stat_qq() +	
  stat_qq_line()	

# sensitivity analysis for H2.1

mod_H2.2_sensitivity = lmer(hypnosis_depth ~ condition + (1|uniqueID), data = data_conf_long_sensitivity)
CI_H2.2_sensitivity = confint(mod_H2.2_sensitivity, level = 0.9)["conditionwhitenoise",]
CI_lb_H2.2_sensitivity = CI_H2.2_sensitivity[1]
CI_ub_H2.2_sensitivity = CI_H2.2_sensitivity[2]

mod_H2.2_sensitivity
CI_H2.2_sensitivity

decision_H2.2_sensitivity = if(SESOI_H2.2 < CI_lb_H2.2_sensitivity){"support_for_null"} else if (0 > CI_ub_H2.2_sensitivity){"support_for_alternative"} else {"inconclusive"} 
decision_H2.2_sensitivity


#######################################
#       Exploratory analysis          #
#######################################

### How convincing is white noise hypnosis?

table(data_eligible_except_for_pain$suspected_control)
20/97

data_eligible_except_for_pain %>% 
  filter(suspected_control == "Yes") %>% 
  group_by(which_real_hypnosis_ENG) %>% 
  summarize(n = n())
7/20

table(data_eligible_except_for_pain$which_real_hypnosis_ENG)
binom.test(x = 62, n = 97, p = 0.5)
62/97

table(data_eligible_except_for_pain$prefered_procedure_ENG)
52/97


### H2 just for pain intensity

## on the pain analysis subsample

mod_H2.1a = lmer(painintensity_change ~ age + participant_gender + height + FPQ_total + condition + mean_water_temp + (1|uniqueID), data = data_eligible_for_pain_intensity_long)
CI_H2.1a = confint(mod_H2.1a, level = 0.9)["conditionwhitenoise",]
CI_lb_H2.1a = CI_H2.1a[1]
CI_ub_H2.1a = CI_H2.1a[2]

CI_H2.1a

data_eligible_for_pain_intensity_long %>%
  summarize(mean(painintensity_change),
            sd(painintensity_change))

data_eligible_for_pain_intensity_long %>%
  group_by(condition) %>% 
  summarize(mean(painintensity_change),
            sd(painintensity_change))

data_eligible_for_pain_intensity_long %>% 
  ggplot() +
  aes(x = condition, y = painintensity_change) +
  geom_violin(aes(fill = condition)) +
  geom_boxplot(width = 0.1) +
  geom_jitter(width = 0.1) +
  theme(text = element_text(size = 20))


## on the confirmatory analysis subsample

mod_H2.1a = lmer(painintensity_change ~ age + participant_gender + height + FPQ_total + condition + mean_water_temp + (1|uniqueID), data = data_conf_long)
CI_H2.1a = confint(mod_H2.1a, level = 0.9)["conditionwhitenoise",]
CI_lb_H2.1a = CI_H2.1a[1]
CI_ub_H2.1a = CI_H2.1a[2]

CI_H2.1a

data_conf_long %>%
  summarize(mean(painintensity_change),
            sd(painintensity_change))

data_conf_long %>%
  group_by(condition) %>% 
  summarize(mean(painintensity_change),
            sd(painintensity_change))

data_conf_long %>% 
  ggplot() +
  aes(x = condition, y = painintensity_change) +
  geom_violin(aes(fill = condition)) +
  geom_boxplot(width = 0.1) +
  geom_jitter(width = 0.1) +
  theme(text = element_text(size = 20))



### effect of hypnotizability on pain reduction

data_conf_long_relaxation_only = data_conf_long %>%
  filter(condition == "relaxation")

data_conf_long_whitenoise_only = data_conf_long %>%
  filter(condition == "whitenoise")

data_eligible_for_pain_intensity_long_relaxation_only = data_eligible_for_pain_intensity_long %>%
  filter(condition == "relaxation")

data_eligible_for_pain_intensity_long_whitenoise_only = data_eligible_for_pain_intensity_long %>%
  filter(condition == "whitenoise")


IV_list = c("hypnotizability_total", "hypnosis_depth", "expectancy_pain_reduction", "expectancy_hypnosis_depth")
DV_list = c("paintolerance_change", "painintensity_change", "hypnosis_depth")

IV_CV_combinations = expand.grid(IV_list, DV_list)
names(IV_CV_combinations) = c("IVs", "DVs")

output_table_conf_long_relaxation_only = data.frame(paintolerance_change = rep(NA, 4), painintensity_change = rep(NA, 4), hypnosis_depth = rep(NA, 4))
row.names(output_table_conf_long_relaxation_only) = IV_list

for(i in 1:nrow(IV_CV_combinations)){
  output_table_conf_long_relaxation_only[as.character(IV_CV_combinations[i,"IVs"]), as.character(IV_CV_combinations[i,"DVs"])] = cor(data_conf_long_relaxation_only[,IV_CV_combinations[i,"IVs"]], data_conf_long_relaxation_only[,IV_CV_combinations[i,"DVs"]], use="complete.obs")
}


output_table_conf_long_whitenoise_only = data.frame(paintolerance_change = rep(NA, 4), painintensity_change = rep(NA, 4), hypnosis_depth = rep(NA, 4))
row.names(output_table_conf_long_whitenoise_only) = IV_list

for(i in 1:nrow(IV_CV_combinations)){
  output_table_conf_long_whitenoise_only[as.character(IV_CV_combinations[i,"IVs"]), as.character(IV_CV_combinations[i,"DVs"])] = cor(data_conf_long_whitenoise_only[,IV_CV_combinations[i,"IVs"]], data_conf_long_whitenoise_only[,IV_CV_combinations[i,"DVs"]], use="complete.obs")
}



output_table_eligible_for_pain_intensity_long_relaxation_only = data.frame(paintolerance_change = rep(NA, 4), painintensity_change = rep(NA, 4), hypnosis_depth = rep(NA, 4))
row.names(output_table_eligible_for_pain_intensity_long_relaxation_only) = IV_list

for(i in 1:nrow(IV_CV_combinations)){
  output_table_eligible_for_pain_intensity_long_relaxation_only[as.character(IV_CV_combinations[i,"IVs"]), as.character(IV_CV_combinations[i,"DVs"])] = cor(data_eligible_for_pain_intensity_long_relaxation_only[,IV_CV_combinations[i,"IVs"]], data_eligible_for_pain_intensity_long_relaxation_only[,IV_CV_combinations[i,"DVs"]], use="complete.obs")
}


output_table_eligible_for_pain_intensity_long_whitenoise_only = data.frame(paintolerance_change = rep(NA, 4), painintensity_change = rep(NA, 4), hypnosis_depth = rep(NA, 4))
row.names(output_table_eligible_for_pain_intensity_long_whitenoise_only) = IV_list

for(i in 1:nrow(IV_CV_combinations)){
  output_table_eligible_for_pain_intensity_long_whitenoise_only[as.character(IV_CV_combinations[i,"IVs"]), as.character(IV_CV_combinations[i,"DVs"])] = cor(data_eligible_for_pain_intensity_long_whitenoise_only[,IV_CV_combinations[i,"IVs"]], data_eligible_for_pain_intensity_long_whitenoise_only[,IV_CV_combinations[i,"DVs"]], use="complete.obs")
}

output_table_conf_long_relaxation_only
output_table_conf_long_whitenoise_only
output_table_eligible_for_pain_intensity_long_relaxation_only
output_table_eligible_for_pain_intensity_long_whitenoise_only

### Some visual analyses of the correlations

data_conf_long %>% 
  ggplot() +
  aes(x = hypnotizability_total, y = paintolerance_change, color = condition) +
  geom_jitter(aes(color = condition)) +
  geom_smooth(method = "lm") +
  theme(text = element_text(size = 20))


data_eligible_for_pain_intensity_long %>%  
  ggplot() +
  aes(x = hypnotizability_total, y = painintensity_change, color = condition) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme(text = element_text(size = 20))


data_eligible_for_pain_intensity_long %>%  
  ggplot() +
  aes(x = hypnotizability_total, y = hypnosis_depth, color = condition) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme(text = element_text(size = 20))

data_conf_long %>% 
  ggplot() +
  aes(x = hypnosis_depth, y = paintolerance_change, color = condition) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme(text = element_text(size = 20))

data_eligible_for_pain_intensity_long %>%  
  ggplot() +
  aes(x = hypnosis_depth, y = painintensity_change, color = condition) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme(text = element_text(size = 20))

data_conf_long %>% 
  ggplot() +
  aes(x = expectancy_pain_reduction, y = paintolerance_change, color = condition) +
  geom_jitter(aes(color = condition)) +
  geom_smooth(method = "lm") +
  theme(text = element_text(size = 20))

data_eligible_for_pain_intensity_long %>%  
  ggplot() +
  aes(x = expectancy_pain_reduction, y = painintensity_change, color = condition) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme(text = element_text(size = 20))

data_conf_long %>% 
  ggplot() +
  aes(x = expectancy_hypnosis_depth, y = paintolerance_change, color = condition) +
  geom_jitter(aes(color = condition)) +
  geom_smooth(method = "lm") +
  theme(text = element_text(size = 20))

data_eligible_for_pain_intensity_long %>%  
  ggplot() +
  aes(x = expectancy_hypnosis_depth, y = painintensity_change, color = condition) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme(text = element_text(size = 20))


#######################################
#             Side effects            #
#######################################


table(data_nontest$felt_negativeeffect)

na.omit(data_nontest$what_negativeeffect)

data_nontest[data_nontest$felt_negativeeffect == "Yes","what_negativeeffect"]

table(data_nontest$prefered_procedure)


data_nontest_recoded = data_nontest %>% 
  mutate(prefered_procedure_onelanguage = recode(prefered_procedure,
                                                 "fehér zaj hipnózis" = "white noise hypnosis",
                                                 "hipnotikus relaxáció" = "hypnotic relaxation"))

table(data_nontest_recoded$prefered_procedure_onelanguage)


chisq.test(table(data_nontest_recoded$prefered_procedure_onelanguage))











## recode experience with hypnosis lecture

# data_nontest$lecture = 0
# data_nontest[grep("Attended a lecture about hypnosis", data_nontest$knowledge_source), "lecture"] = 1




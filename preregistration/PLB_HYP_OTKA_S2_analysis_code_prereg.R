

#######################################
#  PLB-HYP OTKA S2 analysis code plan #
#######################################

#######################################
#              Packages               #
#######################################

library(lme4) # for lmer()
library(tidyverse) # for tidy code

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

# code that loads and cleans data will be added here
# in this code we use a dataset from test sessions of the study
# which is already cleaned

data_pre = read.csv("https://raw.githubusercontent.com/kekecsz/PLB-HYP-OTKA-S2/main/preregistration/OTKA%20-%20Study2_clean_test_data_for_script_testing.csv")

data_pre = data_pre %>% 
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
                            
data_pre = data_pre %>% 
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
         )

data_pre = data_pre %>%
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
                                             (which_real_fieldname == "${e://Field/hypnosis_type_2}" & hypnosis_type_2 == "hypnotic relaxation") ~ "relaxation")
         )






data_pre = data_pre %>%
  mutate(deceived_by_sham = case_when(suspected_control == "No" ~ "Yes",
                                      (suspected_control == "Yes" & which_real_hypnosis_ENG == "whitenoise") ~ "Yes",
                                      (suspected_control == "Yes" & which_real_hypnosis_ENG == "relaxation") ~ "No"))



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

# participants attest that none of these criteria are true for them in the informed consent form
# however, there are some questions in the form that can detect inconsistencies in this report,
# so we exclude those who are ineligible based on these additional questions

data_pre = data_pre %>%
  filter(age > 17) %>% # Are at least 18 years old 
  filter(triedhypnosisbefore == "No") # Have not tried hypnosis before



#######################################
#            Data exploration         #
#######################################

# code for preliminary exploration of the data will be added here
# to get a sense of the distribution of the variables and to
# detect potential data anomalies or errors.

summary(data_pre)



#######################################
#       Confirmatory analysis         #
#######################################

# exclude participants who are not eligible to be entered into the confirmatory analysis

# confirmatory analysis eligibility:

## Only those participants’ data will be included in the confirmatory analysis, who 
## fit the general study eligibility criteria mentioned above, 
## who have a pain tolerance not higher than 110 seconds at baseline, 
## of whom the water temperature in the water tank during the experiment was not higher than 3 Celsius in any of the measurement points, 
## and who were deceived by the sham hypnosis. Being “deceived by sham” is determined by looking at two variables: 1) those we do not report having been suspicious that one of the procedures was not evidence based hypnosis induction will be treated as “deceived by sham”, also, 2) those who report having been suspicious that one of the procedures was not evidence based hypnosis induction, but did not correctly identify which is the real evidence based induction will also be treated as “deceived by sham”. 
## Furthermore, the data of those will be excluded from the confirmatory analysis who have a pain tolerance higher than 110 seconds on both experimental blocks (data of those who have a pain tolerance higher than 110 seconds on only one of the experimental blocks will be retained).

data_conf = data_pre %>% 
  filter(paintolerance_BL < 110.01) %>% ## who have a pain tolerance not higher than 110 seconds at baseline
  filter(watertemp_BL < 3.01) %>% ## water temperature not higher than 3 Celsius 
  filter(X1_watertemp < 3.01) %>% ## water temperature not higher than 3 Celsius
  filter(X2_watertemp < 3.01) %>% ## water temperature not higher than 3 Celsius
  filter(deceived_by_sham == "Yes") %>% ## who were deceived by the sham hypnosis
  filter(X1_paintolerance <110.01 | X2_paintolerance <110.01) ## the data of those will be excluded from the confirmatory analysis who have a pain tolerance higher than 110 seconds on both experimental blocks


# convert data to long format to be able to use linear mixed model

data_long = data_conf %>% 
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



data_long = data_long %>%
  mutate(block_number = as.numeric(recode(condition_paintolerance,
                               "X1_paintolerance_change" = "1",
                               "X2_paintolerance_change" = "2")))

data_long = data_long %>%
  mutate(condition = case_when((group == "relaxation - whitenoise" & block_number == 1) ~ "relaxation",
                               (group == "whitenoise - relaxation" & block_number == 2) ~ "relaxation",
                               (group == "relaxation - whitenoise" & block_number == 2) ~ "whitenoise",
                               (group == "whitenoise - relaxation" & block_number == 1) ~ "whitenoise"))


# check if recoding was correct
cbind(data_long$hypnosis_type_1, data_long$hypnosis_type_2, data_long$condition)



## H1: We hypothesize that sham hypnosis will evoke comparable expectations to conventional hypnosis in the following measures.

### H1.1 in expected pain reduction

SESOI_H1.1 = -10

mod_H1.1 = lmer(expectancy_pain_reduction ~ condition + (1|uniqueID), data = data_long)
CI_H1.1 = confint(mod_H1.1, level = 0.9)["conditionwhitenoise",]
CI_lb_H1.1 = CI_H1.1[1]
CI_ub_H1.1 = CI_H1.1[2]

CI_H1.1

decision_H1.1 = if(SESOI_H1.1 < CI_lb_H1.1){"support_for_null"} else if (0 > CI_ub_H1.1){"support_for_alternative"} else {"inconclusive"} 
decision_H1.1


### H1.2 in expected hypnosis depth

SESOI_H1.2 = -1

mod_H1.2 = lmer(expectancy_hypnosis_depth ~ condition + (1|uniqueID), data = data_long)
CI_H1.2 = confint(mod_H1.2, level = 0.9)["conditionwhitenoise",]
CI_lb_H1.2 = CI_H1.2[1]
CI_ub_H1.2 = CI_H1.2[2]

CI_H1.2

decision_H1.2 = if(SESOI_H1.2 < CI_lb_H1.2){"support_for_null"} else if (0 > CI_ub_H1.2){"support_for_alternative"} else {"inconclusive"} 
decision_H1.2


## H2: We will test the hypothesis that conventional hypnosis will evoke greater hypnotic effects than sham hypnosis in the following measures.

### H2.1 Improved pain tolerance compared to baseline (improvement in the conventional hypnosis condition vs. the sham condition)

SESOI_H2.1 = -10

mod_H2.1 = lmer(paintolerance_change ~ age + participant_gender + height + FPQ_total + condition + mean_water_temp + (1|uniqueID), data = data_long)
CI_H2.1 = confint(mod_H2.1, level = 0.9)["conditionwhitenoise",]
CI_lb_H2.1 = CI_H2.1[1]
CI_ub_H2.1 = CI_H2.1[2]

CI_H2.1

decision_H2.1 = if(SESOI_H2.1 < CI_lb_H2.1){"support_for_null"} else if (0 > CI_ub_H2.1){"support_for_alternative"} else {"inconclusive"} 
decision_H2.1


### H2.2 in Hypnosis depth

SESOI_H2.2 = -1

mod_H2.2 = lmer(hypnosis_depth ~ condition + (1|uniqueID), data = data_long)
CI_H2.2 = confint(mod_H2.2, level = 0.9)["conditionwhitenoise",]
CI_lb_H2.2 = CI_H2.2[1]
CI_ub_H2.2 = CI_H2.2[2]

CI_H2.2

decision_H2.2 = if(SESOI_H2.2 < CI_lb_H2.2){"support_for_null"} else if (0 > CI_ub_H2.2){"support_for_alternative"} else {"inconclusive"} 
decision_H2.2


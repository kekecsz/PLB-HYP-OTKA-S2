
###################################################
#  PLB-HYP OTKA S2 sample size planning code plan #
###################################################

#######################################
#              Packages               #
#######################################

library(MASS)
library(tidyverse)
library(lme4)
library(pbapply)
library(BayesFactor)

#######################################
#         Custom functions            #
#######################################

# this is the data generator function

# This function contains assumptions of correlations of the variables that are being simulated
# in the mvrnorm() function. This assumed correlation matrix has been built based on prior research in the field
# our previous studies, the 8 test sessions run before the live experiment (5 of whom would be eligible for 
# confirmatory analysis), and the first 12 live session participants 
# (only 7 of whom would be eligible for confirmatory analysis) without looking at their condition, 
# to keep the preregistration blind to the results.
# The correlation matrix is found here, for easier access: https://github.com/kekecsz/PLB-HYP-OTKA-S2/blob/main/preregistration/PLB_HYP_OTAK_S2_power%20analysis%20correlation%20matrix.xlsx 

simulate_data_with_moderators = function(
  n,
  mean_paintolerance_hyp,
  mean_paintolerance_plb,
  mean_paintolerance_baseline,
  sd_paintolerance,
  sd_age,
  mean_age,
  sd_height,
  mean_height,
  sd_expected_CPT_pain,
  mean_expected_CPT_pain,
  sd_FPQ_total,
  mean_FPQ_total,
  sd_expectancy_pain_reduction,
  mean_expectancy_pain_reduction_hyp,
  mean_expectancy_pain_reduction_plb,
  sd_expectancy_hypnotic_depth,
  mean_expectancy_hypnotic_depth_hyp,
  mean_expectancy_hypnotic_depth_plb,
  sd_hypnotic_depth,
  mean_hypnotic_depth_hyp,
  mean_hypnotic_depth_plb,
  sd_hypnotizability,
  mean_hypnotizability,
  sd_mean_water_temp,
  mean_mean_water_temp
){
  varnames_sim_pre = c(
    "paintolerance_baseline", "paintolerance_hyp", "paintolerance_plb", "age", "female", "height", "expected_CPT_pain", "FPQ_total", "expectancy_pain_reduction_hyp", "expectancy_pain_reduction_plb", "expectancy_hypnotic_depth_hyp", "expectancy_hypnotic_depth_plb", "hypnotic_depth_hyp", "hypnotic_depth_plb", "hypnotizability", "mean_water_temp"
  )
  
  data_sim_pre = as.data.frame(
    mvrnorm(n = n,
            mu = rep(0, length(varnames_sim_pre)),
            Sigma = matrix(c(  1,0.5,0.5,0.2,0.2,0.2,-0.2,-0.2,0,0,0,0,0,0,0,0.25,
                               0.5,1,0.7,0.2,0.2,0.2,-0.2,-0.2,0.5,0.3,0.4,0.2,0.5,0.3,0.4,0.25,
                               0.5,0.7,1,0.2,0.2,0.2,-0.2,-0.2,0.3,0.5,0.2,0.4,0.3,0.5,0.4,0.25,
                               0.2,0.2,0.2,1,0,0,0,0,0,0,0,0,0,0,0,0,
                               0.2,0.2,0.2,0,1,-0.2,0,0,0,0,0,0,0,0,0,0,
                               0.2,0.2,0.2,0,-0.2,1,0,0,0,0,0,0,0,0,0,0,
                               -0.2,-0.2,-0.2,0,0,0,1,0,0,0,0,0,0,0,0,0,
                               -0.2,-0.2,-0.2,0,0,0,0,1,0,0,0,0,0,0,0,0,
                               0,0.5,0.3,0,0,0,0,0,1,0.6,0.5,0.4,0.3,0.3,0.4,0,
                               0,0.3,0.5,0,0,0,0,0,0.6,1,0.4,0.5,0.3,0.3,0.4,0,
                               0,0.4,0.2,0,0,0,0,0,0.5,0.4,1,0.6,0.4,0.3,0.4,0,
                               0,0.2,0.4,0,0,0,0,0,0.4,0.5,0.6,1,0.3,0.4,0.4,0,
                               0,0.5,0.3,0,0,0,0,0,0.3,0.3,0.4,0.3,1,0.6,0.6,0,
                               0,0.3,0.5,0,0,0,0,0,0.3,0.3,0.3,0.4,0.6,1,0.6,0,
                               0,0.4,0.4,0,0,0,0,0,0.4,0.4,0.4,0.4,0.6,0.6,1,0,
                               0.25,0.25,0.25,0,0,0,0,0,0,0,0,0,0,0,0,1
                               
                             
            ), nrow = length(varnames_sim_pre)))
  )
 
  
  
  names(data_sim_pre) = varnames_sim_pre
  
  data_pre = as.data.frame(matrix(NA, nrow = n, ncol = length(varnames_sim_pre)+1))
  names(data_pre) = c("ID", varnames_sim_pre)
  
  data_pre[,"paintolerance_baseline"] = data_sim_pre[,"paintolerance_baseline"]*sd_paintolerance+mean_paintolerance_baseline
  data_pre[,"paintolerance_hyp"] = data_sim_pre[,"paintolerance_hyp"]*sd_paintolerance+mean_paintolerance_hyp
  data_pre[,"paintolerance_plb"] = data_sim_pre[,"paintolerance_plb"]*sd_paintolerance+mean_paintolerance_plb
  data_pre[,"age"] = data_sim_pre[,"age"]*sd_age+mean_age
  data_pre[data_sim_pre[,"female"] < 0,"female"] = 0
  data_pre[data_sim_pre[,"female"] > -0.00000000001,"female"] = 1
  data_pre[,"height"] = data_sim_pre[,"height"]*sd_height+mean_height
  data_pre[,"expected_CPT_pain"] = data_sim_pre[,"expected_CPT_pain"]*sd_expected_CPT_pain+mean_expected_CPT_pain
  data_pre[,"FPQ_total"] = data_sim_pre[,"FPQ_total"]*sd_FPQ_total+mean_FPQ_total
  data_pre[,"ID"] = 1:n
  data_pre[,"expectancy_pain_reduction_hyp"] = data_sim_pre[,"expectancy_pain_reduction_hyp"]*sd_expectancy_pain_reduction+mean_expectancy_pain_reduction_hyp
  data_pre[,"expectancy_pain_reduction_plb"] = data_sim_pre[,"expectancy_pain_reduction_plb"]*sd_expectancy_pain_reduction+mean_expectancy_pain_reduction_plb
  data_pre[,"expectancy_hypnotic_depth_hyp"] = data_sim_pre[,"expectancy_hypnotic_depth_hyp"]*sd_expectancy_hypnotic_depth+mean_expectancy_hypnotic_depth_hyp
  data_pre[,"expectancy_hypnotic_depth_plb"] = data_sim_pre[,"expectancy_hypnotic_depth_plb"]*sd_expectancy_hypnotic_depth+mean_expectancy_hypnotic_depth_plb
  data_pre[,"hypnotic_depth_hyp"] = data_sim_pre[,"hypnotic_depth_hyp"]*sd_hypnotic_depth+mean_hypnotic_depth_hyp
  data_pre[,"hypnotic_depth_plb"] = data_sim_pre[,"hypnotic_depth_plb"]*sd_hypnotic_depth+mean_hypnotic_depth_plb
  data_pre[,"hypnotizability"] = data_sim_pre[,"hypnotizability"]*sd_hypnotizability+mean_hypnotizability
  data_pre[,"mean_water_temp"] = data_sim_pre[,"mean_water_temp"]*sd_mean_water_temp+mean_mean_water_temp

  data_long = data_pre %>% 
    gather(key = "condition", value = "paintolerance", paintolerance_hyp:paintolerance_plb) %>%
    gather(key = "condition_exp_painred", value = "expectancy_pain_reduction", expectancy_pain_reduction_hyp:expectancy_pain_reduction_plb) %>%
    filter(substring(condition, 15, 17) == substring(condition_exp_painred, 27, 29)) %>% 
    gather(key = "condition_exp_hypdepth", value = "expectancy_hypnotic_depth", expectancy_hypnotic_depth_hyp:expectancy_hypnotic_depth_plb) %>%
    filter(substring(condition, 15, 17) == substring(condition_exp_hypdepth, 27, 29)) %>%
    gather(key = "condition_hypdepth", value = "hypnotic_depth", hypnotic_depth_hyp:hypnotic_depth_plb) %>%
    filter(substring(condition, 15, 17) == substring(condition_hypdepth, 16, 18)) %>%
    arrange(ID)
  
  data_long = data_long[, -which(names(data_long) %in% c("condition_exp_painred", "condition_exp_hypdepth", "condition_hypdepth"))]
  
  data_long$paintolerance_change = data_long$paintolerance - data_long$paintolerance_baseline
  
  data_long = data_long %>% 
    mutate(condition = as.factor(recode(condition,
                                        paintolerance_hyp = "hyp",
                                        paintolerance_plb = "plb")),
           ID = as.factor(ID))
  return(data_long)
}



# this is the data analysis function (calling on the data generator function to generate new data each time)

simul = function(
  n,
  mean_paintolerance_hyp,
  mean_paintolerance_plb,
  mean_paintolerance_baseline,
  sd_paintolerance,
  sd_age,
  mean_age,
  sd_height,
  mean_height,
  sd_expected_CPT_pain,
  mean_expected_CPT_pain,
  sd_FPQ_total,
  mean_FPQ_total,
  sd_expectancy_pain_reduction,
  mean_expectancy_pain_reduction_hyp,
  mean_expectancy_pain_reduction_plb,
  sd_expectancy_hypnotic_depth,
  mean_expectancy_hypnotic_depth_hyp,
  mean_expectancy_hypnotic_depth_plb,
  sd_hypnotic_depth,
  mean_hypnotic_depth_hyp,
  mean_hypnotic_depth_plb,
  sd_hypnotizability,
  mean_hypnotizability,
  sd_mean_water_temp,
  mean_mean_water_temp,
  which_hypothesis_to_test,
  SESOI
){

  data = simulate_data_with_moderators(
    n = n,
    mean_paintolerance_hyp = mean_paintolerance_hyp,
    mean_paintolerance_plb = mean_paintolerance_plb,
    mean_paintolerance_baseline = mean_paintolerance_baseline,
    sd_paintolerance = sd_paintolerance,
    sd_age = sd_age,
    mean_age = mean_age,
    sd_height = sd_height,
    mean_height = mean_height,
    sd_expected_CPT_pain = sd_expected_CPT_pain,
    mean_expected_CPT_pain = mean_expected_CPT_pain,
    sd_FPQ_total = sd_FPQ_total,
    mean_FPQ_total = mean_FPQ_total,
    sd_expectancy_pain_reduction = sd_expectancy_pain_reduction,
    mean_expectancy_pain_reduction_hyp = mean_expectancy_pain_reduction_hyp,
    mean_expectancy_pain_reduction_plb = mean_expectancy_pain_reduction_plb,
    sd_expectancy_hypnotic_depth = sd_expectancy_hypnotic_depth,
    mean_expectancy_hypnotic_depth_hyp = mean_expectancy_hypnotic_depth_hyp,
    mean_expectancy_hypnotic_depth_plb = mean_expectancy_hypnotic_depth_plb,
    sd_hypnotic_depth = sd_hypnotic_depth,
    mean_hypnotic_depth_hyp = mean_hypnotic_depth_hyp,
    mean_hypnotic_depth_plb = mean_hypnotic_depth_plb,
    sd_hypnotizability = sd_hypnotizability,
    mean_hypnotizability = mean_hypnotizability,
    sd_mean_water_temp = sd_mean_water_temp,
    mean_mean_water_temp = sd_mean_water_temp
  )
  
  
  if(which_hypothesis_to_test == "H1.1"){
    mod = lmer(expectancy_pain_reduction ~ condition + (1|ID), data = data)
    CI = confint(mod, level = 0.9)["conditionplb",]
    CI_lb = CI[1]
    CI_ub = CI[2]
    
    decision = if(SESOI < CI_lb){"H0"} else if (0 > CI_ub){"H1"} else {"inconclusive"} 
  
    } else if(which_hypothesis_to_test == "H1.2"){
    mod = lmer(expectancy_hypnotic_depth ~ condition + (1|ID), data = data)
    CI = confint(mod, level = 0.9)["conditionplb",]
    CI_lb = CI[1]
    CI_ub = CI[2]
    
    decision = if(SESOI < CI_lb){"H0"} else if (0 > CI_ub){"H1"} else {"inconclusive"} 
    
    } else if(which_hypothesis_to_test == "H2.1"){
      mod = lmer(paintolerance_change ~ age + female + height + FPQ_total + condition + mean_water_temp + (1|ID), data = data)
      CI = confint(mod, level = 0.9)["conditionplb",]
      CI_lb = CI[1]
      CI_ub = CI[2]
      
      decision = if(SESOI < CI_lb){"H0"} else if (0 > CI_ub){"H1"} else {"inconclusive"} 
      
    } else if(which_hypothesis_to_test == "H2.2"){
      mod = lmer(hypnotic_depth ~ condition + (1|ID), data = data)
      CI = confint(mod, level = 0.9)["conditionplb",]
      CI_lb = CI[1]
      CI_ub = CI[2]
      
      decision = if(SESOI < CI_lb){"H0"} else if (0 > CI_ub){"H1"} else {"inconclusive"} 
      
    }

  return(decision)
  
}


#######################################
#             Simulation              #
#######################################


# set number of simulated samples (iterations) here

# Important: This analysis runs for a long time due to
# the mixed models' processing power needs. On an Intel i7 processor
# 10000 iterations run for about 16 hours.
# to decrease run time you can decrease the number of iterations,
# 1000 iterations will provide a decent estimate and it runs for only
# 1.5 hours.

iter=10000


# in the simul() function we set the appropriate parameters to either simulate the alternative hypothesis being true
# or the null hypothesis being true for each hypothesis


## H1: We hypothesize that sham hypnosis will evoke comparable expectations to conventional hypnosis in the following measures.
### H1.1 in expected pain reduction

### H1.1 - alternative true


results_H1.1_alttrue = pbreplicate(iter, simul(n = 45,
                                               mean_paintolerance_hyp = 75,
                                               mean_paintolerance_plb = 75,
                                               mean_paintolerance_baseline = 50,
                                               sd_paintolerance = 32,
                                               sd_age = 2.5,
                                               mean_age = 22,
                                               sd_height = 14,
                                               mean_height = 175,
                                               sd_expected_CPT_pain = 2,
                                               mean_expected_CPT_pain = 5.5,
                                               sd_FPQ_total = 0.5,
                                               mean_FPQ_total = 3,
                                               sd_expectancy_pain_reduction = 20,
                                               mean_expectancy_pain_reduction_hyp = 50,
                                               mean_expectancy_pain_reduction_plb = 40,
                                               sd_expectancy_hypnotic_depth = 2,
                                               mean_expectancy_hypnotic_depth_hyp = 5,
                                               mean_expectancy_hypnotic_depth_plb = 5,
                                               sd_hypnotic_depth = 2.5,
                                               mean_hypnotic_depth_hyp = 5,
                                               mean_hypnotic_depth_plb = 5,
                                               sd_hypnotizability = 2,
                                               mean_hypnotizability = 6,
                                               sd_mean_water_temp = 0.5,
                                               mean_mean_water_temp = 1.5,
                                               which_hypothesis_to_test = "H1.1",
                                               SESOI = -10))

### H1.1 - null true

results_H1.1_nulltrue = pbreplicate(iter, simul(n = 45,
                                                mean_paintolerance_hyp = 75,
                                                mean_paintolerance_plb = 75,
                                                mean_paintolerance_baseline = 50,
                                                sd_paintolerance = 32,
                                                sd_age = 2.5,
                                                mean_age = 22,
                                                sd_height = 14,
                                                mean_height = 175,
                                                sd_expected_CPT_pain = 2,
                                                mean_expected_CPT_pain = 5.5,
                                                sd_FPQ_total = 0.5,
                                                mean_FPQ_total = 3,
                                                sd_expectancy_pain_reduction = 20,
                                                mean_expectancy_pain_reduction_hyp = 50,
                                                mean_expectancy_pain_reduction_plb = 50,
                                                sd_expectancy_hypnotic_depth = 2,
                                                mean_expectancy_hypnotic_depth_hyp = 5,
                                                mean_expectancy_hypnotic_depth_plb = 5,
                                                sd_hypnotic_depth = 2.5,
                                                mean_hypnotic_depth_hyp = 5,
                                                mean_hypnotic_depth_plb = 5,
                                                sd_hypnotizability = 2,
                                                mean_hypnotizability = 6,
                                                sd_mean_water_temp = 0.5,
                                                mean_mean_water_temp = 1.5,
                                                which_hypothesis_to_test = "H1.1",
                                                SESOI = -10))



### H1.2 in expected hypnosis depth

### H1.2 - alternative true

results_H1.2_alttrue = pbreplicate(iter, simul(n = 45,
                                               mean_paintolerance_hyp = 75,
                                               mean_paintolerance_plb = 75,
                                               mean_paintolerance_baseline = 50,
                                               sd_paintolerance = 32,
                                               sd_age = 2.5,
                                               mean_age = 22,
                                               sd_height = 14,
                                               mean_height = 175,
                                               sd_expected_CPT_pain = 2,
                                               mean_expected_CPT_pain = 5.5,
                                               sd_FPQ_total = 0.5,
                                               mean_FPQ_total = 3,
                                               sd_expectancy_pain_reduction = 20,
                                               mean_expectancy_pain_reduction_hyp = 50,
                                               mean_expectancy_pain_reduction_plb = 50,
                                               sd_expectancy_hypnotic_depth = 2,
                                               mean_expectancy_hypnotic_depth_hyp = 5,
                                               mean_expectancy_hypnotic_depth_plb = 4,
                                               sd_hypnotic_depth = 2.5,
                                               mean_hypnotic_depth_hyp = 5,
                                               mean_hypnotic_depth_plb = 5,
                                               sd_hypnotizability = 2,
                                               mean_hypnotizability = 6,
                                               sd_mean_water_temp = 0.5,
                                               mean_mean_water_temp = 1.5,
                                               which_hypothesis_to_test = "H1.2",
                                               SESOI = -1))



### H1.2 - null true

results_H1.2_nulltrue = pbreplicate(iter, simul(n = 45,
                                                mean_paintolerance_hyp = 75,
                                                mean_paintolerance_plb = 75,
                                                mean_paintolerance_baseline = 50,
                                                sd_paintolerance = 32,
                                                sd_age = 2.5,
                                                mean_age = 22,
                                                sd_height = 14,
                                                mean_height = 175,
                                                sd_expected_CPT_pain = 2,
                                                mean_expected_CPT_pain = 5.5,
                                                sd_FPQ_total = 0.5,
                                                mean_FPQ_total = 3,
                                                sd_expectancy_pain_reduction = 20,
                                                mean_expectancy_pain_reduction_hyp = 50,
                                                mean_expectancy_pain_reduction_plb = 50,
                                                sd_expectancy_hypnotic_depth = 2,
                                                mean_expectancy_hypnotic_depth_hyp = 5,
                                                mean_expectancy_hypnotic_depth_plb = 5,
                                                sd_hypnotic_depth = 2.5,
                                                mean_hypnotic_depth_hyp = 5,
                                                mean_hypnotic_depth_plb = 5,
                                                sd_hypnotizability = 2,
                                                mean_hypnotizability = 6,
                                                sd_mean_water_temp = 0.5,
                                                mean_mean_water_temp = 1.5,
                                                which_hypothesis_to_test = "H1.2",
                                                SESOI = -1))


## H2: We will test the hypothesis that conventional hypnosis will evoke greater hypnotic effects than sham hypnosis in the following measures.

### H2.1 Improved pain tolerance compared to baseline (improvement in the conventional hypnosis condition vs. the sham condition)

### H2.1 - alternative true

results_H2.1_alttrue = pbreplicate(iter, simul(n = 45,
                                  mean_paintolerance_hyp = 75,
                                  mean_paintolerance_plb = 65,
                                  mean_paintolerance_baseline = 50,
                                  sd_paintolerance = 32,
                                  sd_age = 2.5,
                                  mean_age = 22,
                                  sd_height = 14,
                                  mean_height = 175,
                                  sd_expected_CPT_pain = 2,
                                  mean_expected_CPT_pain = 5.5,
                                  sd_FPQ_total = 0.5,
                                  mean_FPQ_total = 3,
                                  sd_expectancy_pain_reduction = 20,
                                  mean_expectancy_pain_reduction_hyp = 50,
                                  mean_expectancy_pain_reduction_plb = 50,
                                  sd_expectancy_hypnotic_depth = 2,
                                  mean_expectancy_hypnotic_depth_hyp = 5,
                                  mean_expectancy_hypnotic_depth_plb = 5,
                                  sd_hypnotic_depth = 2.5,
                                  mean_hypnotic_depth_hyp = 5,
                                  mean_hypnotic_depth_plb = 5,
                                  sd_hypnotizability = 2,
                                  mean_hypnotizability = 6,
                                  sd_mean_water_temp = 0.5,
                                  mean_mean_water_temp = 1.5,
                                  which_hypothesis_to_test = "H2.1",
                                  SESOI = -10))



### H2.1 - null true

results_H2.1_nulltrue = pbreplicate(iter, simul(n = 45,
                                     mean_paintolerance_hyp = 75,
                                     mean_paintolerance_plb = 75,
                                     mean_paintolerance_baseline = 50,
                                     sd_paintolerance = 32,
                                     sd_age = 2.5,
                                     mean_age = 22,
                                     sd_height = 14,
                                     mean_height = 175,
                                     sd_expected_CPT_pain = 2,
                                     mean_expected_CPT_pain = 5.5,
                                     sd_FPQ_total = 0.5,
                                     mean_FPQ_total = 3,
                                     sd_expectancy_pain_reduction = 20,
                                     mean_expectancy_pain_reduction_hyp = 50,
                                     mean_expectancy_pain_reduction_plb = 50,
                                     sd_expectancy_hypnotic_depth = 2,
                                     mean_expectancy_hypnotic_depth_hyp = 5,
                                     mean_expectancy_hypnotic_depth_plb = 5,
                                     sd_hypnotic_depth = 2.5,
                                     mean_hypnotic_depth_hyp = 5,
                                     mean_hypnotic_depth_plb = 5,
                                     sd_hypnotizability = 2,
                                     mean_hypnotizability = 6,
                                     sd_mean_water_temp = 0.5,
                                     mean_mean_water_temp = 1.5,
                                     which_hypothesis_to_test = "H2.1",
                                     SESOI = -10))


### H2.2 in Hypnosis depth

### H2.2 - alternative true

results_H2.2_alttrue = pbreplicate(iter, simul(n = 45,
                                               mean_paintolerance_hyp = 75,
                                               mean_paintolerance_plb = 75,
                                               mean_paintolerance_baseline = 50,
                                               sd_paintolerance = 32,
                                               sd_age = 2.5,
                                               mean_age = 22,
                                               sd_height = 14,
                                               mean_height = 175,
                                               sd_expected_CPT_pain = 2,
                                               mean_expected_CPT_pain = 5.5,
                                               sd_FPQ_total = 0.5,
                                               mean_FPQ_total = 3,
                                               sd_expectancy_pain_reduction = 20,
                                               mean_expectancy_pain_reduction_hyp = 50,
                                               mean_expectancy_pain_reduction_plb = 50,
                                               sd_expectancy_hypnotic_depth = 2,
                                               mean_expectancy_hypnotic_depth_hyp = 5,
                                               mean_expectancy_hypnotic_depth_plb = 5,
                                               sd_hypnotic_depth = 2.5,
                                               mean_hypnotic_depth_hyp = 5,
                                               mean_hypnotic_depth_plb = 4,
                                               sd_hypnotizability = 2,
                                               mean_hypnotizability = 6,
                                               sd_mean_water_temp = 0.5,
                                               mean_mean_water_temp = 1.5,
                                               which_hypothesis_to_test = "H2.2",
                                               SESOI = -1))



### H2.2 - null true

results_H2.2_nulltrue = pbreplicate(iter, simul(n = 45,
                                                mean_paintolerance_hyp = 75,
                                                mean_paintolerance_plb = 75,
                                                mean_paintolerance_baseline = 50,
                                                sd_paintolerance = 32,
                                                sd_age = 2.5,
                                                mean_age = 22,
                                                sd_height = 14,
                                                mean_height = 175,
                                                sd_expected_CPT_pain = 2,
                                                mean_expected_CPT_pain = 5.5,
                                                sd_FPQ_total = 0.5,
                                                mean_FPQ_total = 3,
                                                sd_expectancy_pain_reduction = 20,
                                                mean_expectancy_pain_reduction_hyp = 50,
                                                mean_expectancy_pain_reduction_plb = 50,
                                                sd_expectancy_hypnotic_depth = 2,
                                                mean_expectancy_hypnotic_depth_hyp = 5,
                                                mean_expectancy_hypnotic_depth_plb = 5,
                                                sd_hypnotic_depth = 2.5,
                                                mean_hypnotic_depth_hyp = 5,
                                                mean_hypnotic_depth_plb = 5,
                                                sd_hypnotizability = 2,
                                                mean_hypnotizability = 6,
                                                sd_mean_water_temp = 0.5,
                                                mean_mean_water_temp = 1.5,
                                                which_hypothesis_to_test = "H2.2",
                                                SESOI = -1))

#######################################
#        Simulation results           #
#######################################

# Here are the operational characteristics of the study for each hypothesis separately

mean(results_H1.1_alttrue == "H0") # alpha for alternative being true
mean(results_H1.1_alttrue == "H1") # power for detecting alternative being true
mean(results_H1.1_nulltrue == "H0")  # alpha for null being true
mean(results_H1.1_nulltrue == "H1") # power for detecting null being true

mean(results_H1.2_alttrue == "H0")  # alpha for alternative being true
mean(results_H1.2_alttrue == "H1") # power for detecting alternative being true
mean(results_H1.2_nulltrue == "H0") # alpha for null being true
mean(results_H1.2_nulltrue == "H1") # power for detecting null being true

mean(results_H2.1_alttrue == "H0")  # alpha for alternative being true
mean(results_H2.1_alttrue == "H1") # power for detecting alternative being true
mean(results_H2.1_nulltrue == "H0") # alpha for null being true
mean(results_H2.1_nulltrue == "H1") # power for detecting null being true

mean(results_H2.2_alttrue == "H0")  # alpha for alternative being true
mean(results_H2.2_alttrue == "H1") # power for detecting alternative being true
mean(results_H2.2_nulltrue == "H0") # alpha for null being true
mean(results_H2.2_nulltrue == "H1") # power for detecting null being true

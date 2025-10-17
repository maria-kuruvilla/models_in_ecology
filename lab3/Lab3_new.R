# Goal - Use Bayesian hierarchical models to estimate the effect of 
# sea lice on productivity of pink salmon

install.packages("rstan")
library(rstan)
rstan_options("auto_write" = TRUE)
options(mc.cores = parallel::detectCores())

install.packages("tidybayes")
library(tidybayes)

install.packages("here")
library(here)

install.packages("tidyverse")
library(tidyverse)

install.packages("ggplot2")
library(ggplot2)

# data processing ---------------------------------------------------------


data <- read.csv(here("lab3","Pink_S-R_data.csv"))

# get a glimpse of the data
glimpse(data)

# adding columns to data to calculate survival, 
# louse exposure category, and population (odd/even),
# time period

# Infestations occured in spring 2000-2005; better treatment plans in place by spring 2006
# Fallow management intervention in 2003 (affected fish returning in 2004)



data_new <- data %>% 
  mutate(survival = log(Recruits/Spawners),
         population = paste(River,Odd_Even),
         time_period = case_when(Return_Yr < 2002 ~ factor("before"),
                                 Return_Yr == 2004 ~ factor("fallow"),
                                 Return_Yr >= 2002 & Return_Yr <= 2006 ~ factor("during"),
                                 Return_Yr > 2006 ~ factor("after")),
         exposure = as.factor(ifelse(Area==12, "exposed","unexposed"))) %>% 
  mutate(exposure_category = paste(exposure, time_period),
         row_number = row_number()) %>% 
  # add a column with the start row number and another with the end row number for
  # the time series of each population
  group_by(population) %>%
  mutate(time_series_length = n(),
         start_row = first(row_number),
         end_row = last(row_number)) %>% 
  ungroup()
  
  

glimpse(data_new)

#extract max S for priors
smax_prior= data_new %>%
  group_by(population) %>%
  summarize(m.s=Spawners[which.max(Recruits)],m.r=max(Recruits))


# data visualization ------------------------------------------------------



ggplot(data_new, aes(x = Return_Yr, y = survival, group = population))+
  geom_point(alpha = 0.8, aes(color = exposure_category))+
  geom_line(aes(color = exposure_category), alpha = 0.6, size = 1.2)+
  geom_vline(xintercept = 2002, linetype = "dashed", color = "darkred")+
  facet_wrap(~paste("Area:",Area), ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("exposed before" = "lightblue", "exposed during" = "darkred", 
                              "exposed after" = "darkblue", "exposed fallow" = "violet",
                              "unexposed before" = "lightblue", "unexposed during" = "lightblue",
                              "unexposed after" = "lightblue", "unexposed fallow" = "lightblue")) +
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position = "none"
        ) + 
  labs(title = "Pink Salmon Survival vs Return Year",
       x = "Return Year",
       y = "Survival (log(Recruits/Spawners))",
       color = "Exposure Category")

ggplot(data_new %>% filter(Area==12), aes(x = Spawners, y = survival, 
                                         group = population))+
  geom_point(alpha = 0.8, color = "cadetblue")+
  facet_wrap(~population, scales= "free_x")


# stan model 1a -----------------------------------------------------------

# call model

ricker_1a_model <- stan_model(file = here("lab3","ricker_1a.stan"))

# create data list for model 


data_list_1a <- list(
  N = nrow(data_new),
  J = length(unique(data_new$population)),
  spawners = data_new$Spawners,
  survival = data_new$survival,
  Smax_mean = smax_prior$m.s,
  Smax_sigma = 2*smax_prior$m.s,
  start_row = unique(data_new$start_row),
  end_row = unique(data_new$end_row)
)

# sample from model
ricker_1a_model_sampling <- rstan::sampling(ricker_1a_model,
                                            data = data_list_1a,
                                            chains = 4,
                                            iter = 2000,
                                            warmup = 1000)
                                            

# plot some model results

# check trace plots to look for caterpillar like sampling plots

bayesplot::mcmc_trace(ricker_1a_model_sampling, pars = "alpha_j[10]")
bayesplot::mcmc_trace(ricker_1a_model_sampling, pars = "b[1]")
bayesplot::mcmc_trace(ricker_1a_model_sampling, pars = "Smax[1]")

# plot posterior probability distribution
bayesplot::mcmc_areas(ricker_1a_model_sampling, pars = c("alpha_j[1]"))
bayesplot::mcmc_areas(ricker_1a_model_sampling, pars = c("Smax[11]"))

# plot all the alphas from area 12
population_numbers <- data_new %>%
  # group_by(population) %>%
  #make row for population number
  distinct(population, .keep_all = TRUE) %>%
  mutate(population_n = row_number()) %>% 
  select(population, population_n, Area) 

area_12 <- population_numbers %>% 
    filter(Area == 12)

area_7 <- population_numbers %>% 
  filter(Area == 7)


alpha_j_area_12 <- paste0("alpha_j[",area_7$population_n,"]")

bayesplot::mcmc_areas(ricker_1a_model_sampling,
                      pars = alpha_j_area_12,
                      prob = 0.8) +
  labs(title = "Posterior Distributions of alpha for all populations of Area 12",
       x = "Alpha",
       y = "Density")



# data_processing for 1b --------------------------------------------------

# make dummy variable for exposure category

exposure_category_matrix <- model.matrix(~0 + exposure_category, 
                                         data = data_new)                                    
#make data list for 1b

data_list_1b <- list(
  N = nrow(data_new),
  K = ncol(exposure_category_matrix),
  J = length(unique(data_new$population)),
  spawners = data_new$Spawners,
  survival = data_new$survival,
  exposure_category_matrix = exposure_category_matrix,
  Smax_mean = smax_prior$m.s,
  Smax_sigma = 2*smax_prior$m.s,
  start_row = unique(data_new$start_row),
  end_row = unique(data_new$end_row)
)

# stan model 1b -----------------------------------------------------------

# call model

ricker_1b_model <- stan_model(file = here("lab3","ricker_1b.stan"))

# sample from model

ricker_1b_model_sampling <- rstan::sampling(ricker_1b_model,
                                            data = data_list_1b,
                                            chains = 4,
                                            iter = 2000,
                                            warmup = 1000)

bayesplot::mcmc_trace(ricker_1b_model_sampling, pars = "alpha_j[10]")
bayesplot::mcmc_trace(ricker_1b_model_sampling, pars = "b[1]")
bayesplot::mcmc_trace(ricker_1b_model_sampling, pars = "Smax[1]")


# plot posterior probability distribution
bayesplot::mcmc_areas(ricker_1b_model_sampling, pars = c("alpha_j[1]"))
bayesplot::mcmc_areas(ricker_1b_model_sampling, pars = c("alpha_j[11]"))
bayesplot::mcmc_areas(ricker_1b_model_sampling, pars = c("Smax[11]"))
bayesplot::mcmc_areas(ricker_1b_model_sampling, pars = c("beta[6]"))

exposure_effects_list <- paste0("beta[",1:8,"]")

bayesplot::mcmc_areas(ricker_1b_model_sampling,
                      pars = exposure_effects_list,
                      prob = 0.8) +
  labs(title = "Posterior Distributions of effects of exposure",
       x = "Alpha",
       y = "Density")



# model 3a ----------------------------------------------------------------

# include a random effect of year

#need then number of years

years = max(data_new$Return_Yr) - min(data_new$Return_Yr) + 1

data_list_3a <- list(
  N = nrow(data_new),
  J = length(unique(data_new$population)),
  years = years,
  ii = as.numeric(factor(data_new$Return_Yr)),  #index of brood years
  spawners = data_new$Spawners,
  survival = data_new$survival,
  Smax_mean = smax_prior$m.s,
  Smax_sigma = 2*smax_prior$m.s,
  start_row = unique(data_new$start_row),
  end_row = unique(data_new$end_row)
)

# call model

ricker_3a_model <- stan_model(file = here("lab3","ricker_3a.stan"))

# sample from model

ricker_3a_model_sampling <- rstan::sampling(ricker_3a_model,
                                            data = data_list_3a,
                                            chains = 4,
                                            iter = 2000,
                                            warmup = 1000)

bayesplot::mcmc_trace(ricker_3a_model_sampling, pars = "alpha")
bayesplot::mcmc_trace(ricker_3a_model_sampling, pars = "b[1]")
bayesplot::mcmc_trace(ricker_3a_model_sampling, pars = "Smax[1]")

# plot posterior probability distribution
bayesplot::mcmc_areas(ricker_3a_model_sampling, pars = c("alpha"))

#plot theta

year_effects_list <- paste0("theta_year[",1:years,"]")

bayesplot::mcmc_areas(ricker_3a_model_sampling,
                      pars = year_effects_list,
                      prob = 0.8) +
  labs(title = "Posterior Distributions of year effects",
       x = "Year effect",
       y = "Density")

#extract alpha, theta_year 

post_3a = ricker_3a_model_sampling %>% 
  rstan::extract(pars = c("alpha","theta_year","Smax"), permuted = TRUE)

glimpse(post_3a$alpha)

# plot median theta_year as a time series with credible intervals



theta_year_df <- data.frame(Return_Yr = min(data_new$Return_Yr):max(data_new$Return_Yr),
                            theta_year = apply(post_3a$theta_year, 2, median),
                            theta_year_lower = apply(post_3a$theta_year, 2, quantile, probs = 0.025),
                            theta_year_upper = apply(post_3a$theta_year, 2, quantile, probs = 0.975)
                            )

#plot

ggplot(theta_year_df, aes(x = Return_Yr, y = theta_year))+
  geom_line(color = "#4E654D", size = 1.2)+
  geom_ribbon(aes(ymin = theta_year_lower, ymax = theta_year_upper), 
              alpha = 0.5, fill = "#9CC69B")+
  theme_classic()+
  labs(title = "Median Year Effects from Model 3a",
       x = "Return Year",
       y = "Theta Year Effect")







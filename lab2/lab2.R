# lab 2

library(here)
library(tidyverse)
library(deSolve)

#data processing

# read .txt file

p_aurelia <- read.table(here("lab2","P_aurelia.txt"), header = F)
p_caudatum <- read.table(here("lab2","P_caudatum.txt"), header = F)

colnames(p_aurelia) <- c("day", "number")
colnames(p_caudatum) <- c("day", "number")

p_aurelia$species <- "p_aurelia"
p_caudatum$species <- "p_caudatum"

monoculture_data <- rbind(p_aurelia, p_caudatum)

competition_data <- read.table(here("lab2","Competition - Copy.txt"), header = F)

colnames(competition_data) <- c("day","number_p_aurelia","number_p_caudatum")

competition_data_long <- competition_data %>%
  pivot_longer(cols = c("number_p_aurelia", "number_p_caudatum"),
               names_to = "species",
               values_to = "number",
               names_prefix = "number_") %>% 
  mutate(treatment = "competition")

monoculture_data$treatment <- "monoculture"

data <- rbind(monoculture_data, competition_data_long)

n0 <- c(n=data %>% filter(day==0, species == "p_aurelia", treatment == "monoculture") %>% select(number))
tf <- length(p_aurelia[,1])

#calculate growth

growth = numeric(tf-1)

for(i in seq_along(growth)){
  growth[i] = (p_aurelia$number[i+1]-p_aurelia$number[i])/((p_aurelia$day[i+1] - p_aurelia$day[i])*p_aurelia$number[i])
}

abundance = p_aurelia[-length(p_aurelia[,1]),] %>% 
  mutate(growth = growth)

#plot
ggplot(abundance)+
  geom_point(aes(x=number, y=growth), alpha = 0.2, size = 2)+
  labs(x="Abundance", y="Growth rate")+
  theme_classic()







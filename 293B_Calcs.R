library(coda)
library(rjags)
library(BEST)
library(ggplot2)
library(tidyverse)

setwd("D:/Dissertation")

my_data <- read_csv('293B_Rdata.csv')

#Data Organization
Exp_Group <- my_data[my_data$Group == 'Exp',]
Con_Group <- my_data[my_data$Group == 'Con',]

Exp_hGPA <- Exp_Group[['HS_GPA']]
Exp_hGPA <- na.omit(Exp_hGPA)
Con_hGPA <- Con_Group[['HS_GPA']]
Con_hGPA <- na.omit(Con_hGPA)

Exp_SAT <- Exp_Group[['M_SAT']]
Exp_SAT <- na.omit(Exp_SAT)
Con_SAT <- Con_Group[['M_SAT']]
Con_SAT <- na.omit(Con_SAT)

Exp_ACT <- Exp_Group[['M_ACT']]
Exp_ACT <- na.omit(Exp_ACT)
Con_ACT <- Con_Group[['M_ACT']]
Con_ACT <- na.omit(Con_ACT)

Exp_cGPA <- Exp_Group[['cGPA']]
Exp_cGPA <- na.omit(Exp_cGPA)
Con_cGPA <- Con_Group[['cGPA']]
Con_cGPA <- na.omit(Con_cGPA)

Exp_s23GPA <- Exp_Group[['s23GPA']]
Exp_s23GPA <- na.omit(Exp_s23GPA)
Con_s23GPA <- Con_Group[['s23GPA']]
Con_s23GPA <- na.omit(Con_s23GPA)

Exp_dMAI <- Exp_Group[['dMAI']]
Exp_dMAI <- na.omit(Exp_dMAI)

#hGPA

muM_hGPA <- 3
muSD_hGPA <- 1

priors_hGPA <- list(muM = muM_hGPA, muSD = muSD_hGPA)

BESThGPA <- BESTmcmc(Exp_hGPA, Con_hGPA, priors=priors_hGPA, parallel=FALSE)

plot(BESThGPA)
#plotAll(BESThGPA)

#SAT

muM_SAT <- 550
muSD_SAT <- 30

priors_SAT <- list(muM = muM_SAT, muSD = muSD_SAT)

BESTSAT <- BESTmcmc(Exp_SAT, Con_SAT, priors=priors_SAT, parallel=FALSE)

plot(BESTSAT)
#plotAll(BESTSAT)

#ACT

muM_ACT <- 22
muSD_ACT <- 2

priors_ACT <- list(muM = muM_ACT, muSD = muSD_ACT)

BESTACT <- BESTmcmc(Exp_ACT, Con_ACT, priors=priors_ACT, parallel=FALSE)

plot(BESTACT)
#plotAll(BESTACT)

#cGPA

muM_cGPA <- 3
muSD_cGPA <- 1

priors_cGPA <- list(muM = muM_cGPA, muSD = muSD_cGPA)

BESTcGPA <- BESTmcmc(Exp_cGPA, Con_cGPA, priors=priors_cGPA, parallel=FALSE)

plot(BESTcGPA)
#plotAll(BESTcGPA)

#s23GPA

muM_s23GPA <- 3
muSD_s23GPA <- 1

priors_s23GPA <- list(muM = muM_s23GPA, muSD = muSD_s23GPA)

BESTs23GPA <- BESTmcmc(Exp_s23GPA, Con_s23GPA, priors=priors_s23GPA, parallel=FALSE)

plot(BESTs23GPA)
#plotAll(BESTs23GPA)

muM_dMAI <- 52
muSD_dMAI <- 495

priors_dMAI <- list(muM = muM_dMAI, muSD = muSD_dMAI)

BESTdMAI <- BESTmcmc(Exp_dMAI, priors=priors_dMAI, parallel=FALSE)

plot(BESTdMAI)
#plotAll(BESTdMAI)

#___________________________________________________________________
#Proportion Testing

library(tidyverse)    # ggplot, dplyr, and friends
library(gt)           # Fancy tables
library(glue)         # Easier string interpolation
library(scales)       # Nicer labeling functions
library(ggmosaic)     # Mosaic plots with ggplot
library(ggpattern)    # Pattern fills in ggplot
library(patchwork)    # Combine plots nicely
library(parameters)   # Extract model parameters as data frames
library(cmdstanr)     # Run Stan code from R
library(brms)         # Nice frontend for Stan
library(tidybayes)    # Manipulate MCMC chains in a tidy way
library(likert)       # Contains the pisaitems data
library(jtools)

options(brms.backend = "cmdstanr")

clrs_saguaro <- NatParksPalettes::natparks.pals("Saguaro")
clr_Exp <- clrs_saguaro[1]
clr_Con <- clrs_saguaro[6]
clr_Diff <- clrs_saguaro[4]

CHAINS <- 4
ITER <- 3000
WARMUP <- 1500
BAYES_SEED <- 1231

# Male % 

M_counts <- my_data %>%
  group_by(Group, Male) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

M_values <- M_counts %>%
  filter(Male == 1)

M_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = M_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  #need to up the warm up and iterations because there seemed to be a
  #convergence issue
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


Mp1 <- M_model %>% 
  epred_draws(newdata = M_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of male students in each group",
       y = NULL) +
  theme_nice()

Mp2 <- M_model %>% 
  epred_draws(newdata = M_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "pp." isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the pointrange doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(Mp1 / plot_spacer() / Mp2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

M_delta <- M_model %>% 
  epred_draws(newdata = M_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)

M_counts
M_delta

# Pass Math

Pmath_counts <- my_data %>%
  group_by(Group, Math_Pass) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

Pmath_values <- Pmath_counts %>%
  filter(Math_Pass == 1)

Pmath_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = Pmath_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  #need to up the warm up and iterations because there seemed to be a
  #convergence issue
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


Pmathp1 <- Pmath_model %>% 
  epred_draws(newdata = Pmath_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of students who pass Trig",
       y = NULL) +
  theme_nice()

Pmathp2 <- Pmath_model %>% 
  epred_draws(newdata = Pmath_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "pp." isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the pointrange doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(Pmathp1 / plot_spacer() / Pmathp2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

Pmath_delta <- Pmath_model %>% 
  epred_draws(newdata = Pmath_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)

Pmath_counts
Pmath_delta

# A in Math

Amath_counts <- my_data %>%
  group_by(Group, Math_A) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

Amath_values <- Amath_counts %>%
  filter(Math_A == 1)

Amath_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = Amath_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  #need to up the warm up and iterations because there seemed to be a
  #convergence issue
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


Amathp1 <- Amath_model %>% 
  epred_draws(newdata = Amath_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of students who receive an A in Trig",
       y = NULL) +
  theme_nice()

Amathp2 <- Amath_model %>% 
  epred_draws(newdata = Amath_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "pp." isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the pointrange doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(Amathp1 / plot_spacer() / Amathp2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

Amath_delta <- Amath_model %>% 
  epred_draws(newdata = Amath_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)

Amath_counts
Amath_delta

# Pass chem

chem_data <- my_data[my_data$Took_Chem == 1,]

Pchem_counts <- chem_data %>%
  group_by(Group, Chem_Pass) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

Pchem_values <- Pchem_counts %>%
  filter(Chem_Pass == 1)

Pchem_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = Pchem_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  #need to up the warm up and iterations because there seemed to be a
  #convergence issue
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


Pchemp1 <- Pchem_model %>% 
  epred_draws(newdata = Pchem_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of students who pass Chem",
       y = NULL) +
  theme_nice()

Pchemp2 <- Pchem_model %>% 
  epred_draws(newdata = Pchem_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "pp." isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the pointrange doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(Pchemp1 / plot_spacer() / Pchemp2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

Pchem_delta <- Pchem_model %>% 
  epred_draws(newdata = Pchem_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)

Pchem_counts
Pchem_delta

# A chem

Achem_counts <- chem_data %>%
  group_by(Group, Chem_A) %>%
  summarize(n = n()) %>%
  mutate(total = sum(n),
         prop = n/sum(n)) %>%
  ungroup()

Achem_values <- Achem_counts %>%
  filter(Chem_A == 1)

Achem_model <- brm(
  bf(n | trials(total) ~ 0 + Group),
  data = Achem_values,
  family = beta_binomial(link = "identity", link_phi= "identity"),
  prior = c(prior(beta(1, 1), class = "b", dpar = "mu", lb = 0, ub = 1),
            prior(exponential(0.001), class = "phi", lb = 0)),
  #need to up the warm up and iterations because there seemed to be a
  #convergence issue
  chains = CHAINS, warmup = WARMUP, iter = ITER, seed = BAYES_SEED,
  refresh = 0
)


Achemp1 <- Achem_model %>% 
  epred_draws(newdata = Achem_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ggplot(aes(x = .epred_prop, y = Group, fill = Group)) +
  stat_halfeye() +
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.015)) +
  scale_fill_manual(values = c(clr_Con, clr_Exp)) +
  guides(fill = "none") +
  labs(x = "Proportion of students who receive an A in Chem",
       y = NULL) +
  theme_nice()

Achemp2 <- Achem_model %>% 
  epred_draws(newdata = Achem_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>% 
  ggplot(aes(x = .epred_prop)) +
  stat_halfeye(fill = clr_Diff) +
  # Multiply axis limits by 0.5% so that the right "pp." isn't cut off
  scale_x_continuous(labels = label_percent(), expand = c(0, 0.005)) +
  labs(x = "Percentage point difference in proportions",
       y = NULL) +
  # Make it so the pointrange doesn't get cropped
  coord_cartesian(clip = "off") + 
  theme_nice() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

(Achemp1 / plot_spacer() / Achemp2) + 
  plot_layout(heights = c(0.785, 0.03, 0.185))

Achem_delta <- Achem_model %>% 
  epred_draws(newdata = Achem_values) %>% 
  mutate(.epred_prop = .epred / total) %>% 
  ungroup() %>% 
  mutate(Group = fct_relevel(Group, "Exp")) %>% 
  compare_levels(.epred_prop, by = Group,
                 comparison = "pairwise") %>%
  summarize(median = median_qi(.epred_prop, .width = 0.95),
            p_gt_0 = sum(.epred_prop > 0) / n()) %>% 
  unnest(median)

Achem_counts
Achem_delta

#_______________________________________________________________________
#Metacognition Questions from Exams

#comparison of metacognition pre and post test
metaacc_data <- read_csv('meta_acc.csv')

#Easy Content Question
q1_dpost <- metaacc_data[['d_post_acc1']]
q1_dpre <- metaacc_data[['d_pre_acc1']]

#Hard Content Question
q2_dpost <- metaacc_data[['d_post_acc2']]
q2_dpre <- metaacc_data[['d_pre_acc2']]

muM_dacc <- 7
muSD_dacc <- 21

priors_dacc <- list(muM = muM_dacc, muSD = muSD_dacc)

BESTq1dpost <- BESTmcmc(q1_dpost, priors=priors_dacc, parallel=FALSE)
BESTq1dpre <- BESTmcmc(q1_dpre, priors=priors_dacc, parallel=FALSE)

BESTq2dpost <- BESTmcmc(q2_dpost, priors=priors_dacc, parallel=FALSE)
BESTq2dpre <- BESTmcmc(q2_dpre, priors=priors_dacc, parallel=FALSE)

plot(BESTq1dpost)
plot(BESTq1dpre)
plot(BESTq2dpost)
plot(BESTq2dpre)

#comparison of metacognition pre test (Midterm) and final exam 
#to examine impact of student preparation on their calibration

metastudy_data <- read_csv('studymeta.csv')

#Easy Content
q1_dpostS <- metastudy_data[['dpost1']]
q1_dpreS <- metastudy_data[['dpre1']]

#Difficult Content
q2_dpostS <- metastudy_data[['dpost2']]
q2_dpreS <- metastudy_data[['dpre2']]

muM_dacc <- -3
muSD_dacc <- 19

priors_dstudy <- list(muM = muM_dacc, muSD = muSD_dacc)

BESTq1dpostS <- BESTmcmc(q1_dpostS, priors=priors_dstudy, parallel=FALSE)
BESTq1dpreS <- BESTmcmc(q1_dpreS, priors=priors_dstudy, parallel=FALSE)

BESTq2dpostS <- BESTmcmc(q2_dpostS, priors=priors_dstudy, parallel=FALSE)
BESTq2dpreS <- BESTmcmc(q2_dpreS, priors=priors_dstudy, parallel=FALSE)

plot(BESTq1dpostS)
plot(BESTq1dpreS)
plot(BESTq2dpostS)
plot(BESTq2dpreS)
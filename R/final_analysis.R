# Final analysis script

# packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(forcats)
library(purrr)
library(broom)
library(emmeans)
library(lme4)
library(car)
library(MASS)
library(insight)
library(ggeffects)

# read in city dat
final_city_level_dat <- readRDS("Data/city_dat.RDS")

# summarize some stuff for results
length(unique(final_city_level_dat$city))
length(unique(final_city_level_dat$COMMON_NAME))
sum(final_city_level_dat$number_obs)

# read in region dat
final_region_level_dat <- readRDS("Data/region_dat.RDS")

# summarize some stuff for results
length(unique(final_region_level_dat$city))
length(unique(final_region_level_dat$COMMON_NAME))
final_region_level_dat %>%
  dplyr::select(region, COMMON_NAME, number_obs, urban_score) %>%
  distinct() %>%
  .$number_obs %>%
  sum()
sum(final_region_level_dat$number_obs)



# Now look at individual cities
# now look at individual cities
ggplot(final_city_level_dat, aes(x=UrbanToleranceIndex, y=urban_score, color=region))+
  geom_point()+
  scale_y_log10()+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  facet_wrap(~city, scales="free", ncol=5)+
  theme(strip.text = element_text(size = 8))+
  theme(axis.text=element_text(color="black"))+
  guides(color=FALSE)+
  geom_smooth(method="lm")+
  xlab("Urban Tolerance Abundance Index")+
  ylab("Citizen science urbanness measure")

#ggsave("Figures/individual_city_level_scores.png", width=8.5, height=8, units="in")

# apply a linear regression to every city
city_regressions <- final_city_level_dat %>%
  nest(-city) %>%
  mutate(
    fit=map(data, ~ lm(log10(urban_score) ~ sqrt(UrbanToleranceIndex+8.1765), data=.x)),
    tidied=map(fit, tidy),
    glanced=map(fit, glance)
  )

regression_summaries <- city_regressions %>%
  unnest(glanced) %>%
  dplyr::select(1, 5:15) %>%
  left_join(., final_city_level_dat %>%
              dplyr::select(city, total_checklists) %>%
              distinct())

#write_csv(regression_summaries, "regression_summaries_with_local_scores.csv")

# summarize these city-specific regressions for the paper
# does it make sense to do it this way?
min(regression_summaries$r.squared)
max(regression_summaries$r.squared)
mean(regression_summaries$r.squared)

# now look at individual cities
ggplot(final_region_level_dat, aes(x=UrbanToleranceIndex, y=urban_score, group=city, color=region))+
  geom_point()+
  scale_y_log10()+
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  facet_wrap(~city, scales="free", ncol=5)+
  theme(strip.text = element_text(size = 8))+
  theme(axis.text=element_text(color="black"))+
  guides(color=FALSE)+
  geom_smooth(method="lm")+
  xlab("Urban tolerance abundance index")+
  ylab("Citizen science urbanness measure")

#ggsave("Figures/individual_city_level_scores_using_regional_scores.png", width=8.5, height=8, units="in")

# test the relationship between city and regional scores
# apply a linear regression to every city
city_regressions_regional_scores <- final_region_level_dat %>%
  nest(-city) %>%
  mutate(
    fit=map(data, ~ lm(log10(urban_score) ~ UrbanToleranceIndex, data=.x)),
    tidied=map(fit, tidy),
    glanced=map(fit, glance)
  )

regression_summaries_regional_scores <- city_regressions_regional_scores %>%
  unnest(glanced) %>%
  dplyr::select(1, 5:15)

#write_csv(regression_summaries_regional_scores, "regression_summaries_with_regional_scores.csv")

# summarize these city-specific regressions for the paper
min(regression_summaries_regional_scores$r.squared)
max(regression_summaries_regional_scores$r.squared)
mean(regression_summaries_regional_scores$r.squared)


# does the R2 change between the two approaches
regression_summary <- regression_summaries_regional_scores %>%
  left_join(regression_summaries, by="city") %>%
  dplyr::select(city, r.squared.x, r.squared.y) %>%
  rename(local_correlation=r.squared.y,
         regional_correlation=r.squared.x)

#write_csv(regression_summary, "regression_summary.csv")

# Look at the regression summaries a bit more
# first the r2 versus the level of urbanization at a city-level
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# categorical model analysis
# simple version
# not sure how to do this with a random effect for city...
# but I think it is probably fine without this
# as this really is simple a simplified version of the continuous plot
final_city_level_dat %>%
  group_by(categorical_index) %>%
  summarize(mean=mean(urban_score),
            sd=sd(urban_score))

final_city_level_dat %>%
  lm(log10(urban_score) ~ categorical_index, data=.) -> mod

mod %>%
  summary()

model_emm <- emmeans(mod, "categorical_index")
eff_size(model_emm, sigma=sigma(mod), edf=1000)

# categorical model analysis
# simple version
# not sure how to do this with a random effect for city...
# but I think it is probably fine without this
# as this really is simple a simplified version of the continuous plot
final_region_level_dat %>%
  group_by(categorical_index) %>%
  summarize(mean=mean(urban_score),
            sd=sd(urban_score))

final_region_level_dat %>%
  lm(log10(urban_score) ~ categorical_index, data=.) -> mod

mod %>%
  summary()

model_emm <- emmeans(mod, "categorical_index")
eff_size(model_emm, sigma=sigma(mod), edf=1000)

####
### Bayesian Phylogenetic Models (MCMCglmm)
####

library(MCMCglmm)
library(ape)

# Open the city data:
# Open and prune phylogenetic tree:
tree <- read.nexus("Data/AllBirdsEricson1_summary.tre")

##
### Repeat City-level analyses including phylogenetic effects:
##

# We convert to dataframe (tibble does not work well with MCMCglmm)
names(final_city_level_dat)
dat1 <- data.frame(final_city_level_dat[,c(1,2,4,7,13,17)])
#dat1$UrbanToleranceIndex <- dat1$UrbanToleranceIndex+8.1765
# Species names need to be called "animal"
names(dat1)[2] <- "animal"

# Prune bird tree to our species
tree1 <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% dat1$animal])

# To control number of iterations, burnin, etc.
Nnitt=1010000
Nthin=1000
Nburnin=10000
(Nnitt-Nburnin)/Nthin #Sample

prior1 <-list(R = list(V = 1, nu = 0.002),G = list(G1 = list(V = 1, nu = 0.002),G2 = list(V = 1, nu = 0.002)))

mod1a <- MCMCglmm(log10(urban_score) ~ UrbanToleranceIndex, random=~animal+city, data=dat1, pedigree=tree1, family="gaussian", verbose=T,prior=prior1,nitt=Nnitt,burnin=Nburnin,thin=Nthin)
summary(mod1a)
mod1a$DIC

mod1a_int_only <- MCMCglmm(log10(urban_score) ~ 1, random=~animal+city, data=dat1, pedigree=tree1, family="gaussian", verbose=T,prior=prior1,nitt=Nnitt,burnin=Nburnin,thin=Nthin)
summary(mod1a_int_only)

# explanation of phylogenetic signal
# animal
plot(mod1a_int_only$VCV[,1]/(rowSums(mod1a_int_only$VCV)))
mean(mod1a_int_only$VCV[,1]/(rowSums(mod1a_int_only$VCV)))
Rmisc::CI(mod1a_int_only$VCV[,1]/(rowSums(mod1a_int_only$VCV)), 0.95)
# city
plot(mod1a_int_only$VCV[,2]/(rowSums(mod1a_int_only$VCV)))
mean(mod1a_int_only$VCV[,2]/(rowSums(mod1a_int_only$VCV)))
Rmisc::CI(mod1a_int_only$VCV[,2]/(rowSums(mod1a_int_only$VCV)), 0.95)

# Check for autocorrrelation of samples
autocorr.diag(mod1a$Sol) # Fixed effects
autocorr.diag(mod1a$VCV) # Random effects

rnd.eff <- apply(mod1a$VCV,2,mean) #Random effects
rnd.eff[1]/sum(rnd.eff) # Effect of phylogeny
rnd.eff[2]/sum(rnd.eff) # Effect of city
rnd.eff[3]/sum(rnd.eff) # Remaining variance (Error)


# Test library ggeffects package
marginal.effects.city <- ggpredict(mod1a, terms="UrbanToleranceIndex")

city_marg <- ggplot(marginal.effects.city, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Urban tolerance abundance index")+
  ylab("Citizen science urbanness score (log10)")+
  ggtitle("a)")

city_marg



##
### Repeat Region-level analyses including phylogenetic effects:
##

# We convert to dataframe (tibble does not work well with MCMCglmm)
names(final_region_level_dat)
dat2 <- data.frame(final_region_level_dat[,c(1,2,7,13,15)])
#dat2$UrbanToleranceIndex <- dat2$UrbanToleranceIndex+8.1765
# Species names need to be called "animal"
names(dat2)[2] <- "animal"

# Prune bird tree to our species
tree2 <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% dat2$animal])

prior2 <-list(R = list(V = 1, nu = 0.002),G = list(G1 = list(V = 1, nu = 0.002),G2 = list(V = 1, nu = 0.002)))

mod2a <- MCMCglmm(log10(urban_score) ~ UrbanToleranceIndex, random=~animal+city, data=dat2, pedigree=tree2, family="gaussian", verbose=T,prior=prior2,nitt=Nnitt,burnin=Nburnin,thin=Nthin)
summary(mod2a)

mod2a_int_only <- MCMCglmm(log10(urban_score) ~ 1, random=~animal+city, data=dat2, pedigree=tree2, family="gaussian", verbose=T,prior=prior2,nitt=Nnitt,burnin=Nburnin,thin=Nthin)
summary(mod2a_int_only)

# explanation of phylogenetic signal
# animal
plot(mod2a_int_only$VCV[,1]/(rowSums(mod2a_int_only$VCV)))
mean(mod2a_int_only$VCV[,1]/(rowSums(mod2a_int_only$VCV)))
Rmisc::CI(mod2a_int_only$VCV[,1]/(rowSums(mod2a_int_only$VCV)), 0.95)
# city
plot(mod2a_int_only$VCV[,2]/(rowSums(mod2a_int_only$VCV)))
mean(mod2a_int_only$VCV[,2]/(rowSums(mod2a_int_only$VCV)))
Rmisc::CI(mod2a_int_only$VCV[,2]/(rowSums(mod2a_int_only$VCV)), 0.95)

# Check for autocorrrelation of samples
autocorr.diag(mod2a$Sol) # Fixed effects
autocorr.diag(mod2a$VCV) # Random effects

rnd.eff <- apply(mod2a$VCV,2,mean) #Random effects
rnd.eff[1]/sum(rnd.eff) # Effect of phylogeny
rnd.eff[2]/sum(rnd.eff) # Effect of city
rnd.eff[3]/sum(rnd.eff) # Remaining variance (Error)

marginal.effects.region <- ggpredict(mod2a, terms="UrbanToleranceIndex")

region_marg <- ggplot(marginal.effects.region, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Urban tolerance abundance index")+
  ylab("Citizen science urbanness score (log10)")+
  ggtitle("b)")

region_marg

city_marg + region_marg + plot_layout(ncol=1)

#ggsave("Figures/marginal_effects_plots.png", height=7, width=4.5, units="in")


# plot the relationship between local level urban score and 
# the local level urban tolerance index
# first on a continuous scale
continuous_plot_city <- ggplot()+
  geom_point(data=final_city_level_dat, aes(x=UrbanToleranceIndex, y=urban_score), color="#377EB8")+
  scale_y_log10()+
  #geom_line(data=marginal.effects.city, aes(x=x, y=10^predicted), size=4)+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  guides(color=FALSE)+
  geom_abline(slope = summary(mod1a)$sol[2], intercept = summary(mod1a)$sol[1], color="#E41A1C", size=1)+
  geom_abline(slope = summary(mod1a)$sol[4], intercept = summary(mod1a)$sol[3], color="gray10", size=1, linetype="dashed")+
  geom_abline(slope = summary(mod1a)$sol[6], intercept = summary(mod1a)$sol[5], color="gray10", size=1, linetype="dashed")+
  xlab("Local-scale abundance-based urban tolerance")+
  ylab("Citizen science urbanness score")

continuous_plot_city
#ggsave("Figures/continuous_plot_city.png", width=5, height=5, units="in")


# then on a categorical scale
categorical_plot_city <- ggplot(final_city_level_dat, aes(x=factor(categorical_index, level=c("UrbanAbsent", "WildIncrease", "UrbanIncrease", "WildAbsent")),
                                                          y=urban_score, fill=categorical_index))+
  geom_violin(width=1)+
  geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  scale_fill_brewer(palette ="Set1")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  guides(color=FALSE)+
  xlab("")+
  ylab("Citizen science urbanness score")+
  #coord_flip()+
  guides(fill=FALSE)

categorical_plot_city
#ggsave("Figures/categorical_plot_city.png", dpi=600, width=5, height=5, units="in")


# make plots as above
continuous_plot_region <- ggplot(final_region_level_dat, aes(x=UrbanToleranceIndex, y=urban_score))+
  geom_point(color="#377EB8")+
  scale_y_log10()+
  theme_bw()+
  #geom_line(data=marginal.effects.region, aes(x=x, y=10^predicted), size=4)+
  theme(axis.text=element_text(color="black"))+
  guides(color=FALSE)+
  geom_abline(slope = summary(mod2a)$sol[2], intercept = summary(mod2a)$sol[1], color="#E41A1C", size=1)+
  geom_abline(slope = summary(mod2a)$sol[4], intercept = summary(mod2a)$sol[3], color="gray10", size=1, linetype="dashed")+
  geom_abline(slope = summary(mod2a)$sol[6], intercept = summary(mod2a)$sol[5], color="gray10", size=1, linetype="dashed")+
  xlab("Local-scale abundance-based urban tolerance")+
  ylab("Citizen science urbanness score")

continuous_plot_region
#ggsave("Figures/continuous_plot_region.png", width=5, height=5, units="in")

categorical_plot_region <- ggplot(final_region_level_dat, aes(x=factor(categorical_index, level=c("UrbanAbsent", "WildIncrease", "UrbanIncrease", "WildAbsent")),
                                                              y=urban_score, fill=categorical_index))+
  geom_violin(width=1)+
  geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  scale_fill_brewer(palette ="Set1")+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  guides(color=FALSE)+
  xlab("")+
  ylab("Citizen science urbanness score")+
  #coord_flip()+
  guides(fill=FALSE)

categorical_plot_region
#ggsave("Figures/continuous_plot_region.png", width=5, height=5, units="in")

# make a figure 2 by putting both city and region level scores together
continuous_plot_city + ggtitle("a) Regional-level urban scores (N=771 species)") + categorical_plot_city + ylab("") +
  continuous_plot_region + ggtitle("b) Continental-level urban scores (N=934 species)") + categorical_plot_region + ylab("")+
  plot_layout(ncol=2, nrow=2)
#ggsave("Figures/figure_2.png", height=7, width=8.5, units="in")

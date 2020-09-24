library(lodr)
library(tidyverse)
library(xtable)

# Set sample size and sim number
n <- 100
sims <- 100

# Set containers for sim fits
sim_fits <- list()
sim_fits[[1]] <- list()
sim_fits[[2]] <- list()
sim_fits[[3]] <- list()

names(sim_fits) <- c("lodr", "sub_lod", "complete_case")

# Run analyses
for(i in 1:sims){
  # Sim data
  print(i)
  set.seed(i)
  n <- 100
  x0 <- rep(1, n); x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  y <- 1*x0+1*x1+1*x2+1*x3+rnorm(n)
  x1 <- ifelse(x1<0, 0, x1); x2 <- ifelse(x2<0, 0, x2); x3 <- ifelse(x3<0, 0, x3)
  sim_data <- data.frame(y,x1,x2,x3)
  
  sim_fits[["lodr"]][[i]] <- lod_lm(data=lod_data_ex, frmla=y~x1+x2+x3, lod=c(0,0),
                               var_LOD=c("x2", "x3"), nSamples=100, boots=25)
  
  sim_fits[["sub_lod"]][[i]] <- lm(data=lod_data_ex, formula=y~x1+x2+x3)
  
  data_complete <- sim_data %>%
    mutate(x2=ifelse(x2==0, NA, x2),
           x3=ifelse(x3==0, NA, x3)) %>%
    filter(complete.cases(.))
  sim_fits[["complete_case"]][[i]] <- lm(data=data_complete, formula=y~x1+x2+x3)
}

# Compile into Mean+-SE plot
results_list <- list()
for(i in 1:sims){
  print(i)
  data_frame_lodr <- rownames_to_column(data.frame(summary(sim_fits[["lodr"]][[i]])$coefficients),
                                        var = "covariate") %>%
    mutate(sim=i,
           method="lodr",
           df=summary(sim_fits[["lodr"]][[i]])$df[2]) %>%
    select(covariate, Estimate, Std..Error, df, sim, method)
  
  data_frame_sub_lod <- rownames_to_column(data.frame(summary(sim_fits[["sub_lod"]][[i]])$coefficients),
                                        var = "covariate") %>%
    mutate(sim=i,
           method="sub_lod",
           df=summary(sim_fits[["sub_lod"]][[i]])$df[2]) %>%
    select(covariate, Estimate, Std..Error, df, sim, method)
  
  data_frame_complete_case <- rownames_to_column(data.frame(summary(sim_fits[["complete_case"]][[i]])$coefficients),
                                           var = "covariate") %>%
    mutate(sim=i,
           method="complete_case",
           df=summary(sim_fits[["complete_case"]][[i]])$df[2]) %>%
    select(covariate, Estimate, Std..Error, df, sim, method)
  
  results_list[[i]] <- data.frame(rbind(data_frame_lodr, data_frame_sub_lod, 
                                        data_frame_complete_case))
}

# Single data frame
results_df <- do.call("rbind", results_list)

# Calc. mean estimate, SE for each covariate for each method
mean_results_df <- results_df %>%
  group_by(covariate, method) %>%
  summarise(mean_est = mean(Estimate),
            mean_se = mean(Std..Error),
            mean_df = mean(df)) %>%
  mutate(method=factor(method, levels=c("lodr", "sub_lod", "complete_case"))) %>%
  arrange(covariate, method)

# Plot 
ggplot(data=mean_results_df, mapping=aes(color=covariate))+
  geom_point(aes(x=method, y=mean_est),position=position_dodge(.7), size=3)+
  geom_errorbar(aes(x=method, ymin=mean_est-mean_se, ymax=mean_est+mean_se),
                position=position_dodge(.7),
                width=0.7)+
  labs(x="Analysis Method", y="Estimate+-SE", 
       title = paste0("Mean estimate and standard errors (SE) of regression analysis\nof simulated data across ", sims, " simulations."),
       subtitle = "Dashed line indicates true parameter values",
       color="Covariate")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme_bw()

# Table
mean_results_df_format <- mean_results_df %>%
  mutate(est_bias = round(mean_est-1,3)) %>%
  select(method, covariate, mean_est, est_bias, mean_se, mean_df) %>%
  arrange(method, covariate)

print(xtable(mean_results_df_format), include.rownames=FALSE)
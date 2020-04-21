########################################################
# ANALYSE OUTPUT FROM SIMULATION, QUESTION 1
########################################################
# DOUBLE TREATMENT FREQUENCY AFTER SKIPPING 1 Year
########################################################
# MODERATE PREVALENCE SETTINGS
########################################################

createdf <- function(foreachResults){
  
  # length(foreachResults)
  # length(foreachResults[[1]])
  # length(foreachResults[[1]]$results)
  # length(foreachResults[[1]]$results[[1]])
  # names(foreachResults[[1]]$results[[1]])
  # head(foreachResults[[1]]$results[[1]]$worms)
  # head(foreachResults[[1]]$results[[2]]$meanIntensitySAC[[1]])
  
  iteration <- c()
  times <- c()
  worms <- c()
  L <- c()
  prevKKSAC <- c()
  prevTrue2 <- c()
  meanIntensitySAC <- c()
  meanIntensityAll <- c()
  prevHISAC <- c()
  prevMHISAC <- c()
  
  outputyears <- c(73, 85, 121)
  
  for(i in 1:(length(foreachResults)-1))
  {
    for(j in outputyears)
    {
      iteration <- c(iteration, i)
      times <- c(times, foreachResults[[i]]$results[[j]]$time)
      prevKKSAC <- c(prevKKSAC, foreachResults[[i]]$results[[j]]$prevKKSAC)
      prevTrue2 <- c(prevTrue2, foreachResults[[i]]$results[[j]]$prevTrue2)
      meanIntensitySAC <- c(meanIntensitySAC, foreachResults[[i]]$results[[j]]$meanIntensitySAC[[1]])
      meanIntensityAll <- c(meanIntensityAll, foreachResults[[i]]$results[[j]]$meanIntensityAll[[1]])
	prevHISAC <- c(prevHISAC, foreachResults[[i]]$results[[j]]$prevHISAC)
	prevMHISAC <- c(prevMHISAC, foreachResults[[i]]$results[[j]]$prevMHISAC)
    }
  }
  
  df <- data.frame(iteration=iteration, time=times, prevKKSAC=prevKKSAC, prevTrue2=prevTrue2,
                   meanIntensitySAC = 24*meanIntensitySAC, meanIntensityAll = 24*meanIntensityAll,
			 prevHISAC=prevHISAC, prevMHISAC=prevMHISAC)

  df$prevMISAC <- df$prevMHISAC - df$prevHISAC
  df$prevLISAC <- df$prevKKSAC - df$prevMHISAC

  return(df)
}


require(dplyr)
dftosummary <- function(df){
  df %>% 
  group_by(time) %>% 
  summarize(KKSAC_prev = mean(prevKKSAC), 
            KKSAC_prev_lower = quantile(prevKKSAC, probs = 0.025), KKSAC_prev_upper = quantile(prevKKSAC, probs = 0.975),
            True_prev = mean(prevTrue2), 
            True_prev_lower = quantile(prevTrue2, probs = 0.025), True_prev_upper = quantile(prevTrue2, probs = 0.975),
            SAC_meanIntensity = mean(meanIntensitySAC), 
            SAC_meanIntensity_lower = quantile(meanIntensitySAC, probs = 0.025), SAC_meanIntensity_upper = quantile(meanIntensitySAC, probs = 0.975),
            All_meanIntensity = mean(meanIntensityAll), 
            All_meanIntensity_lower = quantile(meanIntensityAll, probs = 0.025), All_meanIntensity_upper = quantile(meanIntensityAll, probs = 0.975),
		prevHISAC = mean(prevHISAC),
		prevHISAC_lower = quantile(prevHISAC, probs = 0.025), prevHISAC_upper = quantile(prevHISAC, probs = 0.975),
		prevMHISAC = mean(prevMHISAC),
		prevMHISAC_lower = quantile(prevMHISAC, probs = 0.025), prevMHISAC_upper = quantile(prevMHISAC, probs = 0.975),
		prevMISAC = mean(prevMISAC),
		prevMISAC_lower = quantile(prevMISAC, probs = 0.025), prevMISAC_upper = quantile(prevMISAC, probs = 0.975),
		prevLISAC = mean(prevLISAC),
		prevLISAC_lower = quantile(prevLISAC, probs = 0.025), prevLISAC_upper = quantile(prevLISAC, probs = 0.975)
  )
}


path <- "D:\\STH\\ModellingConsortium\\NTDs-Covid-2019\\NTD-Covid-1-skip-year\\"

df0 <- createdf(get(load(paste0(path, "results_hw_R0_1.5_k0.35_skip0_seed234.RData"))))
df2 <- createdf(get(load(paste0(path, "results_hw_R0_1.5_k0.35_skip2_seed234.RData"))))
df3 <- createdf(get(load(paste0(path, "results_hw_R0_1.5_k0.35_skip3_seed234.RData"))))
df4 <- createdf(get(load(paste0(path, "results_hw_R0_1.5_k0.35_skip4_seed234.RData"))))
df5 <- createdf(get(load(paste0(path, "results_hw_R0_1.5_k0.35_skip5_seed234.RData"))))


df0_summary <-dftosummary(df0)
df2_summary <-dftosummary(df2)
df3_summary <-dftosummary(df3)
df4_summary <-dftosummary(df4)
df5_summary <-dftosummary(df5)

df0_summary <- df0_summary[c(1,3), ]
df2_summary <- df2_summary[2:3, ]
df3_summary <- df3_summary[2:3, ]
df4_summary <- df4_summary[2:3, ]
df5_summary <- df5_summary[2:3, ]

df_summary <- rbind(df0_summary, df2_summary, df3_summary, df4_summary, df5_summary)
Year_skipped <- rep(c("None", 2, 3, 4, 5), each = 2)
Time_label <- rep(c("1yr since stop", "10yrs since start"), 5)
df_summary <- cbind(Year_skipped, Time_label, df_summary)
df_summary$Year_skipped <- relevel(factor(df_summary$Year_skipped), "None")
df_summary$Time_label <- factor(df_summary$Time_label, levels=c("1yr since stop", "10yrs since start"))
#df_summary$time <- factor(df_summary$time, labels = c("1yr since stop", "10yrs since start"))

require(ggplot2)
require(RColorBrewer)
 theme_set(
   theme_minimal() +
     theme(legend.position = "right"))

p1 <- ggplot(data = df_summary, aes(x = Year_skipped , y = KKSAC_prev, color = factor(Year_skipped))) +
        geom_point() +
        scale_y_continuous(limits = c(0, 1)) +
        geom_errorbar(aes(x = Year_skipped, ymin=KKSAC_prev_lower, ymax=KKSAC_prev_upper)) +
        scale_color_brewer(palette = "Set1") +
        facet_grid(cols = vars(Time_label)) +
        theme(legend.position = "none") +
        labs(x = "Year skipped", y = "Kato-Katz prevalence in SAC")  
        theme(axis.title = element_text(size=18), axis.text = element_text(size=18), text = element_text(size=18))
ggsave(file=paste0(path, "ModKKprevSAC.jpeg"), plot=p1, dpi=300)

        
p2 <- ggplot(data = df_summary, aes(x = Year_skipped , y = True_prev, color = factor(Year_skipped))) +
  geom_point() +
  scale_y_continuous(limits = c(0, 1)) +
  geom_errorbar(aes(x = Year_skipped, ymin=True_prev_lower, ymax=True_prev_upper)) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(cols = vars(Time_label)) +
  theme(legend.position = "none") +
  labs(x = "Year skipped", y = "True prevalence in all") +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18), text = element_text(size=18))
ggsave(file=paste0(path, "Modtrueprev.jpeg"), plot=p2, dpi=300)


p3 <- ggplot(data = df_summary, aes(x = Year_skipped , y = SAC_meanIntensity, color = factor(Year_skipped))) +
  geom_point() +
  scale_y_continuous(limits = c(0, 900)) +
  geom_errorbar(aes(x = Year_skipped, ymin=SAC_meanIntensity_lower, ymax=SAC_meanIntensity_upper)) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(cols = vars(Time_label)) +
  theme(legend.position = "none") +
  labs(x = "Year skipped", y = "Mean Intensity in SAC (epg)") + coord_cartesian(ylim=c(0, 400)) +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18), text = element_text(size=18))
ggsave(file=paste0(path, "ModIntensitySAC.jpeg"), plot=p3, dpi=300)


p4 <- ggplot(data = df_summary, aes(x = Year_skipped , y = All_meanIntensity, color = factor(Year_skipped))) +
  geom_point() +
  scale_y_continuous(limits = c(0, 900)) +
  geom_errorbar(aes(x = Year_skipped, ymin=All_meanIntensity_lower, ymax=All_meanIntensity_upper)) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(cols = vars(Time_label)) +
  theme(legend.position = "none") +
  labs(x = "Year skipped", y = "Mean Intensity in all (epg)") +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18), text = element_text(size=18))
ggsave(file=paste0(path, "ModIntensityAll.jpeg"), plot=p4, dpi=300)

require(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)


require(xlsx)
write.xlsx(df_summary, file=paste0(path, "results_summary_q2_freq.xlsx"), sheetName="Moderate intensity prevalence", row.names=FALSE)


########################################################################################################################################
# STACKED BAR PLOTS OF PREVALENCE OF INFECTION INTENSITY CATEGORIES
########################################################################################################################################

library(tidyr)

df_summary_bar <- df_summary[, c("Year_skipped", "Time_label", "prevHISAC", "prevMISAC", "prevLISAC")]
df_summary_bar_long <- gather(df_summary_bar, intensity, value, prevHISAC:prevLISAC, factor_key=TRUE)

labs <- c("Control scenario", "Skip 2nd year", "Skip 3rd year", "Skip 4th year", "Skip 5th year")
levels(df_summary_bar_long$Year_skipped) <- labs
levels(df_summary_bar_long$Time_label) <- c("1yr \nsince stop", "10yrs \nsince start")


theme_set(
   theme_grey() +
     theme(legend.position = "right"))


b1 <- ggplot(data = df_summary_bar_long) + geom_bar(aes(x=Time_label, y=value, fill=intensity), stat="identity", position = position_stack(reverse = TRUE)) 
b1 <- b1 + scale_fill_manual(values=c("red", "orange", "yellow"), name="Intensity of infection", labels=c("Heavy", "Moderate", "Light"))
b1 <- b1 + facet_wrap(. ~ Year_skipped, nrow=2)
b1 <- b1 + ggtitle("20-50% baseline prevalence") + ylab("Prevalence of Kato-Katz hookworm \ninfection in SAC (stacked %)") + xlab("")
b1 <- b1 + theme(plot.title=element_text(size=20), axis.text=element_text(size=18), 
			axis.title.y=element_text(size=18), legend.title=element_text(size=18), legend.text=element_text(size=18),
			strip.text.x=element_text(size=18))
print(b1)
ggsave(file=paste0(path, "barplot_moderate.jpeg"), plot=b1, dpi=300, width=12, height=9)



library(ggplot2)

#################################################################################
# HISTOGRAMS OF WORM BURDENS BY AGE CLASS
#################################################################################

path.as <- "D:\\STH\\ModellingConsortium\\NTDs-Covid-2019\\Ascaris\\tests\\"
path.hw <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Untreated\\"
path.tr <- "D:\\STH\\ModellingConsortium\\NTDs-Covid-2019\\Trichuris\\tests\\"
outpath <- "D:\\STH\\ModellingConsortium\\NTDs-Covid-2019\\ModelComparison\\"

#################################################################################
# FUNCTION TO READ IN RESULTS AND ALLOCATE OUTPUT
#################################################################################

getWormData <- function(path, file, ageGroup=c(5,14))
{
	foreachResults <- get(load(paste0(path, file)))

	worms <- c()
	femaleWorms <- c()

	for(i in 1:(length(foreachResults)-1))
	{
		for(j in 1:length(foreachResults[[i]]$results))
		{
			birthDates <- foreachResults[[i]]$results[[j]]$hosts$birthDate
			ages <- foreachResults[[i]]$results[[j]]$time - birthDates
			hostsInAgeGroup <- which((ages >= ageGroup[1]) & (ages <= ageGroup[2]))

			wormsInAgeGroup <- sum(foreachResults[[i]]$results[[j]]$worms$total[hostsInAgeGroup])
			femaleWormsInAgeGroup <- sum(foreachResults[[i]]$results[[j]]$worms$female[hostsInAgeGroup])
			worms <- c(worms, wormsInAgeGroup)
			femaleWorms <- c(femaleWorms, femaleWormsInAgeGroup)
		}
	}
	return(list(worms=worms, femaleWorms=femaleWorms))
}

getFFWormData <- function(path, file, ageGroup=c(5, 14))
{
	foreachResults <- get(load(paste0(path, file)))

	ffWorms <- c()
	for(i in 1:(length(foreachResults)-1))
	{
		for(j in 1:length(foreachResults[[i]]$results))
		{
			birthDates <- foreachResults[[i]]$results[[j]]$hosts$birthDate
			ages <- foreachResults[[i]]$results[[j]]$time - birthDates
			hostsInAgeGroup <- which((ages >= ageGroup[1]) & (ages <= ageGroup[2]))

			wormsInAgeGroup <- foreachResults[[i]]$results[[j]]$worms[hostsInAgeGroup, ]
			
			fw <- wormsInAgeGroup$female > 0
			mw <- (wormsInAgeGroup$total - wormsInAgeGroup$female) > 0
			ffw.true <- ifelse((fw + mw) == 2, 1, 0)
			ffw <- sum(wormsInAgeGroup$female * ffw.true)
			ffWorms <- c(ffWorms, ffw)
		}
	}
	return(ffWorms)
}


getPrevData <- function(path, file)
{
	foreachResults <- get(load(paste0(path, file)))
	
	iteration <- c()
	times <- c()

	prevKKSAC <- c()
	prevKKAll <- c()
	prevTrue2 <- c()

	meanIntensitySAC <- c()
	meanIntensityAll <- c()

	for(i in 1:(length(foreachResults)-1))
	{
		for(j in 1:length(foreachResults[[i]]$results))
		{
			iteration <- c(iteration, i)
			times <- c(times, foreachResults[[i]]$results[[j]]$time)

			prevKKSAC <- c(prevKKSAC, foreachResults[[i]]$results[[j]]$prevKKSAC)
			prevKKAll <- c(prevKKAll, foreachResults[[i]]$results[[j]]$prevKKAll)
			prevTrue2 <- c(prevTrue2, foreachResults[[i]]$results[[j]]$prevTrue2)

			meanIntensitySAC <- c(meanIntensitySAC, foreachResults[[i]]$results[[j]]$meanIntensitySAC[[1]])
			meanIntensityAll <- c(meanIntensityAll, foreachResults[[i]]$results[[j]]$meanIntensityAll[[1]])
		}
	}
	df <- data.frame(iteration=iteration, times=times, prevKKSAC=prevKKSAC,
				prevKKAll=prevKKAll, prevTrue=prevTrue2, 
				meanIntensitySAC=meanIntensitySAC, meanIntensityAll=meanIntensityAll)
	
	return(df)
}


#################################################################################
# ASCARIS, HIGH PREVALENCE, >50%
#################################################################################

file.as.high <- "results_as_R0_3_k0.8_skip0_seed234.RData"

df.as <- getPrevData(path.as, file.as.high) 

as.worms.SAC <- getWormData(path.as, file.as.high, c(5,14))
as.ffWorms.SAC <- getFFWormData(path.as, file.as.high, c(5,14))

as.worms.adults <- getWormData(path.as, file.as.high, c(15,100))
as.ffWorms.adults <- getFFWormData(path.as, file.as.high, c(15,100))

df.as$total.worms.SAC <- as.worms.SAC$worms
df.as$female.worms.SAC <- as.worms.SAC$femaleWorms
df.as$fert.female.worms.SAC <- as.ffWorms.SAC

df.as$total.worms.adults <- as.worms.adults$worms
df.as$female.worms.adults <- as.worms.adults$femaleWorms
df.as$fert.female.worms.adults <- as.ffWorms.adults



p.as.fem.SAC.high.hist <- ggplot(df.as) + geom_histogram(aes(x=female.worms.SAC, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			labs(subtitle=paste0("Mean: ", round(mean(df.as$female.worms.SAC)), ", 95% CI: ", 
				quantile(df.as$female.worms.SAC, c(0.025, 0.975))[1], ", ",
				quantile(df.as$female.worms.SAC, c(0.025, 0.975))[2]))
print(p.as.fem.SAC.high.hist)
ggsave(paste0(outpath, "hist.as.fem.SAC.high.png"), plot=p.as.fem.SAC.high.hist, dpi=300)


p.as.fert.SAC.high.hist <- ggplot(df.as) + geom_histogram(aes(x=fert.female.worms.SAC, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			labs(subtitle=paste0("Mean: ", round(mean(df.as$fert.female.worms.SAC)), ", 95% CI: ", 
				round(quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[1]), ", ",
				round(quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[2])))
print(p.as.fert.SAC.high.hist)
ggsave(paste0(outpath, "hist.as.fert.SAC.high.png"), plot=p.as.fert.SAC.high.hist, dpi=300)


p.as.fem.adults.high.hist <- ggplot(df.as) + geom_histogram(aes(x=female.worms.adults, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			labs(subtitle=paste0("Mean: ", round(mean(df.as$female.worms.adults)), ", 95% CI: ", 
				quantile(df.as$female.worms.adults, c(0.025, 0.975))[1], ", ",
				quantile(df.as$female.worms.adults, c(0.025, 0.975))[2]))
print(p.as.fem.adults.high.hist)
ggsave(paste0(outpath, "hist.as.fem.adults.high.png"), plot=p.as.fem.adults.high.hist, dpi=300)


p.as.fert.adults.high.hist <- ggplot(df.as) + geom_histogram(aes(x=fert.female.worms.adults, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			labs(subtitle=paste0("Mean: ", round(mean(df.as$fert.female.worms.adults)), ", 95% CI: ", 
				round(quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[1]), ", ",
				round(quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[2])))
print(p.as.fert.adults.high.hist)
ggsave(paste0(outpath, "hist.as.fert.adults.high.png"), plot=p.as.fert.adults.high.hist, dpi=300)



p.as.prevKK.high.hist <- ggplot(df.as) + geom_histogram(aes(x=prevKKSAC, y=(..count..)/sum(..count..)), fill="tan2", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("KK prevalence in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$prevKKSAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$prevKKSAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			labs(subtitle=paste0("Mean: ", signif(mean(df.as$prevKKSAC), digits=3), ", 95% CI: ", 
				signif(quantile(df.as$prevKKSAC, c(0.025, 0.975))[1], digits=3), ", ",
				signif(quantile(df.as$prevKKSAC, c(0.025, 0.975))[2], digits=3)))
print(p.as.prevKK.high.hist)
ggsave(paste0(outpath, "hist.as.KK.high.png"), plot=p.as.prevKK.high.hist, dpi=300)


p.as.prevTrue.high.hist <- ggplot(df.as) + geom_histogram(aes(x=prevTrue, y=(..count..)/sum(..count..)), fill="aquamarine1", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("True prevalence in population") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$prevTrue, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$prevTrue, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			labs(subtitle=paste0("Mean: ", signif(mean(df.as$prevTrue), digits=3), ", 95% CI: ", 
				signif(quantile(df.as$prevTrue, c(0.025, 0.975))[1], digits=3), ", ",
				signif(quantile(df.as$prevTrue, c(0.025, 0.975))[2], digits=3)))
print(p.as.prevTrue.high.hist)
ggsave(paste0(outpath, "hist.as.prevTrue.high.png"), plot=p.as.prevTrue.high.hist, dpi=300)



#################################################################################
# ASCARIS, MODERATE PREVALENCE, 20-50%
#################################################################################

file.as.mod <- "results_as_R0_2.5_k0.2_skip0_seed234.RData"

df.as <- getPrevData(path.as, file.as.mod) 

as.worms.SAC <- getWormData(path.as, file.as.mod, c(5,14))
as.ffWorms.SAC <- getFFWormData(path.as, file.as.mod, c(5,14))

as.worms.adults <- getWormData(path.as, file.as.mod, c(15,100))
as.ffWorms.adults <- getFFWormData(path.as, file.as.mod, c(15,100))

df.as$total.worms.SAC <- as.worms.SAC$worms
df.as$female.worms.SAC <- as.worms.SAC$femaleWorms
df.as$fert.female.worms.SAC <- as.ffWorms.SAC

df.as$total.worms.adults <- as.worms.adults$worms
df.as$female.worms.adults <- as.worms.adults$femaleWorms
df.as$fert.female.worms.adults <- as.ffWorms.adults



p.as.fem.SAC.mod.hist <- ggplot(df.as) + geom_histogram(aes(x=female.worms.SAC, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, 20-50%, k=0.2, R0=2.5") +
			labs(subtitle=paste0("Mean: ", round(mean(df.as$female.worms.SAC)), ", 95% CI: ", 
				quantile(df.as$female.worms.SAC, c(0.025, 0.975))[1], ", ",
				quantile(df.as$female.worms.SAC, c(0.025, 0.975))[2]))
print(p.as.fem.SAC.mod.hist)
ggsave(paste0(outpath, "hist.as.fem.SAC.mod.png"), plot=p.as.fem.SAC.mod.hist, dpi=300)


p.as.fert.SAC.mod.hist <- ggplot(df.as) + geom_histogram(aes(x=fert.female.worms.SAC, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, 20-50%, k=0.2, R0=2.5") +
			labs(subtitle=paste0("Mean: ", round(mean(df.as$fert.female.worms.SAC)), ", 95% CI: ", 
				round(quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[1]), ", ",
				round(quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[2])))
print(p.as.fert.SAC.mod.hist)
ggsave(paste0(outpath, "hist.as.fert.SAC.mod.png"), plot=p.as.fert.SAC.mod.hist, dpi=300)


p.as.fem.adults.mod.hist <- ggplot(df.as) + geom_histogram(aes(x=female.worms.adults, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, 20-50%, k=0.2, R0=2.5") +
			labs(subtitle=paste0("Mean: ", round(mean(df.as$female.worms.adults)), ", 95% CI: ", 
				quantile(df.as$female.worms.adults, c(0.025, 0.975))[1], ", ",
				quantile(df.as$female.worms.adults, c(0.025, 0.975))[2]))
print(p.as.fem.adults.mod.hist)
ggsave(paste0(outpath, "hist.as.fem.adults.mod.png"), plot=p.as.fem.adults.mod.hist, dpi=300)


p.as.fert.adults.mod.hist <- ggplot(df.as) + geom_histogram(aes(x=fert.female.worms.adults, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, 20-50%, k=0.2, R0=2.5") +
			labs(subtitle=paste0("Mean: ", round(mean(df.as$fert.female.worms.adults)), ", 95% CI: ", 
				round(quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[1]), ", ",
				round(quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[2])))
print(p.as.fert.adults.mod.hist)
ggsave(paste0(outpath, "hist.as.fert.adults.mod.png"), plot=p.as.fert.adults.mod.hist, dpi=300)



p.as.prevKK.mod.hist <- ggplot(df.as) + geom_histogram(aes(x=prevKKSAC, y=(..count..)/sum(..count..)), fill="tan2", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("KK prevalence in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$prevKKSAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$prevKKSAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, 20-50%, k=0.2, R0=2.5") +
			labs(subtitle=paste0("Mean: ", signif(mean(df.as$prevKKSAC), digits=3), ", 95% CI: ", 
				signif(quantile(df.as$prevKKSAC, c(0.025, 0.975))[1], digits=3), ", ",
				signif(quantile(df.as$prevKKSAC, c(0.025, 0.975))[2], digits=3)))
print(p.as.prevKK.mod.hist)
ggsave(paste0(outpath, "hist.as.KK.mod.png"), plot=p.as.prevKK.mod.hist, dpi=300)


p.as.prevTrue.mod.hist <- ggplot(df.as) + geom_histogram(aes(x=prevTrue, y=(..count..)/sum(..count..)), fill="aquamarine1", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("True prevalence in population") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.as$prevTrue, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$prevTrue, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, 20-50%, k=0.2, R0=2.5") +
			labs(subtitle=paste0("Mean: ", signif(mean(df.as$prevTrue), digits=3), ", 95% CI: ", 
				signif(quantile(df.as$prevTrue, c(0.025, 0.975))[1], digits=3), ", ",
				signif(quantile(df.as$prevTrue, c(0.025, 0.975))[2], digits=3)))
print(p.as.prevTrue.mod.hist)
ggsave(paste0(outpath, "hist.as.prevTrue.mod.png"), plot=p.as.prevTrue.mod.hist, dpi=300)




#################################################################################
# HOOKWORM, HIGH PREVALENCE, >50%
#################################################################################


getPrevDataHw <- function(path, file)
{
	foreachResults <- get(load(paste0(path, file)))
	
	iteration <- c()
	times <- c()

	prevKK <- c()
	
	for(i in 1:(length(foreachResults)-1))
	{
		for(j in 1:length(foreachResults[[i]]$results))
		{
			iteration <- c(iteration, i)
			times <- c(times, foreachResults[[i]]$results[[j]]$time)

			prevKK <- c(prevKK, foreachResults[[i]]$results[[j]]$prevKK)
		}
	}
	df <- data.frame(iteration=iteration, times=times, prevKKSAC=prevKK)
	
	return(df)
}


#################################################################################


file.hw.high <- "results_hw_R0_9.5_k0.35_seed234.RData"

df.hw <- getPrevDataHw(path.hw, file.hw.high) 

hw.worms.SAC <- getWormData(path.hw, file.hw.high, c(5,14))
hw.ffWorms.SAC <- getFFWormData(path.hw, file.hw.high, c(5,14))

hw.worms.adults <- getWormData(path.hw, file.hw.high, c(15,100))
hw.ffWorms.adults <- getFFWormData(path.hw, file.hw.high, c(15,100))

df.hw$total.worms.SAC <- hw.worms.SAC$worms
df.hw$female.worms.SAC <- hw.worms.SAC$femaleWorms
df.hw$fert.female.worms.SAC <- hw.ffWorms.SAC

df.hw$total.worms.adults <- hw.worms.adults$worms
df.hw$female.worms.adults <- hw.worms.adults$femaleWorms
df.hw$fert.female.worms.adults <- hw.ffWorms.adults



p.hw.fem.SAC.high.hist <- ggplot(df.hw) + geom_histogram(aes(x=female.worms.SAC, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.hw$female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.hw$female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Hookworm, >50%, k=0.35, R0=9.5") +
			labs(subtitle=paste0("Mean:", round(mean(df.hw$female.worms.SAC)), ", 95% CI: ", 
				quantile(df.hw$female.worms.SAC, c(0.025, 0.975))[1], ",",
				round(quantile(df.hw$female.worms.SAC, c(0.025, 0.975))[2])))
print(p.hw.fem.SAC.high.hist)
ggsave(paste0(outpath, "hist.hw.fem.SAC.high.png"), plot=p.hw.fem.SAC.high.hist, dpi=300)


p.hw.fert.SAC.high.hist <- ggplot(df.hw) + geom_histogram(aes(x=fert.female.worms.SAC, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.hw$fert.female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.hw$fert.female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Hookworm, >50%, k=0.35, R0=9.5") +
			labs(subtitle=paste0("Mean:", round(mean(df.hw$fert.female.worms.SAC)), ", 95% CI: ", 
				round(quantile(df.hw$fert.female.worms.SAC, c(0.025, 0.975))[1]), ",",
				round(quantile(df.hw$fert.female.worms.SAC, c(0.025, 0.975))[2])))
print(p.hw.fert.SAC.high.hist)
ggsave(paste0(outpath, "hist.hw.fert.SAC.high.png"), plot=p.hw.fert.SAC.high.hist, dpi=300)


p.hw.fem.adults.high.hist <- ggplot(df.hw) + geom_histogram(aes(x=female.worms.adults, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.hw$female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.hw$female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Hookworm, >50%, k=0.35, R0=9.5") +
			labs(subtitle=paste0("Mean:", round(mean(df.hw$female.worms.adults)), ", 95% CI: ", 
				quantile(df.hw$female.worms.adults, c(0.025, 0.975))[1], ",",
				quantile(df.hw$female.worms.adults, c(0.025, 0.975))[2]))
print(p.hw.fem.adults.high.hist)
ggsave(paste0(outpath, "hist.hw.fem.adults.high.png"), plot=p.hw.fem.adults.high.hist, dpi=300)


p.hw.fert.adults.high.hist <- ggplot(df.hw) + geom_histogram(aes(x=fert.female.worms.adults, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.hw$fert.female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.hw$fert.female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Hookworm, >50%, k=0.35, R0=9.5") +
			labs(subtitle=paste0("Mean:", round(mean(df.hw$fert.female.worms.adults)), ", 95% CI: ", 
				round(quantile(df.hw$fert.female.worms.adults, c(0.025, 0.975))[1]), ",",
				round(quantile(df.hw$fert.female.worms.adults, c(0.025, 0.975))[2])))
print(p.hw.fert.adults.high.hist)
ggsave(paste0(outpath, "hist.hw.fert.adults.high.png"), plot=p.hw.fert.adults.high.hist, dpi=300)



p.hw.prevKK.high.hist <- ggplot(df.hw) + geom_histogram(aes(x=prevKKSAC, y=(..count..)/sum(..count..)), fill="tan2", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("KK prevalence in SAC") + ylab("Frequency") + 
			geom_vline(xintercept=quantile(df.hw$prevKKSAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.hw$prevKKSAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Hookworm, >50%, k=0.35, R0=9.5") +
			labs(subtitle=paste0("Mean: ", signif(mean(df.hw$prevKKSAC), digits=3), ", 95% CI: ", 
				signif(quantile(df.hw$prevKKSAC, c(0.025, 0.975))[1], digits=3), ",",
				signif(quantile(df.hw$prevKKSAC, c(0.025, 0.975))[2], digits=3)))
print(p.hw.prevKK.high.hist)
ggsave(paste0(outpath, "hist.hw.KK.high.png"), plot=p.hw.prevKK.high.hist, dpi=300)



#################################################################################
# HOOKWORM, MODERATE PREVALENCE, 20-50%
#################################################################################

file.hw.mod <- "results_hw_R0_1.5_k0.35_seed234.RData"

df.hw <- getPrevDataHw(path.hw, file.hw.mod) 

hw.worms.SAC <- getWormData(path.hw, file.hw.mod, c(5,14))
hw.ffWorms.SAC <- getFFWormData(path.hw, file.hw.mod, c(5,14))

hw.worms.adults <- getWormData(path.hw, file.hw.mod, c(15,100))
hw.ffWorms.adults <- getFFWormData(path.hw, file.hw.mod, c(15,100))

df.hw$total.worms.SAC <- hw.worms.SAC$worms
df.hw$female.worms.SAC <- hw.worms.SAC$femaleWorms
df.hw$fert.female.worms.SAC <- hw.ffWorms.SAC

df.hw$total.worms.adults <- hw.worms.adults$worms
df.hw$female.worms.adults <- hw.worms.adults$femaleWorms
df.hw$fert.female.worms.adults <- hw.ffWorms.adults



p.hw.fem.SAC.mod.hist <- ggplot(df.hw) + geom_histogram(aes(x=female.worms.SAC, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.hw$female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.hw$female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Hookworm, 20-50%, k=0.35, R0=1.5") +
			labs(subtitle=paste0("Mean:", round(mean(df.hw$female.worms.SAC)), ", 95% CI: ", 
				quantile(df.hw$female.worms.SAC, c(0.025, 0.975))[1], ", ",
				quantile(df.hw$female.worms.SAC, c(0.025, 0.975))[2]))
print(p.hw.fem.SAC.mod.hist)
ggsave(paste0(outpath, "hist.hw.fem.SAC.mod.png"), plot=p.hw.fem.SAC.mod.hist, dpi=300)


p.hw.fert.SAC.mod.hist <- ggplot(df.hw) + geom_histogram(aes(x=fert.female.worms.SAC, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.hw$fert.female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.hw$fert.female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Hookworm, 20-50%, k=0.35, R0=1.5") +
			labs(subtitle=paste0("Mean:", round(mean(df.hw$fert.female.worms.SAC)), ", 95% CI: ", 
				round(quantile(df.hw$fert.female.worms.SAC, c(0.025, 0.975))[1]), ", ",
				round(quantile(df.hw$fert.female.worms.SAC, c(0.025, 0.975))[2])))
print(p.hw.fert.SAC.mod.hist)
ggsave(paste0(outpath, "hist.hw.fert.SAC.mod.png"), plot=p.hw.fert.SAC.mod.hist, dpi=300)


p.hw.fem.adults.mod.hist <- ggplot(df.hw) + geom_histogram(aes(x=female.worms.adults, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.hw$female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.hw$female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Hookworm, 20-50%, k=0.35, R0=1.5") +
			labs(subtitle=paste0("Mean:", round(mean(df.hw$female.worms.adults)), ", 95% CI: ", 
				quantile(df.hw$female.worms.adults, c(0.025, 0.975))[1], ", ",
				quantile(df.hw$female.worms.adults, c(0.025, 0.975))[2]))
print(p.hw.fem.adults.mod.hist)
ggsave(paste0(outpath, "hist.hw.fem.adults.mod.png"), plot=p.hw.fem.adults.mod.hist, dpi=300)


p.hw.fert.adults.mod.hist <- ggplot(df.hw) + geom_histogram(aes(x=fert.female.worms.adults, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.hw$fert.female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.hw$fert.female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Hookworm, 20-50%, k=0.35, R0=1.5") +
			labs(subtitle=paste0("Mean:", round(mean(df.hw$fert.female.worms.adults)), ", 95% CI: ", 
				round(quantile(df.hw$fert.female.worms.adults, c(0.025, 0.975))[1]), ", ",
				round(quantile(df.hw$fert.female.worms.adults, c(0.025, 0.975))[2])))
print(p.hw.fert.adults.mod.hist)
ggsave(paste0(outpath, "hist.hw.fert.adults.mod.png"), plot=p.hw.fert.adults.mod.hist, dpi=300)



p.hw.prevKK.mod.hist <- ggplot(df.hw) + geom_histogram(aes(x=prevKKSAC, y=(..count..)/sum(..count..)), fill="tan2", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("KK prevalence in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.hw$prevKKSAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.hw$prevKKSAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Hookworm, 20-50%, k=0.35, R0=1.5") +
			labs(subtitle=paste0("Mean:", signif(mean(df.hw$prevKKSAC), digits=3), ", 95% CI: ", 
				signif(quantile(df.hw$prevKKSAC, c(0.025, 0.975))[1], digits=3), ", ",
				signif(quantile(df.hw$prevKKSAC, c(0.025, 0.975))[2], digits=3)))
print(p.hw.prevKK.mod.hist)
ggsave(paste0(outpath, "hist.hw.KK.mod.png"), plot=p.hw.prevKK.mod.hist, dpi=300)





#################################################################################
# TRICHURIS, HIGH PREVALENCE, >50%
#################################################################################


file.tr.high <- "results_tr_R0_1.8_k0.38_skip0_seed234.RData"

df.tr <- getPrevData(path.tr, file.tr.high) 

tr.worms.SAC <- getWormData(path.tr, file.tr.high, c(5,14))
tr.ffWorms.SAC <- getFFWormData(path.tr, file.tr.high, c(5,14))

tr.worms.adults <- getWormData(path.tr, file.tr.high, c(15,100))
tr.ffWorms.adults <- getFFWormData(path.tr, file.tr.high, c(15,100))

df.tr$total.worms.SAC <- tr.worms.SAC$worms
df.tr$female.worms.SAC <- tr.worms.SAC$femaleWorms
df.tr$fert.female.worms.SAC <- tr.ffWorms.SAC

df.tr$total.worms.adults <- tr.worms.adults$worms
df.tr$female.worms.adults <- tr.worms.adults$femaleWorms
df.tr$fert.female.worms.adults <- tr.ffWorms.adults




p.tr.fem.SAC.high.hist <- ggplot(df.tr) + geom_histogram(aes(x=female.worms.SAC, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Trichuris, >50%, k=0.38, R0=1.2") +
			labs(subtitle=paste0("Mean: ", round(mean(df.tr$female.worms.SAC)), ", 95% CI: ", 
				quantile(df.tr$female.worms.SAC, c(0.025, 0.975))[1], ", ",
				quantile(df.tr$female.worms.SAC, c(0.025, 0.975))[2]))
print(p.tr.fem.SAC.high.hist)
ggsave(paste0(outpath, "hist.tr.fem.SAC.high.png"), plot=p.tr.fem.SAC.high.hist, dpi=300)


p.tr.fert.SAC.high.hist <- ggplot(df.tr) + geom_histogram(aes(x=fert.female.worms.SAC, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$fert.female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$fert.female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Trichuris, >50%, k=0.38, R0=1.2") +
			labs(subtitle=paste0("Mean: ", round(mean(df.tr$fert.female.worms.SAC)), ", 95% CI: ", 
				round(quantile(df.tr$fert.female.worms.SAC, c(0.025, 0.975))[1]), ", ",
				round(quantile(df.tr$fert.female.worms.SAC, c(0.025, 0.975))[2])))
print(p.tr.fert.SAC.high.hist)
ggsave(paste0(outpath, "hist.tr.fert.SAC.high.png"), plot=p.tr.fert.SAC.high.hist, dpi=300)


p.tr.fem.adults.high.hist <- ggplot(df.tr) + geom_histogram(aes(x=female.worms.adults, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Trichuris, >50%, k=0.38, R0=1.2") +
			labs(subtitle=paste0("Mean: ", round(mean(df.tr$female.worms.adults)), ", 95% CI: ", 
				quantile(df.tr$female.worms.adults, c(0.025, 0.975))[1], ", ",
				quantile(df.tr$female.worms.adults, c(0.025, 0.975))[2]))
print(p.tr.fem.adults.high.hist)
ggsave(paste0(outpath, "hist.tr.fem.adults.high.png"), plot=p.tr.fem.adults.high.hist, dpi=300)


p.tr.fert.adults.high.hist <- ggplot(df.tr) + geom_histogram(aes(x=fert.female.worms.adults, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$fert.female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$fert.female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			labs(subtitle=paste0("Mean: ", round(mean(df.tr$fert.female.worms.adults)), ", 95% CI: ", 
				round(quantile(df.tr$fert.female.worms.adults, c(0.025, 0.975))[1]), ", ",
				round(quantile(df.tr$fert.female.worms.adults, c(0.025, 0.975))[2])))
print(p.tr.fert.adults.high.hist)
ggsave(paste0(outpath, "hist.tr.fert.adults.high.png"), plot=p.tr.fert.adults.high.hist, dpi=300)



p.tr.prevKK.high.hist <- ggplot(df.tr) + geom_histogram(aes(x=prevKKSAC, y=(..count..)/sum(..count..)), fill="tan2", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("KK prevalence in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$prevKKSAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$prevKKSAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			labs(subtitle=paste0("Mean: ", signif(mean(df.tr$prevKKSAC), digits=3), ", 95% CI: ", 
				signif(quantile(df.tr$prevKKSAC, c(0.025, 0.975))[1], digits=3), ", ",
				signif(quantile(df.tr$prevKKSAC, c(0.025, 0.975))[2], digits=3)))
print(p.tr.prevKK.high.hist)
ggsave(paste0(outpath, "hist.tr.KK.high.png"), plot=p.tr.prevKK.high.hist, dpi=300)


p.tr.prevTrue.high.hist <- ggplot(df.tr) + geom_histogram(aes(x=prevTrue, y=(..count..)/sum(..count..)), fill="aquamarine1", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("True prevalence in population") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$prevTrue, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$prevTrue, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			labs(subtitle=paste0("Mean: ", signif(mean(df.tr$prevTrue), digits=3), ", 95% CI: ", 
				signif(quantile(df.tr$prevTrue, c(0.025, 0.975))[1], digits=3), ", ",
				signif(quantile(df.tr$prevTrue, c(0.025, 0.975))[2], digits=3)))
print(p.tr.prevTrue.high.hist)
ggsave(paste0(outpath, "hist.tr.prevTrue.high.png"), plot=p.tr.prevTrue.high.hist, dpi=300)



#################################################################################
# ASCARIS, MODERATE PREVALENCE, 20-50%
#################################################################################

file.tr.mod <- "results_tr_R0_1.5_k0.12_skip0_seed234.RData"

df.tr <- getPrevData(path.tr, file.tr.mod) 

tr.worms.SAC <- getWormData(path.tr, file.tr.mod, c(5,14))
tr.ffWorms.SAC <- getFFWormData(path.tr, file.tr.mod, c(5,14))

tr.worms.adults <- getWormData(path.tr, file.tr.mod, c(15,100))
tr.ffWorms.adults <- getFFWormData(path.tr, file.tr.mod, c(15,100))

df.tr$total.worms.SAC <- tr.worms.SAC$worms
df.tr$female.worms.SAC <- tr.worms.SAC$femaleWorms
df.tr$fert.female.worms.SAC <- tr.ffWorms.SAC

df.tr$total.worms.adults <- tr.worms.adults$worms
df.tr$female.worms.adults <- tr.worms.adults$femaleWorms
df.tr$fert.female.worms.adults <- tr.ffWorms.adults



p.tr.fem.SAC.mod.hist <- ggplot(df.tr) + geom_histogram(aes(x=female.worms.SAC, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Trichuris, 20-50%, k=0.12, R0=1.5") +
			labs(subtitle=paste0("Mean: ", round(mean(df.tr$female.worms.SAC)), ", 95% CI: ", 
				quantile(df.tr$female.worms.SAC, c(0.025, 0.975))[1], ", ",
				quantile(df.tr$female.worms.SAC, c(0.025, 0.975))[2]))
print(p.tr.fem.SAC.mod.hist)
ggsave(paste0(outpath, "hist.tr.fem.SAC.mod.png"), plot=p.tr.fem.SAC.mod.hist, dpi=300)


p.tr.fert.SAC.mod.hist <- ggplot(df.tr) + geom_histogram(aes(x=fert.female.worms.SAC, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$fert.female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$fert.female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Trichuris, 20-50%, k=0.12, R0=1.5") +
			labs(subtitle=paste0("Mean: ", round(mean(df.tr$fert.female.worms.SAC)), ", 95% CI: ", 
				round(quantile(df.tr$fert.female.worms.SAC, c(0.025, 0.975))[1]), ", ",
				round(quantile(df.tr$fert.female.worms.SAC, c(0.025, 0.975))[2])))
print(p.tr.fert.SAC.mod.hist)
ggsave(paste0(outpath, "hist.tr.fert.SAC.mod.png"), plot=p.tr.fert.SAC.mod.hist, dpi=300)


p.tr.fem.adults.mod.hist <- ggplot(df.tr) + geom_histogram(aes(x=female.worms.adults, y=(..count..)/sum(..count..)), fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Trichuris, 20-50%, k=0.12, R0=1.5") +
			labs(subtitle=paste0("Mean: ", round(mean(df.tr$female.worms.adults)), ", 95% CI: ", 
				quantile(df.tr$female.worms.adults, c(0.025, 0.975))[1], ", ",
				quantile(df.tr$female.worms.adults, c(0.025, 0.975))[2]))
print(p.tr.fem.adults.mod.hist)
ggsave(paste0(outpath, "hist.tr.fem.adults.mod.png"), plot=p.tr.fem.adults.mod.hist, dpi=300)


p.tr.fert.adults.mod.hist <- ggplot(df.tr) + geom_histogram(aes(x=fert.female.worms.adults, y=(..count..)/sum(..count..)), fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in adults") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$fert.female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$fert.female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Trichuris, 20-50%, k=0.12, R0=1.5") +
			labs(subtitle=paste0("Mean: ", round(mean(df.tr$fert.female.worms.adults)), ", 95% CI: ", 
				round(quantile(df.tr$fert.female.worms.adults, c(0.025, 0.975))[1]), ", ",
				round(quantile(df.tr$fert.female.worms.adults, c(0.025, 0.975))[2])))
print(p.tr.fert.adults.mod.hist)
ggsave(paste0(outpath, "hist.tr.fert.adults.mod.png"), plot=p.tr.fert.adults.mod.hist, dpi=300)



p.tr.prevKK.mod.hist <- ggplot(df.tr) + geom_histogram(aes(x=prevKKSAC, y=(..count..)/sum(..count..)), fill="tan2", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("KK prevalence in SAC") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$prevKKSAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$prevKKSAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Trichuris, 20-50%, k=0.12, R0=1.5") +
			labs(subtitle=paste0("Mean: ", signif(mean(df.tr$prevKKSAC), digits=3), ", 95% CI: ", 
				signif(quantile(df.tr$prevKKSAC, c(0.025, 0.975))[1], digits=3), ", ",
				signif(quantile(df.tr$prevKKSAC, c(0.025, 0.975))[2], digits=3)))
print(p.tr.prevKK.mod.hist)
ggsave(paste0(outpath, "hist.tr.KK.mod.png"), plot=p.tr.prevKK.mod.hist, dpi=300)


p.tr.prevTrue.mod.hist <- ggplot(df.tr) + geom_histogram(aes(x=prevTrue, y=(..count..)/sum(..count..)), fill="aquamarine1", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("True prevalence in population") + ylab("Frequency") +
			geom_vline(xintercept=quantile(df.tr$prevTrue, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.tr$prevTrue, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Trichuris, 20-50%, k=0.12, R0=1.5") +
			labs(subtitle=paste0("Mean: ", signif(mean(df.tr$prevTrue), digits=3), ", 95% CI: ", 
				signif(quantile(df.tr$prevTrue, c(0.025, 0.975))[1], digits=3), ", ",
				signif(quantile(df.tr$prevTrue, c(0.025, 0.975))[2], digits=3)))
print(p.tr.prevTrue.mod.hist)
ggsave(paste0(outpath, "hist.tr.prevTrue.mod.png"), plot=p.tr.prevTrue.mod.hist, dpi=300)








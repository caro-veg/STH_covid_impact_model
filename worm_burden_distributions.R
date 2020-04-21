library(ggplot2)

#################################################################################
# HISTOGRAMS OF WORM BURDENS BY AGE CLASS
#################################################################################

path.as <- "D:\\STH\\ModellingConsortium\\NTDs-Covid-2019\\Ascaris\\tests\\"
path.hw <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Untreated\\correctedTimeStep\\"
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



p.as.fem.SAC.high.hist <- ggplot(df.as, aes(x=female.worms.SAC)) + geom_histogram(fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in SAC") +
			geom_vline(xintercept=quantile(df.as$female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=2300, y=900, size=6,
				label=paste("Mean:", round(mean(df.as$female.worms.SAC)), "\n95% CI:\n", 
				quantile(df.as$female.worms.SAC, c(0.025, 0.975))[1], ",",
				quantile(df.as$female.worms.SAC, c(0.025, 0.975))[2]))
print(p.as.fem.SAC.high.hist)
ggsave(paste0(outpath, "hist.as.fem.SAC.high.png"), plot=p.as.fem.SAC.high.hist, dpi=300)


p.as.fert.SAC.high.hist <- ggplot(df.as, aes(x=fert.female.worms.SAC)) + geom_histogram(fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in SAC") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=2290, y=900, size=6,
				label=paste("Mean:", round(mean(df.as$fert.female.worms.SAC)), "\n95% CI:\n", 
				round(quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[1]), ",",
				round(quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[2])))
print(p.as.fert.SAC.high.hist)
ggsave(paste0(outpath, "hist.as.fert.SAC.high.png"), plot=p.as.fert.SAC.high.hist, dpi=300)


p.as.fem.adults.high.hist <- ggplot(df.as, aes(x=female.worms.adults)) + geom_histogram(fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in adults") +
			geom_vline(xintercept=quantile(df.as$female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=1740, y=900, size=6,
				label=paste("Mean:", round(mean(df.as$female.worms.adults)), "\n95% CI:\n", 
				quantile(df.as$female.worms.adults, c(0.025, 0.975))[1], ",",
				quantile(df.as$female.worms.adults, c(0.025, 0.975))[2]))
print(p.as.fem.adults.high.hist)
ggsave(paste0(outpath, "hist.as.fem.adults.high.png"), plot=p.as.fem.adults.high.hist, dpi=300)


p.as.fert.adults.high.hist <- ggplot(df.as, aes(x=fert.female.worms.adults)) + geom_histogram(fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in adults") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=1740, y=900, size=6,
				label=paste("Mean:", round(mean(df.as$fert.female.worms.adults)), "\n95% CI:\n", 
				round(quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[1]), ",",
				round(quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[2])))
print(p.as.fert.adults.high.hist)
ggsave(paste0(outpath, "hist.as.fert.adults.high.png"), plot=p.as.fert.adults.high.hist, dpi=300)



p.as.prevKK.high.hist <- ggplot(df.as, aes(x=prevKKSAC)) + geom_histogram(fill="tan2", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("KK prevalence in SAC") +
			geom_vline(xintercept=quantile(df.as$prevKKSAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$prevKKSAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=0.68, y=900, size=6,
				label=paste("Mean:", signif(mean(df.as$prevKKSAC), digits=3), "\n95% CI:\n", 
				signif(quantile(df.as$prevKKSAC, c(0.025, 0.975))[1], digits=3), ",",
				signif(quantile(df.as$prevKKSAC, c(0.025, 0.975))[2], digits=3)))
print(p.as.prevKK.high.hist)
ggsave(paste0(outpath, "hist.as.KK.high.png"), plot=p.as.prevKK.high.hist, dpi=300)


p.as.prevTrue.high.hist <- ggplot(df.as, aes(x=prevTrue)) + geom_histogram(fill="aquamarine1", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("True prevalence in population") +
			geom_vline(xintercept=quantile(df.as$prevTrue, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$prevTrue, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=0.795, y=900, size=6,
				label=paste("Mean:", signif(mean(df.as$prevTrue), digits=3), "\n95% CI:\n", 
				signif(quantile(df.as$prevTrue, c(0.025, 0.975))[1], digits=3), ",",
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



p.as.fem.SAC.mod.hist <- ggplot(df.as, aes(x=female.worms.SAC)) + geom_histogram(fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in SAC") +
			geom_vline(xintercept=quantile(df.as$female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=1000, y=6000, size=6,
				label=paste("Mean:", round(mean(df.as$female.worms.SAC)), "\n95% CI:\n", 
				quantile(df.as$female.worms.SAC, c(0.025, 0.975))[1], ",",
				quantile(df.as$female.worms.SAC, c(0.025, 0.975))[2]))
print(p.as.fem.SAC.mod.hist)
ggsave(paste0(outpath, "hist.as.fem.SAC.mod.png"), plot=p.as.fem.SAC.mod.hist, dpi=300)


p.as.fert.SAC.mod.hist <- ggplot(df.as, aes(x=fert.female.worms.SAC)) + geom_histogram(fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in SAC") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=1000, y=5000, size=6,
				label=paste("Mean:", round(mean(df.as$fert.female.worms.SAC)), "\n95% CI:\n", 
				round(quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[1]), ",",
				round(quantile(df.as$fert.female.worms.SAC, c(0.025, 0.975))[2])))
print(p.as.fert.SAC.mod.hist)
ggsave(paste0(outpath, "hist.as.fert.SAC.mod.png"), plot=p.as.fert.SAC.mod.hist, dpi=300)


p.as.fem.adults.mod.hist <- ggplot(df.as, aes(x=female.worms.adults)) + geom_histogram(fill="mistyrose", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total female worms in adults") +
			geom_vline(xintercept=quantile(df.as$female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=800, y=6000, size=6,
				label=paste("Mean:", round(mean(df.as$female.worms.adults)), "\n95% CI:\n", 
				quantile(df.as$female.worms.adults, c(0.025, 0.975))[1], ",",
				quantile(df.as$female.worms.adults, c(0.025, 0.975))[2]))
print(p.as.fem.adults.mod.hist)
ggsave(paste0(outpath, "hist.as.fem.adults.mod.png"), plot=p.as.fem.adults.mod.hist, dpi=300)


p.as.fert.adults.mod.hist <- ggplot(df.as, aes(x=fert.female.worms.adults)) + geom_histogram(fill="lavender", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("Total fertilised female worms in adults") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=750, y=8000, size=6,
				label=paste("Mean:", round(mean(df.as$fert.female.worms.adults)), "\n95% CI:\n", 
				round(quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[1]), ",",
				round(quantile(df.as$fert.female.worms.adults, c(0.025, 0.975))[2])))
print(p.as.fert.adults.mod.hist)
ggsave(paste0(outpath, "hist.as.fert.adults.mod.png"), plot=p.as.fert.adults.mod.hist, dpi=300)



p.as.prevKK.mod.hist <- ggplot(df.as, aes(x=prevKKSAC)) + geom_histogram(fill="tan2", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("KK prevalence in SAC") +
			geom_vline(xintercept=quantile(df.as$prevKKSAC, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$prevKKSAC, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=0.18, y=6000, size=6,
				label=paste("Mean:", signif(mean(df.as$prevKKSAC), digits=3), "\n95% CI:\n", 
				signif(quantile(df.as$prevKKSAC, c(0.025, 0.975))[1], digits=3), ",",
				signif(quantile(df.as$prevKKSAC, c(0.025, 0.975))[2], digits=3)))
print(p.as.prevKK.mod.hist)
ggsave(paste0(outpath, "hist.as.KK.mod.png"), plot=p.as.prevKK.mod.hist, dpi=300)


p.as.prevTrue.mod.hist <- ggplot(df.as, aes(x=prevTrue)) + geom_histogram(fill="aquamarine1", colour="black", alpha=0.6) + 
			theme_classic(base_size=16, base_line_size = 16/(22 * 16/11)) +
			xlab("True prevalence in population") +
			geom_vline(xintercept=quantile(df.as$prevTrue, c(0.025, 0.975))[1], colour="red") +
			geom_vline(xintercept=quantile(df.as$prevTrue, c(0.025, 0.975))[2], colour="red") +
			ggtitle("Ascaris, >50%, k=0.8, R0=3") +
			geom_text(x=0.25, y=8000, size=6,
				label=paste("Mean:", signif(mean(df.as$prevTrue), digits=3), "\n95% CI:\n", 
				signif(quantile(df.as$prevTrue, c(0.025, 0.975))[1], digits=3), ",",
				signif(quantile(df.as$prevTrue, c(0.025, 0.975))[2], digits=3)))
print(p.as.prevTrue.mod.hist)
ggsave(paste0(outpath, "hist.as.prevTrue.mod.png"), plot=p.as.prevTrue.mod.hist, dpi=300)






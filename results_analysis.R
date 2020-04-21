########################################################
# ANALYSE OUTPUT FROM SIMULATION
########################################################

path <- "D:\\STH\\ModellingConsortium\\NTDs-Covid-2019\\NTD-Covid-1-skip-year\\"
filename <- "results_hw_R0_1.5_k0.35_skip_seed234.RData"

foreachResults <- get(load(paste0(path, filename)))


length(foreachResults)
length(foreachResults[[1]])
length(foreachResults[[1]]$results)
length(foreachResults[[1]]$results[[1]])
names(foreachResults[[1]]$results[[1]])
head(foreachResults[[1]]$results[[1]]$worms)

scenario <- c()
iteration <- c()
times <- c()

prevKKSAC <- c()
prevKKAll <- c()
prevTrue2 <- c()

meanIntensitySAC <- c()
meanIntensityAll <- c()

outputIndices <- c(73, 85, 121)

for(i in 1:(length(foreachResults)-1))
{
	for(j in outputIndices)
	{
		scenario <- c("baseline")
		iteration <- c(iteration, i)
		times <- c(times, foreachResults[[i]]$results[[j]]$time)

		prevKKSAC <- c(prevKKSAC, foreachResults[[i]]$results[[j]]$prevKKSAC)
		prevKKAll <- c(prevKKAll, foreachResults[[i]]$results[[j]]$prevKKAll)
		prevTrue2 <- c(prevTrue2, foreachResults[[i]]$results[[j]]$prevTrue2)

		meanIntensitySAC <- c(meanIntensitySAC, foreachResults[[i]]$results[[j]]$meanIntensitySAC[[1]])
		meanIntensityAll <- c(meanIntensityAll, foreachResults[[i]]$results[[j]]$meanIntensityAll[[1]])
	}
}


df <- data.frame(scenario <- scenario,
			iteration=iteration, 
			time=times, 
			prevKKSAC=prevKKSAC, 
			prevKKAll=prevKKAll, 
			prevTrue2=prevTrue2,
			meanIntensitySAC=meanIntensitySAC,
			meanIntensityAll=meanIntensityAll)







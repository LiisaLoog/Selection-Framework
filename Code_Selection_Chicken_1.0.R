################################################################################################################ 
# Liisa Loog 2017
# Code for analysing selection using ancient DNA data.

# Reference:
# Loog et al. (2017) Inferring allele frequency trajectories from ancient DNA indicates 
# that selection on a chicken gene coincided with changes in medieval husbandry practices. 
# Molecular Biology and Evolution. doi:10.1093/molbev/msx142.

################################################################################################################ 
################################################ Load Libraries ################################################ 
################################################################################################################ 

library(deSolve)

################################################################################################################
################################################ Define function for ODE solver ################################
################################################################################################################ 

func <- function(t,y, params){
	s <- params[1]
	m <- params[2]
	f.ex <- params[3]
	t.sel.start <- params[4]
	t.gf.start <- params[5]
	t.gf.stop <- params[6]

	if( t < t.sel.start){ 
		s<- 0 
	}  
	if(t< t.gf.start){ 
		m <- 0
	}
	if(t > t.gf.stop){ 
		m <- 0
	}
	return(list(s*y*y*(1-y)+ m*(f.ex-y))) # selection on a recessive allele + gene flow 
}

################################################################################################################
################################################ Load data ##################################################### 
################################################################################################################

data <- read.csv("TSHR_times_genotypes.csv", head = T, sep = ",")
#data <- read.csv("BCDO2_times_genotypes.csv", head = T, sep = ",")
nsamples <- length(data[,1])

################################################################################################################
################################################ Define allele frequencies #####################################
################################################################################################################

freq.D <- as.numeric(data[,4]) ### frequency of ancestral allele
freq.A <- as.numeric(data[,5]) ### frequency of derived allele

################################################################################################################
################################################ Define sample ages ############################################ 
################################################################################################################

sample.t.lower <- as.numeric(data[,6]) ### Lower age bundary (in years AD)
sample.t.upper <- as.numeric(data[,7]) ### Upper age bundary (in years AD)
sample.t.mean <- as.numeric(data[,8]) ### Mean age (in years AD) if not using age ranges

################################################################################################################
################################################ Define fixed parameters ####################################### 
################################################################################################################

#time.period <- c(min(sample.t.lower),max(sample.t.upper)) # Time period to consider (in years AD)
time.period <- c(10,2000) # Time period to consider (in years AD)
time.step <- 10
gf.period <- c(1750, 1990) # Time period for gene flow from Asia (in years AD)
total.gf <- 0.15/(gf.period[2] - gf.period[1] ) # Fraction of influx/ gene flow per time unit (year)
f.external <- 0.99 # Frequncy of the derived allele in the Asian chiken population.
ranges <- T # F or T - Whether to use age ranges (T) or point estimates for ages (F)
t.sim.start <- min(sample.t.lower)


################################################################################################################
#################################### Define ranges for free parameters #########################################
################################################################################################################

s.values <- 10^seq(from = -4, to = -0.0, by = 0.040) # Selection coefficients (on a log scale)
sel.start <- seq(time.period[1],time.period[2],time.step) # Selection starting times (in years AD)
f.anc <- array(0:100)/100 # Ancestral allele frequencies for the derived allele 

################################################################################################################
############ Set up an array to record the maximum liklihood allele frequecy trajectory ########################
################################################################################################################
 
curve.sample.times <- seq(time.period[1], time.period[2], length.out = 200)
ML.curve <- matrix(-Inf, nrow =100 , ncol = length(curve.sample.times))
max.log.L.cube <- -Inf

################################################################################################################
################## Set up a 3D array to record likelihoods for all parameter combinations #####################
################################################################################################################

results.cube <-array(NA,dim=c(length(s.values),length(sel.start),length(f.anc)))

################################################################################################################
################################################ Do a parameter sweep ##########################################
################################################################################################################

	start.time <- Sys.time()

for(i in 1:length(s.values)){
	print(paste("i = ", i))
	for(j in 1:length(sel.start)){
		print(paste("j = ", j))
		for(k in 1:length(f.anc)){
#			print(paste("k = ", k))

			selection.strength <- s.values[i]
			t.sel.start <- sel.start[j]
			f.ancestral <- f.anc[k]
			t.gf.start <- gf.period[1]
			t.gf.stop <- gf.period[2]
			influx <- total.gf
			f.ex <- f.external
			params <- c(selection.strength, influx, f.ex, t.sel.start, t.gf.start, t.gf.stop)
			log.L <- 0

			if (ranges == F) {
				f.t <- lsoda(y = f.ancestral, times= sample.t.mean , func =func, parms =params,  rtol = 1e-8, atol = 1e-8)
				#f.t <- f.t[2:length(f.t[,1]),]
				# Cap the likelhoods to the range [1e-8,1-1e-8] to prevent NaNs in the likelihood calculation
				f.t[f.t > 1-1E-8] = 1-1E-8
				f.t[f.t < 1E-8] = 1E-8 
				log.La <- freq.D*log(f.t[,2])+freq.A*log(1-f.t[,2])
				log.La[is.na(log.La )] <- 0
				log.L <- sum(log.La)
				results.cube[i,j,k] <- log.L
			} else {
				# Calculate average data likelihood for each sample over nsamples uniformly
				# spaced points in the date range for the sample.			
				for(m in seq(1: nsamples)){
					#print(paste("m = ", m))
					# generate sampling points in the date range for the sample
					sample.times = sample.t.lower[m] + (sample.t.upper[m] - sample.t.lower[m]) *  (seq(from = 1, to =10, by = 1)-0.5)/10	
					# add starting time of the curve for the ODE solver
					sample.times <- c(t.sim.start, sample.times) 
					# generate the allele frequency curve
					f.t <- lsoda(y = f.ancestral, times=sample.times, func =func, parms =params,  rtol = 1e-8, atol = 1e-8) 

					f.t <- f.t[2:length(f.t[,1]),] # Remove frequency from initial point
					# Cap the likelhoods to the range [1e-8,1-1e-8] to prevent NaNs in the likelihood calculation
					f.t[f.t > 1-1E-8] = 1-1E-8
					f.t[f.t < 1E-8] = 1E-8 
					log.La <- freq.D[m]*log(f.t[,2])+freq.A[m]*log(1-f.t[,2]) # calculate the likelihood of the allele frequency curve
					# Simply adding likelihoods cause numerical problems due to underflow.
					# I solve this by scaling the likelihoods before addition and reversing the scaling for the result.					
					max.log.La <- max(log.La)
					log.La <- log(mean(exp(log.La-max.log.La))) + max.log.La
					log.L <- log.L+log.La 
				}
			results.cube[i,j,k] <- log.L
			}

################################################################################################################
############################ Record liklihood-weighted allele frequecy trajectory ##############################
################################################################################################################

			max.log.L.cube <- max(max.log.L.cube,  log.L)
			f.t.curve <- lsoda(y=f.ancestral, times=curve.sample.times, func=func, parms=params,  rtol=1e-8, atol=1e-8)
			# Sample allele frequencies uniformly in time
			for(i.curve.time in 1:length(ML.curve[1,])){
				index <- min(floor(f.t.curve[i.curve.time,2]*100) +1,100)
				ML.curve[index,i.curve.time] <- log(exp(ML.curve[index,i.curve.time] - max.log.L.cube) +  exp(log.L - max.log.L.cube)) + max.log.L.cube 
			}	
		}
	}
}
stop.time <- Sys.time()
print(stop.time - start.time)

sum(is.na(results.cube))
#save.image(file = "chicken.TSHR.Rdata") / save.image(file = "chicken.BCDO2.Rdata") # Save R.image

################################################################################################################
############################# Calculate & plot posterior marginal likelihoods ##################################
################################################################################################################

############################# Marginal likelihoods for selection coefficients ##############################

new.results.cube <- exp(results.cube-max(results.cube)) 
ML.s <- apply(new.results.cube,1,mean)
plot(log10(s.values), ML.s , type = "l", main ="marginal likelihood of selection strength", xlab = "log10 of selection strength", ylab= "marginal likelihood")

############################# Marginal likelihoods for starting times of selection ##########################

new.results.cube <- exp(results.cube-max(results.cube))
ML.sel.start <- apply(new.results.cube,2,mean)
plot(sel.start, ML.sel.start, type = "l",  main ="marginal likelihood of starting time of selection", xlab = "starting time of selection", ylab= "marginal likelihood")

############################# Marginal likelihoods for ancestral allele frequencies ##########################

new.results.cube <- exp(results.cube-max(results.cube))
ML.f.anc <- apply(new.results.cube,3,mean)
plot(f.anc, ML.f.anc, type = "l", main = "marginal likelihood of ancestral allele frequency", xlab = "ancestral (selected) allele frequency", ylab= "marginal likelihood")

################################################################################################################
############################ Normalize & plot the allele frequecy trajectory ###################################
################################################################################################################

ML.curve <-  exp(ML.curve-max(ML.curve)) 
ML.curve[ML.curve > 0.2] <- 0.2 # Set a cut-off  
image(t(ML.curve), xaxt="n",  yaxt="n")
axis(1, at = seq(0,1, length.out=11), labels = pretty(curve.sample.times, n = 11))
axis(2, at = seq(0,1, length.out=11), labels = 0:10/10)

################################################################################################################
############################################## THE END #########################################################
################################################################################################################
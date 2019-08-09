library(adaptivetau)
library(countrycode)
attach(mtcars)

# add country codes to prevalence data
prevalences <- read.csv("/Users/christinecho/Documents/Annual-data/PEI/IHME-GBD_2017_DATA-385b67b1-1/IHME-GBD_2017_DATA-385b67b1-1.csv")
location_code <- countrycode(prevalences$location_name, origin = "country.name", destination = "iso3c")
prevalences <- cbind(prevalences, location_code)
allcountries <- as.character(aggregate_mean_predictions_df$country)
par(mfrow = c(4, 2))
par(mar = c(1,2,1,1))

# run simulations for all countries
for(i in 1:length(allcountries)) {
  # skip if population data or prevalence data not available
  if (!(allcountries[i] %in% pop_data$iso3c) | !(allcountries[i] %in% prevalences$location_code)) next
  
  # initial values
  N = subset(pop_data, iso3c == allcountries[i])$value
  W = as.integer(subset(prevalences, location_code == allcountries[i])$val)
  X = 0
  S = as.integer(N-W-X)
  epsilon = 0.05
  rt = 1/7
  rx = 1/13
  rw = 1/14
  beta = N/S*(epsilon*rt+(1-epsilon)*rw)
  rnovo = 1^-6
  
  init.values = c(
    S = S,   # susceptible humans
    W = W,   # infected wild type humans
    X = X)   # infected resistant humans
  
  # transitions
  transitions = list(c(S = +1, W = -1), # infected wild type to susceptible
                     c(S = -1, W = +1), # susceptible to infected wild type
                     c(W = -1, X = +1), # infected wild type to infected resistant
                     c(S = +1, X = -1), # infected resistant to susceptible
                     c(S = -1, X = +1)) # susceptible to infected resistant
  
  # function to calculate rates given variables and params. 
  rates <- function(x, params, t) {
    return(c(x["W"]*(params$epsilon*params$rt + (1 - params$epsilon)*params$rw), 
             params$beta*x["S"]*x["W"]/params$N,
             params$epsilon*x["W"]*params$rnovo,
             x["X"]*params$rx,
             params$beta*x["S"]*x["X"]/params$N))
  }
  
  params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = epsilon, rnovo = rnovo, N=N)
  
  # run simulation
  r = ssa.adaptivetau(init.values, transitions, rates, params, tf=100)
  
  # plot
  matplot(r[,"time"], r[,c("S","W", "X")], type='l', xlab='Time',
          ylab='Counts', ylim = c(0, max(r[, "X"], r[, "W"])))
  legend("right", legend=c("S", "W", "X"), lty=1:3, col=1:3)
  title(allcountries[i], line = -2)
}

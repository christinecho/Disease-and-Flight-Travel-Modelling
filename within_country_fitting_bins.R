library(adaptivetau)
library(EasyABC)
attach(mtcars)
library(countrycode)
library(dplyr)
library(DescTools)

# get all countries and don't include points without date
mcr1_df <- read.csv("/Users/christinecho/Documents/Annual-data/PEI/mcr1_edited.csv", stringsAsFactors = FALSE)
#country_codes <- countrycode(mcr1_df$Country, origin = "country.name", destination = "iso3c")
#mcr1_df <- cbind(mcr1_df, country_codes)

# order countries by ascending GDP per capita to get prevalence values
gdp_per_capita <- read.csv("/Users/christinecho/Documents/Annual-data/gdpdata.csv", stringsAsFactors = FALSE)
countrycodes <- countrycode(gdp_per_capita$country, origin = "country.name", destination = "iso3c")
gdp_per_capita <- cbind.data.frame(gdp_per_capita, countrycodes)
mcr1_df <- subset(mcr1_df, !is.na(Year))
allcountries <- unique(mcr1_df$Country.Codes)
prevalence_min <- 1e-4
prevalence_max <- 1e-3
par(mfrow = c(2, 2))
par(mar = c(2,2,1,1))

for(i in 19:length(allcountries)) {
  if (is.na(i)) next
  # compute prevalence, larger prevalence if gdp is lower
  curr_country_code <- as.character(allcountries[i])
  print(curr_country_code)
  if ((!curr_country_code %in% pop_data$iso3c) | (!curr_country_code %in% gdp_per_capita$countrycodes)) next
  N = subset(pop_data, iso3c == curr_country_code)$value
  country_rank <- gdp_per_capita[which(gdp_per_capita$countrycodes == curr_country_code), ]$rank
  country_prevalence_val <- prevalence_min + (prevalence_max - prevalence_min) * (nrow(gdp_per_capita)-country_rank)/nrow(gdp_per_capita)
  
  # bin experimental values by year
  country_outbreaks <- mcr1_df[which(mcr1_df$Country.Codes == curr_country_code), ]
  # will be useful when separating by month
  #country_outbreaks$Formatted_Dates <- as.Date(country_outbreaks$Formatted_Dates, format = c("%m/%d/%y"))
  # country_oubreaks$Year <- as.numeric(country_outbreaks$Year)
  year_freq <- data.frame(table(country_outbreaks$Year))
  time_span <- max(as.numeric(year_freq$Var1))- min(as.numeric(year_freq$Var1)) + 1
  if(time_span == 1) next
  # other initial values
  W = as.integer(country_prevalence_val * N)
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
  
  params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = epsilon, rnovo = rnovo, N=N);
  
  # run simulation
  r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
  
  # plot experimental values
  matplot(r[,"time"], r[,c("S","W", "X")], type='l', xlab='Time',
          ylab='Counts', col = c("black", "red", "blue"), ylim = c(0, max(r[,"W"], r[,"X"])))
  legend("right", legend=c("S", "W", "X"), lty=1:3, col=c("black", "red", "blue"))
  matpoints(as.numeric(year_freq$Var1), year_freq$Freq, pch = 4)
  print(time_span)
  
  # plot initial simulation values
  for (i in 0:time_span) {
    matpoints(i, r[which(abs(r[,"time"]-i)==min(abs(r[, "time"]-i))), ]["X"], pch = 2)
  }
  
  # function to run for abc sequential
  # returns Euclidean distance from simulated and experimental values
  runSimulation <- function(new_epsilon) {
    params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = new_epsilon, rnovo = rnovo, N=N);
    r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
    simulated_values <- c()
    for (i in 0:(time_span-1)) {
      simulated_values <- c(simulated_values, r[which(abs(r[,"time"]-i)==min(abs(r[, "time"]-i))), ]["X"])
    }
    return(calculateDistance(simulated_values, year_freq$Freq))
  }
  
  # returns sum of Euclidean distance between experimental and simulated values
  calculateDistance <- function(x, y){
    dist = 0
    for(i in 1:length(x)) {
      dist = dist + abs(x[i] - y[i])
    }
    return(dist)
  }
  
  # ABC model
  tolerance = c(160, 80, 40, 20, 10, 5, 2.5, 1)/50
  epsilon_sample = list(c("unif", 0, 0.2))
  ABC_Beaumont <- ABC_sequential(method="Beaumont", model=runSimulation, prior=epsilon_sample, 
                                 nb_simul=100, summary_stat_target = (0), tolerance_tab = tolerance, verbose = TRUE)
  
  # plot lines converging onto actual line
  for (i in 1:length(tolerance)) {
    new_epsilon = median(ABC_Beaumont[["intermediary"]][[i]][["posterior"]][,2])
    print(new_epsilon)
    params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = new_epsilon, rnovo = rnovo, N=N);
    r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
    a<- list()
    for (j in 0:time_span) {
      a[j+1] <- r[which(abs(r[,"time"]-j)==min(abs(r[, "time"]-j))), ]["X"]
    }
    matlines(0:time_span, a, col = i+2)
  }
  title(curr_country_code, line = -2)
  text(time_span/2, max(r[,"W"])/2, new_epsilon)
}
library(adaptivetau)
library(EasyABC)
attach(mtcars)
library(countrycode)
library(dplyr)
library(DescTools)
library(RColorBrewer)

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
par(mar = c(1,2,1,1))

for(i in 1:length(allcountries)) {
  if (is.na(i)) next
  # compute prevalence, larger prevalence if gdp is lower
  curr_country_code <- as.character(allcountries[i])
  print(curr_country_code)
  if ((!curr_country_code %in% pop_data$iso3c) | (!curr_country_code %in% gdp_per_capita$countrycodes)) next
  N = subset(pop_data, iso3c == curr_country_code)$value
  country_rank <- gdp_per_capita[which(gdp_per_capita$countrycodes == curr_country_code), ]$rank
  country_prevalence_val <- prevalence_min + (prevalence_max - prevalence_min) * country_rank/nrow(gdp_per_capita)
  
  # bin experimental values by year
  country_outbreaks <- mcr1_df[which(mcr1_df$Country.Codes == curr_country_code), ]
  first_outbreak = min(country_outbreaks$Year)
  years <- seq(min(country_outbreaks$Year), max(country_outbreaks$Year))
  years <- data.frame(years)
  year_frequencies <- as.data.frame(table(country_outbreaks$Year))
  year_frequencies$Var1 = as.numeric(levels(year_frequencies$Var1))[year_frequencies$Var1]
  year_frequencies <- merge(years, year_frequencies, by.x = "years", by.y = "Var1", all=T)
  year_frequencies[is.na(year_frequencies)] <- 0
  time_span <- max(as.numeric(year_frequencies$years))- min(as.numeric(year_frequencies$years)) + 1
  if(time_span == 1) next
  
  # other initial values
  W = as.integer(country_prevalence_val * N)
  X = 0
  S = as.integer(N-W-X)
  epsilon = 0.05
  prop = 1
  rt = 1/7 * prop
  rx = 1/13 * prop
  rw = 1/14 * prop
  beta = N/S*(epsilon*rt+(1-epsilon)*rw)
  rnovo = 1e-6
  set.seed(3)
  init.values = c(
    S = S,   # susceptible humans
    W = W,   # infected wild type humans
    X = X)   # infected resistant humans
  
  params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = epsilon, rnovo = rnovo, N=N)
  
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
  
  # run simulation
  r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
  
  # plot experimental values
  #matplot(r[,"time"], r[,c("S","W", "X")], type='l', xlab='Time',
          #ylab='Counts', col = c("black", "red", "blue"), ylim = c(0, 100))
  #legend("right", legend=c("S", "W", "X"), lty=1:3, col=c("black", "red", "blue"))
  matplot(as.numeric(year_frequencies$years - first_outbreak), year_frequencies$Freq, pch = 4, type = 'l', ylim = c(0,100))
  
  # plot initial simulation values
  #for (i in 0:time_span) {
    #matpoints(i, r[which(abs(r[,"time"]-i)==min(abs(r[, "time"]-i))), ]["X"], pch = 2)
  #}
  
  # function to run for abc sequential
  # returns Euclidean distance from simulated and experimental values
  runSimulation <- function(x) {
    new_epsilon = x[1]
    rproportion = x[2]
    rt = 1/7*rproportion
    rx = 1/13*rproportion
    rw = 1/14*rproportion
    N = subset(pop_data, iso3c == curr_country_code)$value
    beta = N/S*(new_epsilon*rt+(1-new_epsilon)*rw)
    params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = new_epsilon, rnovo = rnovo, N=N);
    
    r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
    simulated_values <- c()
    mostRecent = 0
    for (i in 0:(time_span-1)) {
      closest_time <- r[which(abs(r[,"time"]-i)==min(abs(r[, "time"]-i))), ]["X"]
      if(is.na(closest_time)) {
        simulated_values <- c(simulated_values, mostRecent) #might need to fix this
      } else {
        simulated_values <- c(simulated_values, closest_time)
        mostRecent = closest_time
      }
    }
    return(calculateDistance(simulated_values, year_frequencies$Freq))
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
  tolerance = c(160, 80, 40, 20, 10, 5, 2.5, 1)/100000
  epsilon_sample = c("unif", 0, 1)
  rproportion = c("unif", 0, 1000)
  prior = list(epsilon_sample, rproportion)
  ABC_Beaumont <- ABC_sequential(method="Beaumont", model=runSimulation, prior=prior, 
                                 nb_simul=1000, summary_stat_target=0, tolerance_tab = tolerance, verbose = T)
  
  colors = colorRampPalette(brewer.pal(8,"Reds"))(length(tolerance))
  # plot lines converging onto actual line
  for (k in 1:length(tolerance)) {
    rej <- abc((0), ABC_Beaumont[["intermediary"]][[k]][["posterior"]][,2:3], ABC_Beaumont[["intermediary"]][[k]][["posterior"]][,4], tol = 0.001, method = "rejection")
    print(rej$ss)
    new_epsilon = rej$unadj.values[1]#median(rej$unadj.values[,1])
    new_proportion = rej$unadj.values[2]#median(rej$unadj.values[,2])
    rt = 1/7 * new_proportion
    rw = 1/13 * new_proportion
    rx = 1/14 * new_proportion
    N = subset(pop_data, iso3c == curr_country_code)$value
    beta = N/S*(new_epsilon*rt+(1-new_epsilon)*rw)
    params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = new_epsilon, rnovo = rnovo, N=N);
    r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
    
    #print(r[which(abs(r[,"time"]-k)==min(abs(r[, "time"]-k))), ]["X"])
    matlines(r[,"time"], r[, "X"], col = colors[k])
    #a<- list()
    #for (j in 0:time_span) {
      #a[j+1] <- r[which(abs(r[,"time"]-j)==min(abs(r[, "time"]-j))), ]["X"]
    #}
    #matlines(0:time_span, a, col = colors[k])
  }
  title(curr_country_code, line = -2)
  text(4, 9, new_epsilon)
  text(4, 5, new_proportion)
}
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
mcr1_df$Date <- as.integer(mcr1_df$Date)
mcr1_df <- subset(mcr1_df, !is.na(Date))
allcountries <- unique(mcr1_df$Country.Codes)
prevalence_min <- 1e-4
prevalence_max <- 1e-3
par(mfrow = c(2, 2))
par(mar = c(1,2,1,1))

# save all epsilon values to data frame
epsilon_within_df <- data.frame(Country=character(), Epsilon_Value = numeric())

for(i in 1:length(allcountries)) {
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
  first_outbreak = min(country_outbreaks$Date)
  days <- seq(min(country_outbreaks$Date), max(country_outbreaks$Date))
  days <- data.frame(days)
  day_frequencies <- as.data.frame(table(country_outbreaks$Date))
  day_frequencies$Var1 = as.numeric(levels(day_frequencies$Var1))[day_frequencies$Var1]
  day_frequencies <- merge(days, day_frequencies, by.x = "days", by.y = "Var1", all=T)
  day_frequencies[is.na(day_frequencies)] <- 0
  
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
          ylab='Counts', col = c("black", "red", "blue"), ylim = c(0,max(r[,"W"],r[,"X"])))
  legend("right", legend=c("S", "W", "X"), lty=1:3, col=c("black", "red", "blue"))
  matpoints(day_frequencies$days - first_outbreak, day_frequencies$Freq, pch = 4, cex = 0.1)
  time_span <- max(as.numeric(day_frequencies$days))- min(as.numeric(day_frequencies$days)) + 1
  
  # function to run for abc sequential
  # returns Euclidean distance from simulated and experimental values
  runSimulation <- function(new_epsilon) {
    
    params = list(beta = N/S*(new_epsilon*rt+(1-new_epsilon)*rw), rt = rt, rx = rx, rw = rw, epsilon = new_epsilon, rnovo = rnovo, N=N);
    e = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
    simulated_values <- c()
    for (s in 0:(time_span-1)) {
      simulated_values <- c(simulated_values, e[which(abs(e[,"time"]-s)==min(abs(e[, "time"]-s))), ]["X"])
    }
    return(calculateDistance(simulated_values, day_frequencies$Freq))
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
  tolerance = c(160, 80, 40, 20, 10, 5, 2.5, 1)
  epsilon_sample = list(c("unif", 0, 0.01))
  ABC_Beaumont <- ABC_sequential(method="Beaumont", model=runSimulation, prior=epsilon_sample, 
                                 nb_simul=100, summary_stat_target = (0), tolerance_tab = tolerance, verbose = TRUE)
  
  colors = colorRampPalette(brewer.pal(8,"Reds"))(length(tolerance))
  # plot lines converging onto actual line
  for (k in 1:length(tolerance)) {
    rej <- abc(c(0), ABC_Beaumont[["intermediary"]][[k]][["posterior"]][,2], ABC_Beaumont[["intermediary"]][[k]][["posterior"]][,3], tol = 0.2, method = "rejection")
    new_epsilon = median(rej$unadj.values)
    beta = N/S*(new_epsilon*rt+(1-new_epsilon)*rw)
    params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = new_epsilon, rnovo = rnovo, N=N);
    r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
    #print(r[which(abs(r[,"time"]-k)==min(abs(r[, "time"]-k))), ]["X"])
    matlines(r[,"time"], r[,"X"], col = colors[k])
  }
  epsilon_within_df[i, ] = c(curr_country_code, new_epsilon)
  title(curr_country_code, line = -2)
  text(4, 9, new_epsilon)
}
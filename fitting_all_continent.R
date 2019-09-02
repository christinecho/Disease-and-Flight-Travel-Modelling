library(adaptivetau)
library(EasyABC)
attach(mtcars)
library(countrycode)
library(dplyr)
library(DescTools)
library(RColorBrewer)

mcr1_df <- read.csv("/Users/christinecho/Documents/Annual-data/PEI/mcr1_edited.csv", stringsAsFactors = FALSE)

# load in gdp per capita for prevalence values
gdp_per_capita <- read.csv("/Users/christinecho/Documents/Annual-data/gdpdata.csv", stringsAsFactors = FALSE)
countrycodes <- countrycode(gdp_per_capita$country, origin = "country.name", destination = "iso3c")
gdp_per_capita <- cbind.data.frame(gdp_per_capita, countrycodes)

# clean data set, removing data points without human hosts or dates
mcr1_df <- subset(mcr1_df, !is.na(Year))
mcr1_df <- subset(mcr1_df, Host == "Human") 
mcr1_df <- subset(mcr1_df, Country.Codes %in% gdp_per_capita$countrycodes)

allcountries <- unique(mcr1_df$Country.Codes)
prevalence_min <- 1e-4
prevalence_max <- 1e-3

# add continents column to datasets
mcr1_df$continent = countrycode(mcr1_df$Country, origin = "country.name", destination = "continent")
gdp_per_capita$continent = countrycode(gdp_per_capita$country, origin = "country.name", destination = "continent")

# separate by West vs rest of the world
numGroups = 2
asiaEtcContinents = c("Asia", "Africa")
westContinents = c("Americas", "Europe")

mcr1AsiaEtc = subset(mcr1_df, continent %in% asiaEtcContinents)
mcr1West = subset(mcr1_df, continent %in% westContinents)

# compute average prevalence for two groups
gdpAsiaEtc = subset(gdp_per_capita, continent %in% asiaEtcContinents)
asiaEtcPrev = prevalence_min + (prevalence_max - prevalence_min) * mean(gdpAsiaEtc$rank)/nrow(gdp_per_capita)
gdpWest = subset(gdp_per_capita, continent %in% westContinents)
westPrev = prevalence_min + (prevalence_max - prevalence_min) * mean(gdpWest$rank)/nrow(gdp_per_capita)

# get total populations
Na = subset(pop_data, country %in% mcr1AsiaEtc$Country)
Nw = subset(pop_data, country %in% mcr1West$Country)
  
# bin experimental values by year
firstOutbreakAsia = min(mcr1AsiaEtc$Year)
yearsAsia <- seq(firstOutbreakAsia, max(mcr1AsiaEtc$Year))
yearsAsia <- data.frame(yearsAsia)
year_frequenciesAsia <- as.data.frame(table(mcr1AsiaEtc$Year))
year_frequenciesAsia$Var1 = as.numeric(levels(year_frequenciesAsia$Var1))[year_frequenciesAsia$Var1]
year_frequenciesAsia <- merge(yearsAsia, year_frequenciesAsia, by.x = "yearsAsia", by.y = "Var1", all=T)
year_frequenciesAsia[is.na(year_frequenciesAsia)] <- 0
# remove last two years, might have been some intervention
year_frequenciesAsia = subset(year_frequenciesAsia, years <= 2016)
time_spanAsia <- max(as.numeric(year_frequenciesAsia$years))- min(as.numeric(year_frequenciesAsia$years)) + 1

firstOutbreakWest = min(mcr1West$Year)
yearsWest <- seq(firstOutbreakWest, max(mcr1West$Year))
yearsWest <- data.frame(yearsWest)
year_frequenciesWest <- as.data.frame(table(mcr1West$Year))
year_frequenciesWest$Var1 = as.numeric(levels(year_frequenciesWest$Var1))[year_frequenciesWest$Var1]
year_frequenciesWest <- merge(yearsWest, year_frequenciesWest, by.x = "yearsWest", by.y = "Var1", all=T)
year_frequenciesWest[is.na(year_frequenciesWest)] <- 0
# remove last two years, might have been some intervention
year_frequenciesAsia = subset(year_frequenciesAsia, years <= 2016)
time_spanWest <- max(as.numeric(year_frequenciesAsia$years))- min(as.numeric(year_frequenciesAsia$years)) + 1
  
# initial values
Wa = as.integer(asiaEtcPrev * Na)
Xa = 0
Sa = as.integer(Na-Wa-Xa)

Ww = as.integer(westPrev * Nw)
Xw = 0
Sw = Nw-Ww-Xw

init.values = c(
   S = Sa,   # susceptible humans
   W = Wa,   # infected wild type humans
   X = Xa)   # infected resistant humans

rnovo = 1e-6
tolerance = c(160, 80, 40, 20, 10, 5, 2.5, 1)/100000
epsilon_sample = c("unif", 0, 1)
rproportion = c("unif", 0, 1000)
prior = list(epsilon_sample, rproportion)

# run simulation for Asia and rest of world
ABC_Asia <- ABC_sequential(method="Beaumont", model=runSimulation, prior=prior, 
                               nb_simul=1000, summary_stat_target=0, tolerance_tab = tolerance, verbose = T)
# run simulation for West

rnovo = 1e-5

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
    W = as.integer(country_prevalence_val * N)
    X = 0
    S = as.integer(N-W-X)
    
    init.values = c(
      S = S,   # susceptible humans
      W = W,   # infected wild type humans
      X = X)   # infected resistant humans
    
    beta = N/S*(new_epsilon*rt+(1-new_epsilon)*rw)
    
    params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = new_epsilon, rnovo = rnovo, N=N);
    
    #set.seed(1)
    r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
    #set.seed(NULL)
    simulated_values <- c()
    
    # add points, and account for if a year is the average of two years with data points
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
    if(calculateDistance(simulated_values, year_frequencies$Freq) < 200) {
      #print(calculateDistance(simulated_values, year_frequencies$Freq))
      #matlines(r[,"time"], r[, "X"], col = "pink")
      #print(new_epsilon)
      #print(rproportion)
      #print(tail(r))
      #print(simulated_values)
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
                                 nb_simul=1000, summary_stat_target=0, tolerance_tab = tolerance, verbose = T, use_seed = F)
  
  colors = colorRampPalette(brewer.pal(8,"Reds"))(length(tolerance))
  # plot lines converging onto actual line
  for(k in 1:length(tolerance)) {
    rej <- abc((0), ABC_Beaumont[["intermediary"]][[k]][["posterior"]][,2:3], ABC_Beaumont[["intermediary"]][[k]][["posterior"]][,4], tol = 0.001, method = "rejection")
    eps = rej$unadj.values[1]#median(rej$unadj.values[,1])
    new_proportion = rej$unadj.values[2]#median(rej$unadj.values[,2])
    print(rej$ss)
    print(eps)
    print(new_proportion)
    rt = 1/7 * new_proportion
    rx = 1/13 * new_proportion
    rw = 1/14 * new_proportion
    N = subset(pop_data, iso3c == curr_country_code)$value
    W = as.integer(country_prevalence_val * N)
    X = 0
    S = as.integer(N-W-X)
    init.values = c(
      S = S,   # susceptible humans
      W = W,   # infected wild type humans
      X = X)   # infected resistant humans
    
    beta = N/S*(eps*rt+(1-eps)*rw)
    
    new_params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = eps, rnovo = rnovo, N=N);
    #set.seed(1)
    a = ssa.adaptivetau(init.values, transitions, rates, new_params, tf=time_span)
    
    View(a)
    #print(r[which(abs(r[,"time"]-k)==min(abs(r[, "time"]-k))), ]["X"])
    matlines(a[,"time"], a[, "X"], col = colors[k])
    #a<- list()
    #for (j in 0:time_span) {
    #a[j+1] <- r[which(abs(r[,"time"]-j)==min(abs(r[, "time"]-j))), ]["X"]
    #}
    #matlines(0:time_span, a, col = colors[k])
  }
  title(curr_country_code, line = -2)
  text(4, 9, new_epsilon)
  
}

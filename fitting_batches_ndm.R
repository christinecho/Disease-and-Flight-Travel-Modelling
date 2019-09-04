library(adaptivetau)
library(EasyABC)
attach(mtcars)
library(countrycode)
library(dplyr)
library(zoo)
library(DescTools)
library(RColorBrewer)


# load in gdp per capita for prevalence values
ndm1_df <- read.csv("/Users/christinecho/Documents/Annual-data/PEI/ndm1.csv", stringsAsFactors = FALSE)
gdp_per_capita <- read.csv("/Users/christinecho/Documents/Annual-data/gdpdata.csv", stringsAsFactors = FALSE)
countrycodes <- countrycode(gdp_per_capita$country, origin = "country.name", destination = "iso3c")
gdp_per_capita <- cbind.data.frame(gdp_per_capita, countrycodes)
par(mfrow = c(2, 2))
par(mar = c(2,2,1,1))

# clean data set, removing data points without human hosts or dates
ndm1_df$date = as.Date(ndm1_df$date)
ndm1_df <- subset(ndm1_df, !is.na(date))
ndm1_df <- subset(ndm1_df, Country %in% gdp_per_capita$country)

prevalence_min <- 1e-4
prevalence_max <- 1e-3

# add continents column to datasets
gdp_per_capita$continent = countrycode(gdp_per_capita$country, origin = "country.name", destination = "continent")
pop_data$continent = countrycode(pop_data$country, origin = "country.name", destination = "continent")

# separate by West vs rest of the world
numGroups = 2
asiaEtcContinents = c("Asia", "Africa","Australia")
westContinents = c("America", "Europe")

groups = list(asiaEtcContinents, westContinents)

for (i in 1:2) {
  # compute average prevalence
  continents = groups[[i]]
  ndmOutbreaks = subset(ndm1_df, Continent %in% continents)
  gdp = subset(gdp_per_capita, continent %in% continents)
  prev = prevalence_min + (prevalence_max - prevalence_min) * mean(gdp$rank)/nrow(gdp_per_capita)
  
  # bin experimental values by month
  # note: some pretty crappy code here but does the job, probably need to rewrite this 
  month_frequencies = ndmOutbreaks %>%
  mutate(month = format(date, "%m"), year = format(date, "%Y")) %>%
  group_by(month, year) %>%
  summarise(total = n()) %>%
  arrange(year, month)  
  subset(month_frequencies, !is.na(month) & !is.na(year))
  month_frequencies$Date = as.Date(as.yearmon(paste(month_frequencies$year, month_frequencies$month), "%Y %m"))
  months = seq(min(month_frequencies$Date), max(month_frequencies$Date), by="months")
  month_frequencies = merge(data.frame(months), month_frequencies, by.x="months",by.y="Date",all=T)
  month_frequencies[is.na(month_frequencies)] <- 0
  
  # sparse data, so going to try to fit data to cumulative sums
  cumulativeSums = cumsum(month_frequencies$total)
  
  time_span <- nrow(month_frequencies)-1
  
  # initial values
  rnovo = 1e-6
  rt = 1/7*30
  rx = 1/13*30
  rw = 1/14*30
  
  N = sum(subset(pop_data, continent %in% continents)$value)
  W = as.integer(prev * N)
  X = 0
  S = N-W-X
  
  init.values = c(
    S = S,
    W = W,
    X = X
  )
  print(N)
  print(prev)
  matplot(c(0:(time_span)), cumulativeSums, pch = 4, type = 'l', ylim = c(0,200))
  
  # returns sum of squares between experimental and simulated values
  calculateDistance <- function(x, y){
    dist = 0
    for(i in 1:length(x)) {
      dist = dist + abs(x[i] - y[i])**2
    }
    return(dist)
  }
  
  # function to run for abc sequential
  runSimulation <- function(x) {
    epsilon = x[1]
    rproportion = x[2]
    N = sum(subset(pop_data, continent %in% continents)$value)
    rt = 1/7*rproportion
    rx = 1/13*rproportion
    rw = 1/14*rproportion
    beta = N/S*(epsilon*rt+(1-epsilon)*rw)
    params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = epsilon, rnovo = rnovo, N=N)
    r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
    set.seed(NULL)
    simulated_values <- c()
    
    # smooth spline can only work with four unique values
    if (length(unique(r[,"X"])) <= 10) {
      return(calculateDistance(rep(0,time_span+1), cumulativeSums))
    }
    else {
    print("printing")
    #print(head(r))
    #print(tail(r))
    ss = smooth.spline(x=r[,"time"], y=r[,"X"], all.knots = T, spar = 1.4)
    prediction = stats::predict(ss, x = 0:time_span)
   # print(prediction$y)
    simulated_values = cumsum(prediction$y)
    #print(simulated_values)
    #matplot(r[,"time"], r[,"X"], pch = 4, type = 'l', ylim = c(0,200), col="blue")
    #matlines(c(0:(time_span)), prediction$y, pch = 4, type = 'l', ylim = c(0,200), col=colors[k])
    return(calculateDistance(simulated_values, cumulativeSums))
    }
  }
  
  # fit simulations and fit epsilon and rate proportion
  tolerance = c(160, 80, 40, 20, 10, 5, 2.5, 1)/10000
  epsilon_sample = c("unif", 0, 0.05)
  rproportion = c("unif", 0, 100)
  prior = list(epsilon_sample, rproportion)
  #prior = list(epsilon_sample)
  ABC <- ABC_sequential(method="Beaumont", model=runSimulation, prior=prior, 
                        nb_simul=1000, summary_stat_target=0, tolerance_tab = tolerance, verbose = T)
  
  # plot lines converging onto actual line
  colors = colorRampPalette(brewer.pal(8,"Reds"))(length(tolerance))
  for(k in 1:length(tolerance)) {
    rej <- abc((0), ABC[["intermediary"]][[k]][["posterior"]][,2:3], ABC[["intermediary"]][[k]][["posterior"]][,4], tol = 0.01, method = "rejection")
    #rej <- abc((0), ABC[["intermediary"]][[k]][["posterior"]][,2], ABC[["intermediary"]][[k]][["posterior"]][,3], tol = 0.001, method = "rejection")
    eps = median(rej$unadj.values[1])
    new_proportion = median(rej$unadj.values[2])#median(rej$unadj.values[,2])
    print(mean(rej$ss))
    print(eps)
    print(new_proportion)
    rt = 1/7 * new_proportion
    rx = 1/13 * new_proportion
    rw = 1/14 * new_proportion
    N = sum(subset(pop_data, continent %in% continents)$value)
    beta = N/S*(eps*rt+(1-eps)*rw)
    params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = eps, rnovo = rnovo, N=N)
    
    # run simulations two times in a row with same initial parameters to see if they're similar
    a = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
    print("once")
    print(tail(a))
    a = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
    set.seed(NULL)
    print("twice")
    print(tail(a))
    ss = smooth.spline(x=a[,"time"], y=a[,"X"], all.knots = T, spar = 1.4)
    prediction = stats::predict(ss, x = 0:time_span)
    simulated_values = cumsum(prediction$y)
    print(simulated_values)
    matlines(c(0:(time_span)), simulated_values, pch = 4, type = 'l', ylim = c(0,200), col=colors[k])
  }  
}

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


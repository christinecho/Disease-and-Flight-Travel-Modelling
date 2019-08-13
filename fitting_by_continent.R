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

# group by continent
mcr1_df %>%
mutate(Continent = countrycode(Country, 'country.name', 'continent')) %>%
group_by(Continent)

for(i in 1:length(unique(mcr1_df$Continent))) {
  continent <- unique(mcr1_df$Continent)[i]
  if (is.na(i)) next
  # compute prevalence, larger prevalence if gdp is lower
  
  # bin experimental values by year
  country_outbreaks <- mcr1_df[which(mcr1_df$Continent == continent), ]
  country_prevalence_val <- prevalence_min
  
  print(continent)
  print(nrow(country_outbreaks))
  
  #country_outbreaks$Formatted_Dates <- as.Date(country_outbreaks$Formatted_Dates, format = c("%m/%d/%y"))
  # country_oubreaks$Year <- as.numeric(country_outbreaks$Year)
  time_span <- max(country_outbreaks$Year) - min(country_outbreaks$Year) + 1
  Years <- seq(min(country_outbreaks$Year), max(country_outbreaks$Year))
  Years <- data.frame(Years)
  year_frequencies <- as.data.frame(table(country_outbreaks$Year))
  year_frequencies$Var1 = as.numeric(levels(year_frequencies$Var1)[year_frequencies$Var1])
  year_frequencies <- merge(Years, year_frequencies, by.x = "Years", by.y = "Var1", all=T)
  year_frequencies[is.na(year_frequencies)] <- 0
  days <- (year_frequencies$Years - min(year_frequencies$Years))*365
  
  if(time_span == 1) next
  # other initial values
  allPop <- pop_data[pop_data$iso3c %in% country_outbreaks$Country.Codes, ]
  N = sum(allPop$value)
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
  r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span*365)
  
  # plot experimental values
  matplot(r[,"time"], r[,c("S","W", "X")], type='l', xlab='Time',
          ylab='Counts', col = c("black", "red", "blue"), ylim = c(0, max(r[,"W"], r[,"X"])))
  legend("right", legend=c("S", "W", "X"), lty=1:3, col=c("black", "red", "blue"))
  matpoints(days, year_frequencies$Freq, pch = 4)

  # plot initial simulation values
  for (i in seq(0, (time_span-1)*365, 365)) {
    matpoints(i, r[which(abs(r[,"time"]-i)==min(abs(r[, "time"]-i))), ]["X"], pch = 2)
  }
  
  # function to run for abc sequential
  # returns Euclidean distance from simulated and experimental values
  runSimulation <- function(new_epsilon) {
    #print(new_epsilon)
    params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = new_epsilon, rnovo = rnovo, N=N);
    r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span*365)
    simulated_values <- c()
    
    for (i in seq(0, (time_span-1)*365, 365)) {
      simulated_values <- c(simulated_values, r[which(abs(r[,"time"]-i)==min(abs(r[, "time"]-i))), ]["X"])
    }
    return(calculateDistance(year_frequencies$Freq, simulated_values))
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
  tolerance = c(160, 80, 40, 20, 10, 5, 2.5, 1)*10
  epsilon_sample = list(c("unif", 0, 0.05))
  ABC_Beaumont <- ABC_sequential(method="Beaumont", model=runSimulation, prior=epsilon_sample, 
                                 nb_simul=100, summary_stat_target = (0), tolerance_tab = tolerance, verbose = TRUE)
  
  colors = colorRampPalette(brewer.pal(8,"Reds"))(length(tolerance))
  # plot lines converging onto actual line
  for (k in 1:length(tolerance)) {
    #new_epsilon = median(ABC_Beaumont[["intermediary"]][[i]][["posterior"]][,2])
    rej <- abc(c(0), ABC_Beaumont[["intermediary"]][[k]][["posterior"]][,2], ABC_Beaumont[["intermediary"]][[k]][["posterior"]][,3], tol = 0.2, method = "rejection")
    #print(median(rej$unadj.values))
    #print(new_epsilon)
    params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = median(rej$unadj.values), rnovo = rnovo, N=N);
    r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span*365)
    a<- c()
    for (j in seq(0, (time_span-1)*365, 365)) {
      a <- c(a, r[which(abs(r[,"time"]-j)==min(abs(r[, "time"]-j))), ]["X"])
    }
    matlines(r[,"time"], r[,"X"], col = colors[k]) 
    #matlines(seq(0, (time_span-1)*365, 365), a, col = colors[k])
  }
  rej <- abc(0, ABC_Beaumont$param, ABC_Beaumont$stats, tol = 0.2, method = "rejection")
  title(continent, line = -2)
  print(median(rej$unadj.values))
  text(time_span*365/2, max(r[,"X"])/2+1000, median(rej$unadj.values))
}

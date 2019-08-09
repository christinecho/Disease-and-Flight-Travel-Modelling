library(adaptivetau)
library(EasyABC)

# get all mcr1 data points from China
china_pts <- subset(mcr1, Country.Codes == "CHN" & !is.na(Formatted_Dates))
china_pts$Year <- as.numeric(format(china_pts$Formatted_Dates,'%Y'))
first_pt <- min(china_pts$Year)
time_span <- as.numeric(max(china_pts$Year)) - as.numeric(min(china_pts$Year))
par(mfrow = c(1, 1))

# get prevalences from China, set year values relative to first detection
prevalences <- read.csv("/Users/christinecho/Documents/Annual-data/PEI/IHME-GBD_2017_DATA-ebca34f7-1/IHME-GBD_2017_DATA-ebca34f7-1.csv")
location_code <- countrycode(prevalences$location_name, origin = "country.name", destination = "iso3c")
prevalences <- cbind(prevalences, location_code)
china_prevalences <- subset(prevalences, location_code == "CHN")
year_diff <- as.numeric(china_pts$Year) - first_pt
year_diff <- unique(year_diff)
y_vals <- sapply(unique(china_pts$Year), function(x) china_prevalences[china_prevalences$year == x, ]$val)

# initial values
N = subset(pop_data, iso3c == "CHN")$value
W = as.integer(china_prevalences[china_prevalences$year == first_pt, ]$val) # as.integer(.1*N)
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
matpoints(year_diff, y_vals, pch = 4)

# plot initial simulation values
for (i in 0:time_span) {
  matpoints(i, r[which(abs(r[,"time"]-i)==min(abs(r[, "time"]-i))), ]["W"], pch = 2)
}


# function to run for abc sequential
# returns Euclidean distance from simulated and experimental values
runSimulation <- function(new_epsilon) {
  params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = new_epsilon, rnovo = rnovo, N=N);
  r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
  simulated_values <- c()
  for (i in 0:time_span) {
    simulated_values <- c(simulated_values, r[which(abs(r[,"time"]-i)==min(abs(r[, "time"]-i))), ]["W"])
  }
  return(calculateDistance(simulated_values, y_vals))
}

# returns sum of Euclidean distance between experimental and simulated values
calculateDistance <- function(x, y){
  dist = 0
  for(i in 1:length(x)) {
    #dist = dist + (x[i]-y[i])**2
    dist = dist + abs(x[i] - y[i])
  }
  return(dist)
}

# ABC model
tolerance = c(160, 80, 40, 20, 10, 5, 2.5, 1)/10
epsilon_sample = list(c("unif", 0, 0.1))
ABC_Beaumont <- ABC_sequential(method="Beaumont", model=runSimulation, prior=epsilon_sample, 
                              nb_simul=100, summary_stat_target = (0), tolerance_tab = tolerance, verbose = TRUE)

# plot lines converging onto actual line
for (i in 1:length(tolerance)) {
  new_epsilon = median(ABC_Beaumont[["intermediary"]][[i]][["posterior"]][,2])
  params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = new_epsilon, rnovo = rnovo, N=N);
  r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time_span)
  print(r[which(abs(r[,"time"]-i)==min(abs(r[, "time"]-i))), ]["W"])
  a<- list()
  for (j in 0:time_span) {
    a[j+1] <- r[which(abs(r[,"time"]-j)==min(abs(r[, "time"]-j))), ]["W"]
  }
  matlines(0:time_span, a, col = i+2)
  
}

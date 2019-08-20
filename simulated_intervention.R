# load in all relevant datasets
mcr1_df <- read.csv("/Users/christinecho/Documents/Annual-data/PEI/mcr1_edited.csv", stringsAsFactors = FALSE)
gdp_per_capita <- read.csv("/Users/christinecho/Documents/Annual-data/gdpdata.csv", stringsAsFactors = FALSE)
gdp_per_capita <- cbind.data.frame(gdp_per_capita, countrycodes)
mcr1_df$Date <- as.integer(mcr1_df$Date)
mcr1_df <- subset(mcr1_df, !is.na(Date))
allcountries <- unique(mcr1_df$Country.Codes)
prevalence_min <- 1e-4
prevalence_max <- 1e-3
par(mfrow = c(2, 2))
par(mar = c(1,2,1,1))
outbreaks_per_country <- table(mcr1_df$Country.Codes)
total_outbreaks <- sum(outbreaks_per_country)

for(i in 1:length(allcountries)) {

# run initial simulation, with normal epsilon. Get value from dataframe.
  # other initial values
  W = as.integer(country_prevalence_val * N) #maybe output prevalence values in dataframe?
  S = as.integer(N-W-X)
  epsilon = 0.05 # this value will be from a future dataframe
  rt = 1/7
  rx = 1/13
  rw = 1/14
  beta = N/S*(epsilon*rt+(1-epsilon)*rw)
  Csp = 1e-16
  rnovo = 1^-6
  
  init.values = c(
    S = S,   # susceptible humans
    W = W,   # infected wild type humans
    X = X)   # infected resistant humans
  
# get date of first emergence

# get date where x = w
  # find the minimum value of abs(x-w), and get time from that row
  t_half_normal <- r[which.min(abs(r[,"W"] - r[,"X"])), ][1]
  
# repeat above process with reduced epsilon

# compare t_halves

# for each simulation, add up columns of x and w from t = 0 to t_half
}  
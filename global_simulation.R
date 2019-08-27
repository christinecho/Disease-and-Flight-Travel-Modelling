library(adaptivetau)

mcr1_df = read.csv("/Users/christinecho/Documents/Annual-data/PEI/mcr1_edited.csv", stringsAsFactors = FALSE)
countryData = read.csv("/Users/christinecho/Documents/Annual-data/PEI/countrydata.csv", stringsAsFactors = FALSE)
countries = countryData$iso3c

# clean data set, remove data points with NAs and unavailable population/gdp
mcr1_df = subset(mcr1_df, !is.na(Year) & !is.na(Country.Codes))
mcr1_df = subset(mcr1_df, Country.Codes %in% countries)

# passenger volume matrix
volume = merge(prediction, AirPortInfoWithCountries[, c("NodeName", "codes")], by.x=c("Origin"), by.y=c("NodeName"))
volume = merge(volume, AirPortInfoWithCountries[, c("NodeName", "codes")], by.x=c("Destination"), by.y=c("NodeName"))
names(volume)[names(volume) == "codes.x"] <- "origin_code"
names(volume)[names(volume) == "codes.y"] <- "destination_code"
volume = subset(volume, volume$origin_code %in% countries & volume$destination_code %in% countries)
volume = volume[volume$origin_code != volume$destination_code, ]
volume = with(volume, tapply(PredMu, list(origin_code, destination_code), FUN = sum))
volume[is.na(volume)] = 0
colnames(volume) = paste0("O",1:numCountries)
rownames(volume) = paste0("D",1:numCountries)

# Model Parameters
numCountries = nrow(countryData)
numDomesticTrans = 5
numInternationalTrans = 1
time = max(mcr1_df$Year) - min(mcr1_df$Year) + 1

# Epidemic Parameters
rt = 1/7/365
rx = 1/13/365
rw = 1/14/365
rnovo = 1e-6

beta = rep(NA, numCountries)
prevalence = rep(NA, numCountries)
epsilon = rep(NA, numCountries)
minPrevalence = 1e-4
maxPrevalence = 1e-3
for (i in 1:numCountries) {
  # prevalences based on gdp per capita, gdp and prevalence inversely correlated
  prevalence[i] = prevalence_min + (prevalence_max - prevalence_min) * countryData$rank[i]/nrow(countryData)
  # TODO: Calculate beta and epsilon. Need to ask Thom about this.
}

params = list(beta = beta, rt = rt, rx = rx, rw = rw, epsilon = epsilon, rnovo = rnovo)
#########################################
# Initial Conditions
#########################################
N = countryData$value
W = as.integer(N * prevalence)
S = N-W
X = rep(0, numCountries)

init.values = c(
  S = S,
  W = W,
  X = X
)
#########################################
# Transitions
#########################################
domesticTrans = c("S=+1, W=-1", # infected wild type to susceptible
                  "S=-1, W=+1", # susceptible to infected wild type
                  "W=-1, X=+1", # infected wild type to infected resistant
                  "S=+1, X=-1", # infected resistant to susceptible
                  "S=-1, X=+1" # susceptible to infected resistant
)

# write text files for transitions
# domestic transitions
domesticTransExpress = NULL
for(i in 1:numCountries) {
  countryExpress = rep(NA, numDomesticTrans)
  for(j in 1:numDomesticTrans) {
    InChain = unlist(strsplit(domesticTrans[j], "="))
    InChainlast = InChain[length(InChain)]
    countryExpress[j] = paste0(paste0(InChain[1:(length(InChain)-1)],i,"=", collapse =""), InChainlast)
  }
  domesticTransExpress = c(domesticTransExpress, countryExpress)
}

# international transitions
internationalTransExpress = NULL
for (i in 1:numCountries) {
  countryExpress = rep(NA, numCountries-1)
  count = 1
  for(j in 1:numCountries) {
    if(i==j) next
    str1 = paste0(c("X",i), collapse = "")
    str2 = paste0(c("=+1, X", j), collapse = "")
    str3 = "=-1"
    countryExpress[count] = paste0(c(str1, str2, str3), collapse = "")
    count = count + 1
  }
  internationalTransExpress = c(internationalTransExpress, countryExpress)
}

allTransExpress = c(domesticTransExpress, internationalTransExpress)
#########################################
# Rates
#########################################
domesticRates = c(
  "(x['W']*(params$epsilon['E']*params$rt + (1 - params$epsilon['E'])*params$rw))", 
  "(params$beta['B']*x['S']*x['W']/params$N['N'])",
  "(params$epsilon['E']*x['W']*params$rnovo)",
  "(x['X']*params$rx)",
  "(params$beta['B']*x['S']*x['X']/params$N['N'])"
)
internationalRates = "(params$Csp * x['S']*x['X']/params$N['N']*params$volume['D','O'])"

# create text file for rates
# domestic rates
domesticRatesExpress = NULL
for(i in 1:numCountries) {
  countryExpress = rep(NA, numDomesticTrans)
  for(j in 1:numDomesticTrans) {
    InChain = unlist(strsplit(domesticRates[j], "']"))
    InChainlast = InChain[length(InChain)]
    countryExpress[j] = paste0(paste0(InChain[1:(length(InChain)-1)],i,"']", collapse =""), InChainlast)
  }
  domesticRatesExpress = c(domesticRatesExpress, countryExpress)
}

# international rates
internationalRatesExpress = NULL
for (i in 1:numCountries) {
  countryExpress = rep(NA, numCountries-1)
  count = 1
  for(j in 1:numCountries) {
    if (i == j) next
    InChain = unlist(strsplit(internationalRates, "']"))
    InChainlast = InChain[length(InChain)]
    str = paste0(paste0(InChain[1:(length(InChain)-1)],i,"']", collapse =""), InChainlast)
    InChain = unlist(strsplit(str, "',"))
    InChainlast = InChain[length(InChain)]
    countryExpress[count] = paste0(paste0(InChain[1:(length(InChain)-1)],j,"',", collapse =""), InChainlast)
    count = count + 1
  }
  internationalRatesExpress = c(internationalRatesExpress, countryExpress)
}

allRatesExpress = c(domesticRatesExpress, internationalRatesExpress)


# Write the transitions list
allTransitions = paste0("c(", allTransExpress, "),")
allTransitions[length(allTransitions)] =  substr(allTransitions[length(allTransitions)],1,(nchar(allTransitions[length(allTransitions)])-1))
transChar = c("transitions = list(", allTransitions, ")")
write.table(transChar, file = paste0(getwd(),"/expressions.R"), row.names = F, col.names = F, quote = F)

# Write the rates function
allRates = paste0(allRatesExpress, ",")
allRates[length(allRates)] =  substr(allRates[length(allRates)],1,(nchar(allRates[length(allRates)])-1))
ratesChar = c("rates <- function(x, params, t) {", "return(c(", allRates, "))}")
write.table(ratesChar, file = paste0(getwd(),"/expressions.R"), append = T, row.names = F, col.names = F, quote = F)

# Run simulation
source(paste0(getwd(),"/expressions.R"))
r = ssa.adaptivetau(init.values, transitions, rates, params, tf=time)

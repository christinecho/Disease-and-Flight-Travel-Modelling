mcr1_df = read.csv("/Users/christinecho/Documents/Annual-data/PEI/mcr1_edited.csv", stringsAsFactors = FALSE)
countryData = read.csv("/Users/christinecho/Documents/Annual-data/PEI/countrydata.csv", stringsAsFactors = FALSE)

# Clean data set, remove data points with NAs and unavailable population/gdp
mcr1_df = subset(mcr1_df, !is.na(Year) & !is.na(Country.Codes))
mcr1_df = subset(mcr1_df, Country.Codes %in% countryData$iso3c)

# Model Parameters
numCountries = nrow(countryData)
numDomesticTrans = 5
numInternationalTrans = 1

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

#########################################
# Initial Conditions
#########################################
N = countryData$value
W = N * prevalence
S = N-W

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
    str1 = paste(c("X",i), collapse = "")
    str2 = paste(c("=+1, X", j), collapse = "")
    str3 = "=-1"
    countryExpress[count] = paste(c(str1, str2, str3), collapse = "")
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

internationalRates = "params$Csp * x['S']*x['X']/params$N['N']*params$passengers['A','B']"

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

internationalRatesExpress = NULL
for (i in 1:numCountries) {
  countryExpress = rep(NA, numCountries-1)
  count = 1
  for(j in 1:numCountries) {
    InChain = unlist(strsplit())
  }
}

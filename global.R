# read in all datasets
mcr1_df <- read.csv("/Users/christinecho/Documents/Annual-data/PEI/mcr1_edited.csv", stringsAsFactors = FALSE)
mcr1_df <- subset(mcr1_df, !is.na(mcr1_df$Country.Codes))
mcr1_df <- subset(mcr1_df, !is.na(mcr1_df$Date))
predictions <- mean_predic
countries <-  unique(AirPortInfoWithCountries$codes)
total_outbreaks <- nrow(mcr1_df)
X_entering <- data.frame(Country = character(), X = numeric(), stringsAsFactors = FALSE)
gdp_per_capita <- read.csv("/Users/christinecho/Documents/Annual-data/gdpdata.csv", stringsAsFactors = FALSE)
countrycodes <- countrycode(gdp_per_capita$country, origin = "country.name", destination = "iso3c")
gdp_per_capita <- cbind.data.frame(gdp_per_capita, countrycodes)
prevalence_min <- 1e-4
prevalence_max <- 1e-3

# get sum of people going to each country
merged_df <- merge(prediction, AirPortInfoWithCountries[, c("NodeName", "codes")], by.x=c("Origin"), by.y=c("NodeName"))
merged_df <- merge(merged_df, AirPortInfoWithCountries[, c("NodeName", "codes")], by.x=c("Destination"), by.y=c("NodeName"))
names(merged_df)[names(merged_df) == "codes.x"] <- "origin_code"
names(merged_df)[names(merged_df) == "codes.y"] <- "destination_code"
international <- merged_df[merged_df$origin_code != merged_df$destination_code, ]
sum_predic <- aggregate(PredMu ~ origin_code+destination_code, data = international, sum)
first_outbreak = min(mcr1_df$Date)
last_outbreak = max(mcr1_df$Date)
time_span = last_outbreak - first_outbreak + 1
offset = first_outbreak + 1

# vector of actual number of outbreaks per day
days <- seq(min(mcr1_df$Date), max(mcr1_df$Date))
days <- data.frame(days)
day_frequencies <- as.data.frame(table(mcr1_df$Date))
day_frequencies$Var1 = as.numeric(levels(day_frequencies$Var1))[day_frequencies$Var1]
day_frequencies <- merge(days, day_frequencies, by.x = "days", by.y = "Var1", all=T)
day_frequencies[is.na(day_frequencies)] <- 0

# ABC sequential model
tolerance = c(160, 80, 40, 20, 10, 5, 2.5, 1)*100
Csp_sample = list(c("unif", 0, 1))
ABC_Beaumont <- ABC_sequential(method="Beaumont", model=calculateDiffInX, prior=Csp_sample, 
                               nb_simul=100, summary_stat_target = (0), tolerance_tab = tolerance, verbose = TRUE)

# returns sum of sum of squares
calculateDiffInX <- function(Csp) {
  print(Csp)
  sumSquaresCountries <- c()
  sumofSumSquares <- 0
  for (i in 1:length(allcountries)) {
    if ((!countries[i] %in% pop_data$iso3c) | (!countries[i] %in% gdp_per_capita$countrycodes)) next
    
    # add domestic outbreaks for each day
    X_total <- numeric(time_span)
    domestic_outbreaks = subset(mcr1_df, mcr1_df$Country.Codes == countries[i])
    for (m in domestic_outbreaks$Date) {
      X_total[m - offset] = X_total[m - offset] + 1
    }
      
    incomingX <- numeric(time_span)
      
    # get number of outbreaks for every other country on a given day
      for (j in 1:length(countries)) {
        if (i == j) next
        international_outbreaks = subset(mcr1_df, mcr1_df$Country.Codes == countries[j])
        if (nrow(international_outbreaks) == 0) next
        if ((!countries[j] %in% pop_data$iso3c)) next
        Nj = subset(pop_data, iso3c == countries[j])$value
        international_outbreaks <- table(international_outbreaks$Date)
        dates <- as.integer(rownames(international_outbreaks))
        for (k in 1:length(dates)) {
          if (nrow(subset(sum_predic, sum_predic$destination_code == countries[i] & sum_predic$origin_code == countries[j])) == 0) {
            people_traveling = 0
          } else {
            people_traveling = subset(sum_predic, sum_predic$destination_code == countries[i] & sum_predic$origin_code == countries[j])$PredMu/365
          }
          incomingX[dates[k] - offset] = incomingX[dates[k] - offset] + international_outbreaks[[k]]/Nj * people_traveling
        }        
      }
    
    # calculate incoming outbreaks
    Ni = subset(pop_data, iso3c == countries[i])$value
    country_rank <- gdp_per_capita[which(gdp_per_capita$countrycodes == curr_country_code), ]$rank
    country_prevalence_val <- prevalence_min + (prevalence_max - prevalence_min) * country_rank/nrow(gdp_per_capita)
    S <- Ni - country_prevalence_val*Ni
    incomingX <- incomingX * Csp * S
    X_total = X_total + incomingX
    sumofSumSquares = sumofSumSquares + calculateSumSquares(day_frequencies$Freq, X_total)
  }
  return(sumofSumSquares)
}

# calculates sum of squares between actual and simulated number of outbreaks
calculateSumSquares <- function(actualVector, simulatedVector) {
  sum = 0
  for (i in 1:length(actualVector)) {
    sum = sum + (actualVector[i] - simulatedVector[i])**2
  }
  return(sum)
}
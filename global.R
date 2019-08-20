mcr1_df <- read.csv("/Users/christinecho/Documents/Annual-data/PEI/mcr1_edited.csv", stringsAsFactors = FALSE)
mcr1_df <- subset(mcr1_df, !is.na(mcr1_df$Country.Codes))
mcr1_df <- subset(mcr1_df, !is.na(mcr1_df$Date))
predictions <- mean_predic
countries <-  unique(AirPortInfoWithCountries$codes)
total_outbreaks <- nrow(mcr1_df)
X_entering <- data.frame(Country = character(), X = numeric(), stringsAsFactors = FALSE)

# get sum of people going to each country
merged_df <- merge(prediction, AirPortInfoWithCountries[, c("NodeName", "codes")], by.x=c("Origin"), by.y=c("NodeName"))
merged_df <- merge(merged_df, AirPortInfoWithCountries[, c("NodeName", "codes")], by.x=c("Destination"), by.y=c("NodeName"))
names(merged_df)[names(merged_df) == "codes.x"] <- "origin_code"
names(merged_df)[names(merged_df) == "codes.y"] <- "destination_code"
international <- merged_df[merged_df$origin_code != merged_df$destination_code, ]
sum_predic <- aggregate(PredMu ~ origin_code+destination_code, data = international, sum)
first_outbreak = min(mcr1_df$Date)
last_outbreak = max(mcr1_df$Date)

# vector of actual number of outbreaks per day
days <- seq(min(mcr1_df$Date), max(mcr1_df$Date))
days <- data.frame(days)
day_frequencies <- as.data.frame(table(mcr1_df$Date))
day_frequencies$Var1 = as.numeric(levels(day_frequencies$Var1))[day_frequencies$Var1]
day_frequencies <- merge(days, day_frequencies, by.x = "days", by.y = "Var1", all=T)
day_frequencies[is.na(day_frequencies)] <- 0

# ABC model
# tolerance = c(160, 80, 40, 20, 10, 5, 2.5, 1)*100
#Csp_sample = list(c("unif", 0, 1e-4))
#ABC_Beaumont <- ABC_sequential(method="Beaumont", model=calculateDiffInX, prior=Csp_sample, 
                               #nb_simul=100, summary_stat_target = (0), tolerance_tab = tolerance, verbose = TRUE)

calculateDiffInX <- function(Csp) {
  sumSquaresCountries <- c()
  for (i in 1:1) {
    domestic_outbreaks = subset(mcr1_df, mcr1_df$Country.Codes == countries[i])
    incomingX <- c()
    for (k in first_outbreak:last_outbreak) {
      domestic_outbreak_day = nrow(subset(domestic_outbreaks, domestic_outbreaks$Date == k))
      X = domestic_outbreak_day
      if (nrow(subset(mcr1_df, mcr1_df$Date == k))==0) {
        incomingX[k-first_outbreak+1] <- Csp * X
        next
      }
      for (j in 1:length(countries)) {
        if (i == j) next
        international_outbreaks = subset(mcr1_df, mcr1_df$Country.Codes == countries[j])
        if ((!countries[j] %in% pop_data$iso3c)) next
        Nj = subset(pop_data, iso3c == countries[j])$value
        if (nrow(subset(sum_predic, sum_predic$destination_code == countries[i] & sum_predic$origin_code == countries[j])) == 0) {
          people_traveling = 0
        } else {
          people_traveling = subset(sum_predic, sum_predic$destination_code == countries[i] & sum_predic$origin_code == countries[j])$PredMu/365
        }
        international_outbreak_day = nrow(subset(international_outbreaks, international_outbreaks$Date == k))
        X = X + international_outbreak_day/Nj * people_traveling
      }
      incomingX[k-first_outbreak+1] <- Csp * X
    }
    print(incomingX)
    sumSquaresCountries[i] <- calculateSumSquares(day_frequencies$Freq, incomingX)
  }
  return(sum(sumSquaresCountries))
}

calculateSumSquares <- function(actualVector, simulatedVector) {
  sum = 0
  for (i in 1:length(actualVector)) {
    sum = sum + (actualVector[i] - simulatedVector[i])**2
  }
  return(sum)
}

library(DescTools)
library(rworldmap)

# add country codes
merged_df <- merge(prediction, AirPortInfoWithCountries[, c("NodeName", "codes")], by.x=c("Origin"), by.y=c("NodeName"))
merged_df <- merge(merged_df, AirPortInfoWithCountries[, c("NodeName", "codes")], by.x=c("Destination"), by.y=c("NodeName"))

#rename columns
names(merged_df)[names(merged_df) == "codes.x"] <- "origin_code"
names(merged_df)[names(merged_df) == "codes.y"] <- "destination_code"


# remove in-country flights
international <- merged_df[merged_df$codes.x != merged_df$codes.y, ]

# aggregate sum
sum_predictions_origin <- aggregate(PredMu ~ codes.x, data = international, sum)
sum_predictions_destination <- aggregate(PredMu ~ codes.y, data = international, sum)
aggregate_sum_predictions<- sum_predictions_destination
aggregate_sum_predictions$PredMu <- sum_predictions_origin$PredMu + sum_predictions_destination$PredMu

# aggregate mean
freq_1 <- Freq(international$codes.x)
freq_2 <- Freq(international$codes.y)
freq_2$freq <- freq_1$freq + freq_2$freq
aggregate_mean_predictions_df <- data.frame(freq_2$level, freq_2$freq)
colnames(aggregate_mean_predictions_df) <- c("code", "frequency")
mean_predictions <- sum_predictions$PredMu / aggregate_mean_predictions_df$frequency
aggregate_mean_predictions_df <- cbind(aggregate_mean_predictions_df, mean_predictions)

# heat map
spdf <- joinCountryData2Map(sum_predictions, joinCode="ISO3", nameJoinColumn="codes.y")
mapCountryData(spdf, nameColumnToPlot="PredMu", catMethod="fixedWidth", mapTitle = "Sum of Travel from International Flights")
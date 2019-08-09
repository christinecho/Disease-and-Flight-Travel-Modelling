# add country codes
merged_df <- merge(prediction, AirPortInfoWithCountries[, c("NodeName", "codes")], by.x=c("Origin"), by.y=c("NodeName"))
merged_df <- merge(merged_df, AirPortInfoWithCountries[, c("NodeName", "codes")], by.x=c("Destination"), by.y=c("NodeName"))

#rename columns
names(merged_df)[names(merged_df) == "codes.x"] <- "origin_code"
names(merged_df)[names(merged_df) == "codes.y"] <- "destination_code"


# remove in-country flights
international <- merged_df[merged_df$origin_code != merged_df$destination_code, ]

# aggregate mean predictions
mean_predic <- aggregate(PredMu ~ origin_code+destination_code, data = international, sum)
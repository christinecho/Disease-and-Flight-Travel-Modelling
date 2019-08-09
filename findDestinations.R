findDestinationCountries = function(destination)
{
  # find if destination is in India by looking up OAGName in AirPort Info
  # if so, add prediction to origin
}

indiaData <- subset(AirPortInfoWithCountries, countries=='India')
indiaDestination <- subset(prediction, Destination %in% indiaData$NodeName)
indiaOrigin <- subset(prediction, Origin %in% indiaData$NodeName)

# get total predictions to country from Indian airports
predictionsum <- aggregate(indiaDestination$PredMu~indiaDestination$Origin, FUN = sum)

library(ggmap)
library(plyr)
register_google(key="AIzaSyDEXzboKO2DG9GC-4QBEvTYQdkdEzr5oIg", account_type = "standard")
# cities <- as.data.frame(ndm1$Country)
# cities <- as.data.frame(lapply(cities, as.character), stringsAsFactors = FALSE)
# b <- apply(cities, 1, geocode, )
#cities<-c("Hubei", "Changsha", "Shandong")
# a <- geocode(cities$`mcr$Region`)
a <- data.frame()
for(i in seq(145, nrow(cities))) {
  print(ndm1$Country[i])
  print(geocode(ndm1$Country[i]))
}
#library(rworldmap)
#map <- getMap()
#mapPoints <- ggmap(map) + geom_point(aes(x = AirportInfo$Lon, y = AirportInfo$Lat, size=))
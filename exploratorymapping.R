library(rworldmap)
worldMap <- getMap(resolution = "coarse")
plot(worldMap)

points(AirportInfo$Lon, AirportInfo$Lat, col="gray", cex=0.1)
points(mcr1$Longitude, mcr1$Latitude, col = "red", cex = 0.3, pch = 22)
points(ndm$Longitude, ndm$Latitude, col = "blue", cex = 0.3, pch = 24)
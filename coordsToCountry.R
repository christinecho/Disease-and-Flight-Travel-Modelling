library(sp)
library(rworldmap)

#convert coordinates to countries
coords2country = function(points)
{
  countriesSP <- getMap(resolution = "low")
  pointsSP = SpatialPoints(points, proj4string = CRS(proj4string(countriesSP)))
  
  indices = over(pointsSP, countriesSP)
  
  indices$ADMIN
}

#https://stackoverflow.com/questions/14334970/convert-latitude-and-longitude-coordinates-to-country-name-in-r

#country frequency table
countries<-coords2country(data.frame(lon=AirportInfo$Lon, lat=AirportInfo$Lat))
countries.freq<-as.data.frame(table(countries))
colnames(countries.freq) <- c("country", "value")

#match to world map and plot heat map
matched <- joinCountryData2Map(countries.met, joinCode="NAME", nameJoinColumn="country")
mapCountryData(matched, nameColumnToPlot="value", mapTitle="Country Airport Sample", catMethod = "pretty", colourPalette = "heat")
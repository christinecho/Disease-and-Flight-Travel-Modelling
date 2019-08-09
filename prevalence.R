library(countrycode)

prevalences <- read.csv("/Users/christinecho/Documents/Annual-data/PEI/IHME-GBD_2017_DATA-385b67b1-1/IHME-GBD_2017_DATA-385b67b1-1.csv")

location_code <- countrycode(prevalences$location_name, origin = "country.name", destination = "iso3c")
prevalences <- cbind(prevalences, location_code)


library(ggplot2)
library(ggmap)
library(maps)


mcr1_map = ggplot() + borders("world", colour = "gray85", fill = "gray80") + 
  geom_point(aes(x = Longitude, y = Latitude, color = Formatted_Dates), data = mcr1, alpha = 0.2) + 
  scale_color_gradient(trans = "date", high = 'blue', low = 'red')
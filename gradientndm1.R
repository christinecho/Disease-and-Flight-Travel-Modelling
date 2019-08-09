library(ggplot2)
library(ggmap)
library(maps)

ndm_map = ggplot() + borders("world", colour = "gray85", fill = "gray80") + 
  geom_point(aes(x = Longitude, y = Latitude, color = date), data = ndm, alpha = 0.5) +
  scale_color_gradient(trans = "date", high = 'blue', low = 'red')
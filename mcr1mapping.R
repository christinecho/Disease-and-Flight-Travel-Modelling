library(ggplot2)
library(ggmap)
library(maps)
library(ggthemes)
library(gganimate)
library(animation)


# world <- ggplot() + borders("world", colour = "gray85", fill = "gray80") +
  # theme_map()

# map <- ggplot() + geom_point(data = mcr1, aes(x = Longitude, y = Latitude), colour = "red")
# map
ndm_sorted<-ndm[order(ndm$`Date of Detection`),]


# animates map
world <- ggplot() + borders("world", colour = "gray85", fill = "gray80") + theme_map() 
map <- world + geom_point(data = ndm_sorted, aes(x = Longitude, y = Latitude), colour = "red") +
  transition_states(ndm_sorted$`Date of Detection`) +
  labs(title="Date: {closest_state}")

animate(map)
# anim_save("ndm1.gif", anim)
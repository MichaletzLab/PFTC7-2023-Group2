library(mapview) 
library(sf)

##Make the first map that is just the plots pinned by elevation
# Enter Latlons and elevations
lng = c(28.88308149,28.88270548,28.89511749,28.89388199,28.89680451,28.89565154,28.89267001,28.8924665,28.89438299,28.89257249) 
lat = c(-28.66095603,-28.66243452,-28.67986452,-28.68031898,-28.719115,-28.71870303,-28.72832348,-28.72822399,-28.73774635,-28.73774152) 
names = c("1E", "1W","2E", "2W", "3E", "3W","4E", "4W","5E", "5W") 
elevs = c("2000","2000","2200","2200","2400","2400","2600","2600","2800","2800")

points_df = data.frame(lng, lat, elevs) 

# Convert the data frame to a spatial points data frame 
Elevation = st_as_sf(points_df,  
                      coords = c("lng", "lat"), crs = 4326) 

# Plot the points on a map 
mapview(Elevation, label = points_df$elevs)

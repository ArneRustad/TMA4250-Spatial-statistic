library(fields)
library(MASS)
library(ggplot2)

path = "C:\\Users\\Lene\\Documents\\Skole\\romlig\\TMA4250-Spatial-statistics"

# Problem 1

#read the data

data.complit = read.delim(url("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/complit.dat"), 
                        header = FALSE , sep =" ")
data.seismic = read.delim(url("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/seismic.dat"), 
                          header = FALSE , sep ="\t", col.names = c("seismic"))

image.plot()

matrix.seismic = data.matrix(seismic.data)
grid = expand.grid(x = 1:75, y = 1:75)
df.obs = data.frame(grid, obs = data.seismic$seismic)
df.obs

ggplot(data = df.obs, aes(x = x, y = y, fill = obs)) + geom_tile() + scale_fill_distiller(palette = "Spectral")+ 
  labs(fill = "Seismic value") + ggtitle("Image plot of seismic data") + theme_minimal() + xlim(c(0,75)) + ylim(c(0,75))
ggsave("imageplot_seismic.pdf", width = 5, height = 4, path = path)



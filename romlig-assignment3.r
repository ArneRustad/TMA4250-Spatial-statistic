library(fields)
library(MASS)
library(ggplot2)

path = "C:\\Users\\Lene\\Documents\\Skole\\romlig\\TMA4250-Spatial-statistics"

# Problem 1

#a

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


#b)

phi0 = function(d){
  return(pnorm(d,0.02, 0.06))
}

phi1 = function(d){
  return(pnorm(d,0.08, 0.06))
}

sim.posterior.given.unif.prior = function(d){
  p = phi0(d)/(phi1(d) + phi0(d))    # probability for 0
  u = runif(length(p))
  realization = ifelse(u<p, 0, 1)
  grid = expand.grid(x = 1:75, y = 1:75)
  return(data.frame(grid, realization = factor(realization)))
}


set.seed(1)
for(i in 1:6){
  df.real = sim.posterior.given.unif.prior(df.obs$obs)
  p = ggplot(data = df.real, aes(x = x, y = y, fill = realization)) + geom_tile() + 
    labs(fill = "Lithology") + ggtitle(paste("Realization", i,  "of posterior Mosaic RF")) + theme_minimal() + xlim(c(0,75)) + 
    ylim(c(0,75)) + scale_fill_manual(values = c("#3C8EC1", "#DF5452"),labels = c("Sand - 0", "Shale - 1"))
  print(p)
  ggsave(paste("realization",i,".pdf", sep=""), width = 5, height = 4, path = path)
}


# expectation plot

expectation.mosaic.rf = function(d){
  p = phi1(d)/(phi1(d) + phi0(d))    # probability for 1
  grid = expand.grid(x = 1:75, y = 1:75)
  return(data.frame(grid, expectation = p))
}

df.expectation = expectation.mosaic.rf(df.obs$obs)
ggplot(data = df.expectation, aes(x = x, y = y, fill = expectation)) + geom_tile() + 
  labs(fill = "Lithology") + ggtitle("Expectation of posterior Mosaic RF") + theme_minimal() + xlim(c(0,75)) + ylim(c(0,75)) +
  scale_fill_distiller(palette = "Spectral")
ggsave("posterior_expectation_map.pdf", width = 5, height = 4, path = path)


# variance plot

var.mosaic.rf = function(d){
  p = phi1(d)/(phi1(d) + phi0(d))    # probability for 1
  var = p*(1-p)
  grid = expand.grid(x = 1:75, y = 1:75)
  return(data.frame(grid, variance = var))
}

df.var = var.mosaic.rf(df.obs$obs)
ggplot(data = df.var, aes(x = x, y = y, fill = variance)) + geom_tile() + 
  labs(fill = "Lithology") + ggtitle("Variance of posterior Mosaic RF") + theme_minimal() + xlim(c(0,75)) + ylim(c(0,75)) +
  scale_fill_distiller(palette = "Spectral")
ggsave("posterior_var_map.pdf", width = 5, height = 4, path = path)


# maximum marginal posterior predictor

get.mmap = function(d){
  p = phi1(d)/(phi1(d) + phi0(d))    # probability for 1
  mmap = ifelse(p>= 0.5, 1, 0)
  grid = expand.grid(x = 1:75, y = 1:75)
  return(data.frame(grid, mmap = factor(mmap)))
}

df.mmap = get.mmap(df.obs$obs)

ggplot(data = df.mmap, aes(x = x, y = y, fill = mmap)) + geom_tile() + 
  labs(fill = "Lithology") + ggtitle("Maximum marginal posterior predictor") + theme_minimal() + xlim(c(0,75)) + ylim(c(0,75)) +
  scale_fill_manual(values = c("#3C8EC1", "#DF5452"),labels = c("Sand - 0", "Shale - 1"))
ggsave("mmap_map.pdf", width = 5, height = 4, path = path)

library(fields)
library(MASS)
library(ggplot2)
library(tidyverse)

path = "C:\\Users\\Lene\\Documents\\Skole\\romlig\\TMA4250-Spatial-statistics"

# Problem 1

#a

#read the data

data.complit = read.delim(url("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/complit.dat"), 
                        header = FALSE , sep =" ")
data.seismic = read.delim(url("https://www.math.ntnu.no/emner/TMA4250/2020v/Exercise3/seismic.dat"), 
                          header = FALSE , sep ="\t", col.names = c("seismic"))

matrix.seismic = data.matrix(data.seismic)
grid = expand.grid(x = 1:75, y = 1:75)
df.obs = data.frame(grid, obs = data.seismic$seismic)
df.obs

ggplot(data = df.obs, aes(x = x, y = y, fill = obs)) + geom_tile() + scale_fill_distiller(palette = "Spectral")+ 
  labs(fill = "Seismic value") + ggtitle("Image plot of seismic data") + theme_minimal() + xlim(c(0,75)) + ylim(c(0,75))
ggsave("imageplot_seismic.pdf", width = 5, height = 4, path = path)


#b)

phi0 = function(d){
  return(dnorm(d,0.02, 0.06))
}

phi1 = function(d){
  return(dnorm(d,0.08, 0.06))
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




grid.small = expand.grid(x = 1:66, y = 1:66)

matrix.complit = data.matrix(data.complit)

complit.long = cbind(grid.small, obs = NA)

for(y in 1:66){
  complit.long[(1+(y-1)*66):(y*66),3] = matrix.complit[(67-y),]
}

complit.long.df  = data.frame(complit.long)

ggplot(data = complit.long.df, aes(x = x, y = y, fill = factor(obs))) + geom_tile() + 
  labs(fill = "Lithology") + ggtitle("Observations of the lithology distribution") + theme_minimal() + xlim(c(0,66)) + ylim(c(0,66)) +
  scale_fill_manual(values = c("#3C8EC1", "#DF5452"),labels = c("Sand - 0", "Shale - 1"))
ggsave("obs_distribution.pdf", width = 5, height = 4, path = path)


## pseudo likelihood inference

mmpl.estimator = function(training.image) {
  count.similar = matrix(0, ncol = ncol(training.image), nrow = nrow(training.image))
  count.different = matrix(0, ncol = ncol(training.image), nrow = nrow(training.image))
  for (i in 1:nrow(training.image)) {
    for (j in 1:ncol(training.image)) {
      if (i + 1 <= 66) {
        if (training.image[i,j] == training.image[i+1,j]) count.similar[i,j] = count.similar[i,j] + 1
        else count.different[i,j] = count.different[i,j] + 1
      }
      if (i - 1 >= 1) {
        if (training.image[i,j] == training.image[i-1,j]) count.similar[i,j] = count.similar[i,j] + 1
        else count.different[i,j] = count.different[i,j] + 1
      }
      if (j + 1 <= 66) {
        if (training.image[i,j] == training.image[i,j+1]) count.similar[i,j] = count.similar[i,j] + 1
        else count.different[i,j] = count.different[i,j] + 1
      }
      if (j - 1 >= 1) {
        if (training.image[i,j] == training.image[i,j-1]) count.similar[i,j] = count.similar[i,j] + 1
        else count.different[i,j] = count.different[i,j] + 1
      }
    }
  }
  count.similar.long = as.vector(count.similar)
  count.different.long = as.vector(count.different)
  
  mpl = function(beta, count.similar.long, count.different.long) {
    return(sum(count.similar.long * log(beta) - log(beta^(count.similar.long) + beta^(count.different.long) ) ))
  }
  #beta.hat = optim(1.2, fn = mpl, gr = NULL, training.image)
  
  mpl(1, count.similar.long, count.different.long)
  beta.hat = optim(1.2, fn = mpl, gr = NULL, method = "Brent",
                   count.similar.long, count.different.long,
                   lower = 1, upper = 1e10,
                   control = list(fnscale = -1))
  print(beta.hat$value)
  print(mpl(3, count.similar.long, count.different.long))
  print(mpl(4, count.similar.long, count.different.long))
  return(beta.hat$par)
}

beta.hat = mmpl.estimator(matrix.complit)
beta.hat


map.obs = matrix(NA, nrow = max(df.obs$y), ncol = max(df.obs$x))
n = max(df.obs$y)
for (j in 1:nrow(map.obs)) {
  map.obs[j,] = df.obs$obs[df.obs$y == j]
}
image.plot(map.obs)


# MCMC

gibbs.sim = function(beta, map.obs, n.sweep = 25, start = "random", output = "realization") {
  n = nrow(map.obs)
  n.squared = n^2
  if (start == "random") l.initial = sample(c(0,1), n.squared, replace=TRUE)
  else {
    if (start == "0") l.initial = rep(0, n.squared)
    else {
      if (start == "1") l.initial = rep(1, n.squared)
      else l.initial = start
    }
  }
  l = matrix(l.initial, ncol = n, nrow = n)
  colnames(l) = 1:n
  rownames(l) = 1:n
  sand.proportion = rep(NA, n.sweep+1)
  sand.proportion[1] =  1 - sum(l) / n.squared
  counter.prop = 1
  for (k in 1:(n.sweep*n.squared)) {
    j = sample(1:n, 1)
    i = sample(1:n, 1)
    #print(paste((l[ifelse(j == n, 1, j+1),j] == 1), l[ifelse(j == n, 1, j+1),j]))
    neighbours.equal.1 = (l[ifelse(j == n, 1, j+1),i] == 1) + (l[ifelse(j == 1, n, j-1), i] == 1) +
      (l[j, ifelse(i == n, 1, i+1)] == 1) + (l[j, ifelse(i == 1, n, i - 1)] == 1)
    #print(neighbours.equal.1)
    #print(l[(j-1):(j+1), (i-1):(i+1)])
    #print(paste("Neigbours=1:", neighbours.equal.1, " phi0:", phi0(map.obs[j,i]), " phi1:", phi1(map.obs[j,i])))
    p = phi1(map.obs[j,i]) * beta^(neighbours.equal.1) / 
      ( phi0(map.obs[j,i]) * beta^(4 - neighbours.equal.1) + phi1(map.obs[j, i]) * beta^(neighbours.equal.1))
    #print(p)
    u = runif(1)
    l[j,i] = ifelse(u < p, 1, 0)
    
    if (k %% n.squared == 0) {
      counter.prop = counter.prop + 1
      sand.proportion[counter.prop] = 1 - sum(l) / n.squared
    }
    
  }
  #plot(seq(1,n.sweep+1, 1), sand.proportion)
  if (output == "convergence") {
    return (sand.proportion)
  }
  return(l)
}

realization = gibbs.sim(beta.hat, map.obs, 100)

# function for wide to long format
wide.to.long = function(m) {
  n = nrow(m)
  grid = expand.grid(x = 1:n, y = 1:n)
  m.long = cbind(grid, obs = NA)
  for(y in 1:n){
    m.long[(1+(y-1)*n):(y*n),3] = m[y,]
  }
  return(data.frame(m.long))
}


df.realization.long = wide.to.long(realization)

ggplot(data = df.realization.long, aes(x = x, y = y, fill = factor(obs))) + geom_tile() + 
  labs(fill = "Lithology") + ggtitle("Observations of the lithology distribution") + theme_minimal() + xlim(c(0,75)) + ylim(c(0,75)) +
  scale_fill_manual(values = c("#3C8EC1", "#DF5452"),labels = c("Sand - 0", "Shale - 1"))
# this should not be saved


####### convergence plot
set.seed(1)
n.sweep = 10
df.convergence = data.frame(sweep = seq(0,n.sweep), 
                            sand.proportion = gibbs.sim(beta.hat, map.obs, n.sweep,
                                                              start = "random", output = "convergence"),
                            initial = "Uniform random 1")
df.convergence = rbind(df.convergence,
                       data.frame(sweep = seq(0,n.sweep),
                                  sand.proportion = gibbs.sim(beta.hat, map.obs, n.sweep,
                                                              start = "random", output = "convergence"),
                                  initial = "Uniform random 2"))
df.convergence = rbind(df.convergence,
                       data.frame(sweep = seq(0,n.sweep),
                                  sand.proportion = gibbs.sim(beta.hat, map.obs, n.sweep,
                                                              start = "random", output = "convergence"),
                                  initial = "Uniform random 3"))
df.convergence = rbind(df.convergence,
                       data.frame(sweep = seq(0,n.sweep),
                                  sand.proportion = gibbs.sim(beta.hat, map.obs, n.sweep,
                                            start = "0", output = "convergence"),
                                  initial = "All sand"))
df.convergence = rbind(df.convergence,
                       data.frame(sweep = seq(0,n.sweep),
                                  sand.proportion = gibbs.sim(beta.hat, map.obs, n.sweep,
                                            start = "1", output = "convergence"),
                                  initial = "All shale"))
ggplot(df.convergence, aes(x = sweep, y = sand.proportion, col = initial)) + geom_line() + ylab("Proportion of sand") + xlab("Sweep")



###### plot of posterior realizations


plot.posterior.realizations = function(beta.hat, map.obs, n.realizations = 6, n.sweep.initial = 50, n.sweep) {
  for (i in range) {
    if (i == 1) {
      realization.i =   gibbs.sim(beta.hat, map.obs, n.sweep.initial, start = "random")
    }
    else {
      realization.i =   gibbs.sim(beta.hat, map.obs, n.sweep, start = "random", start = realization.i)
    }
    df.realization.i = wide.to.long(realization.i)
    ggplot(data = df.realization.i, aes(x = x, y = y, fill = factor(obs))) + geom_tile() + 
      labs(fill = "Lithology") + ggtitle("Observations of the lithology distribution") + theme_minimal() + xlim(c(0,75)) + ylim(c(0,75)) +
      scale_fill_manual(values = c("#3C8EC1", "#DF5452"),labels = c("Sand - 0", "Shale - 1"))
    ggsave("fill_inn_here.pdf", width = 5, height = 4, path = path)
}

}
rm(list = ls())
source("GA_Search.R")
library(nnet)

### Comment the below out if not using shell script to run code and just enter
### values for 'fn' and 'dat.dir' manually..
### ------
# Read in of a string argument to indicate segment
args <- commandArgs(TRUE)
fn <- args[1]   #   filename - e.g. 'bank_8_24fh_1.csv'
dat.dir <- args[2]   #  dataset - e.g. 'Bank'
### ------


# Read in data
# ----
out.dir <- paste("~/", dat.dir, "/", sep = "")


#-------------------------------------------------------------------------------
standardize <- function(Z) {
  Z.mean <- mean(Z)
  Z.sd <- sd(Z)
  Z.stand <- (Z - Z.mean) / Z.sd
  return(Z.stand)
}
#-------------------------------------------------------------------------------

# Read in data
# ----
inp.file <- fn
fn <- unlist(strsplit(fn, split = "/"))
fn <- fn[length(fn)]


names.1 <- scan(file = inp.file, what = character(), sep = ",",
                nlines = 1)
dat <- read.table(inp.file, skip = 1, sep = ",", fill = TRUE)


id <- dat[,1]
id <- dat[,1]
dat <- dat[,-1]
names.1 <- names.1[-1]
# store output variable in 'y'
y <- dat[,ncol(dat)]
# store input variables in 'x'
x <- dat[,1:(ncol(dat) - 1)]

names(y) <- names.1[ncol(dat)]
names(x) <- names.1[1:(ncol(dat) - 1)]

ndata <- length(y)

# Standardize the output data with zero mean and unit variance
y.stand <- standardize(y)
y <- y.stand

# Standardize the input data with zero mean and unit variance
z.stand <- apply(x, 2, standardize)
x <- z.stand

### Define objective functions (used to find optimal inputs)
#--------------------------------------------------------
fitness <- function(par, x, y, gp.ids, k = 1) {
# fitness is based on out-of-sample AIC of an ANN. 
# The influence of the penalty term can be adjusted by altering the parameter k.
# For k=1, the result is the standard AIC.
# Uses leave n out cross validation, where n depends on the size of the dataset
# and is ~20% of the data
  aic <- 0
  for(i in unique(gp.ids)) {
    for(try in 1:10) {
      fit <- nnet(x = as.matrix(x[gp.ids != i,par == 1], ncol = sum(par)),
                  y = y[gp.ids != i], size = 1, rang = 0.1, maxit = 1000, 
                  linout = TRUE, trace = FALSE, method = "SANN")
      if(try == 1) {
        fit.best <- fit
      } else if(fit$value < fit.best$value) {
        fit.best <- fit
      }
    }
    pred <- predict(fit.best, 
                    newdata = as.matrix(x[gp.ids == i,par == 1], ncol = sum(par)),
                    type = "raw")
    err <- sum((y[gp.ids == i] - pred)^2)
    n.data <- length(y[gp.ids == i])
    aic.i <- n.data * log(err / n.data) + k * 2 * length(fit.best$wts)
    aic <- aic + aic.i
  }
  return(-aic)
}
#--------------------------------------------------------
fitness.1 <- function(par, x, y, gp.ids) {
# fitness is based on out-of-sample R2 of an ANN. 
# Uses leave n out cross validation, where n depends on the size of the dataset
# and is ~20% of the data
  R2 <- 0
  mean.y <- mean(y)
  for(i in unique(gp.ids)) {
    for(try in 1:10) {
      fit <- nnet(x = as.matrix(x[gp.ids != i,par == 1], ncol = sum(par)),
                  y = y[gp.ids != i], size = 1, rang = 0.1, maxit = 1000, 
                  linout = TRUE, trace = FALSE, method = "SANN")
      if(try == 1) {
        fit.best <- fit
      } else if(fit$value < fit.best$value) {
        fit.best <- fit
      }
    }
    pred <- predict(fit.best, 
                    newdata = as.matrix(x[gp.ids == i,par == 1], ncol = sum(par)),
                    type = "raw")
    err <- sum((y[gp.ids == i] - pred)^2)
    R2.i <- 1 - err/sum((y[gp.ids == i] - mean.y)^2)
    R2 <- R2 + R2.i
  }
  R2 <- R2 / length(unique(gp.ids))
  return(R2)
}
#--------------------------------------------------------
names(x) <- names.1
par <- rep(1, ncol(x))

# Use k-fold cross validation to determine stopping criterion (out-of-sample AIC). 
# Number of folds = 5
# fold 'IDs' saved in gp.ids
rep.lnth <- floor(length(y)/5)
gp.ids <- rep(1:5, each = rep.lnth)
if(length(gp.ids) < length(y)) gp.ids <- c(gp.ids, c(1:5)[1:(length(y) - rep.lnth*5)])
gp.ids <- gp.ids[order(runif(length(gp.ids)))]

# Search for best inputs using "GA_Search.R" and objective function 'fitness' above
tmp <- GA.search(par = par, fn = fitness, x = x, y = y, gp.ids = gp.ids,
                 lower = 0, upper = 1, size.pop = 50, max.pop = 10000,
                 prob.cross = 0.7, prob.mut = 0.02, max.eval = 100000,
                 crossover = "single point", mutation = "random", max.conv = 20)

# save best parameters in variable 'combs' (binary string determining which inputs are included - 1 for yes, 0 for no)
# and their corresponding names in 'test.names'
combs <- tmp$par
test.names <- names(x)[which(combs == 1)]
# write 'optimal' inputs to output file
out.file <- paste(out.dir, fn, "_out_AIC.csv", sep = "")
write(paste("Out-of-sample AIC = ", tmp$value, sep = ""), file = out.file, sep = ",")
write.table(test.names, file = out.file, row.names = FALSE, col.names = FALSE,
            sep = ",", quote = FALSE, append = TRUE)
write("", file = out.file, append = TRUE)
out.dat <- as.data.frame(cbind(tmp$counts, tmp$cpu.time, tmp$elap.time))
names(out.dat) <- c("ModEvals","CPUtime", "ElapsedTime")
write.table(out.dat, file = out.file, append = TRUE, sep = ",",
            na = "", row.names = FALSE, col.names = TRUE)


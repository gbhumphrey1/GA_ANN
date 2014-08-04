#  ~ Genetic Algorithm Module ~
#------------------------------------------------------------------------
GA.search <- function(par, fn, ..., lower = -1e6, upper = 1e6,
                      size.pop = 20, max.pop = 10000,
                      prob.cross = 0.7, prob.mut = 0.02,
                      stepsize.0 = 0.2, stepsize.F = 1e-3, elitism = 1,
                      crossover = c("single point", "uniform ave", "child ave"),
                      mutation = c("random", "stepsize"),
                      max.eval = 10000, tol = 1e-8, max.conv = 10,
                      integer.code = TRUE, reset.rand = FALSE, seed = 837648) {

# PURPOSE: Implements the GA search algorithm to find unconstrained maximum
# of user-provided objective (fitness) function.

# Input arguments:
#  function fitness.fn = evaluates the fitness function;
#  lower = lower bound on best.X;
#  upper = upper bound on best.X;
#  inival = optional initial X for search
#  varNames (optional) = variable names
# Output arguments:
#  best.X = variable set at which fitness is at maximum;
#  best.F = best fitness achieved (at bestX);
#
# COMMENTS:
# This routine is not well optimised for very high-dimensional problems:
# 1. local automatic array u is effectively sized nDim**3
# 2. local automatic arrays point, newF and oldF sized nDim**2
# 3. Since the termination criterion does not use the mathematical definition
#    of a maximum of f(x), ie, does not check df/dx=0 and d2f/dx2<0, the
#    success of the optimisation must be verified by the user. the
#    simplest way is to independently calculate the derivaties of the function
#    at the reported GA mode.
#----


#[1] Initialize
#------------
  cpu.time <- proc.time()
  fn1 <- function(par) fn(par, ...)
  time.0 <- Sys.time()
  nDim <- length(par)
  mPop <- size.pop
  elitism.0 <- elitism

  crossover <- match.arg(NULL, crossover)
  mutation <- match.arg(NULL, mutation)

  if(reset.rand) {
    seed <- seed - 153351
    set.seed(seed)
  }

  best.X <- par
  best.F <- -1e8
  old.best.F <- -1e8
  nEval <- 0
  converge <- 0
  converged <- 1
  g <- 0
  stepsize <- stepsize.0
  u <- array(dim = c(mPop, nDim))
  fitness <- rep(-1e-8, mPop)

  j <- 1
  if(all(lower <= par & par <= upper)) { # user-supplied guess ok
    u[j,] <- par
    fitness[j] <- my.fitness(u[j,], fn1) # get fitness
    nEval <- nEval + 1
    if(fitness[j] > -Inf) {
      j <- j + 1
    } else {
      cat("Poor initial solution - solution not used\n")
    }
  } else {
    cat("Infeasible initial solution - solution not used\n")
  }

  not.ok <- j:mPop
  gen <- length(not.ok)
  while(gen > 0) {
    if(nEval > max.eval) { # looks like it's hard to even get a
                           # feasible sample. too hard basket!
      cat("Could not find feasible sample...\n")
      cat("stopped after", max.eval, "iterations\n")
      return(list(par = best.X, value = best.F, counts = nEval,
                  convergence = converged))
    }
    if(integer.code) {
      tmp <- round(runif(nDim*gen, lower, upper))
    } else {
      tmp <- runif(nDim*gen, lower, upper)
    }
    dim(tmp) <- c(gen, nDim)
    u[not.ok,] <- tmp
    fitness[not.ok] <- apply(u[not.ok,], 1, my.fitness, fn1)
    nEval <- nEval + gen
    ok <- fitness > -Inf & !is.na(fitness)
    not.ok <- which(!ok)
    gen <- length(not.ok)
  }

  point <- (1:mPop)[order(fitness, decreasing = TRUE)]
  best.pt <- point[1]
  best.X <- u[best.pt,]
  best.F <- fitness[best.pt]
  old.best.F <- best.F
  df <- 0

  cat("initial ", "value", -best.F, "\n")

# Start of GA loop process
# ---
  g <- 0
  while(converged > 0) {  # main GA loop, iterate until converges or max generations exceeded.

    g <- g + 1  # g represents current generation, loops till max generation reached

  # Elitism ("elitism" determines how many members are replaced)
    if(elitism > 0) {
      for(i in 1:elitism) {
        u[point[mPop - i + 1],] <- best.X
        fitness[point[mPop - i + 1]] <- best.F
      }
    }

  # The GA algorithm evolves the population, hopefully increasing the fitness of chromosomes.
  # On output, u = evolved samples, new.F = corresponding fitnesses.
    tmp <- GA.evolve(fitness, u, fn1, prob.cross, prob.mut, crossover,
                     mutation, stepsize, lower, upper, integer.code)
    fitness <- tmp$fit
    u <- tmp$child
    nEval <- nEval + mPop

    point <- (1:mPop)[order(fitness, decreasing = TRUE)]
    if(fitness[point[1]] > best.F) {
      best.pt <- point[1]
      best.F <- fitness[best.pt]
      best.X <- u[best.pt,]
    }
  # Check for convergence (when best fitness does not improve maxConverge times)
  # save best estimate on exit (even if maxEval exceeded)
    df <- best.F - old.best.F

    if(abs(df) <= max(abs(old.best.F*tol), abs(tol))) {
#      elitism <- 0
      converge <- converge + 1  # looks like the scheme is converging...
      if(converge == 1) convergence <- nEval
      if(mutation == "stepsize" & converge >= max.conv) {
        if(stepsize == stepsize.F) {
          print("GA converged")
          converged <- 0
          cat("final ", "value", -best.F, "\n")
          cat("converged\n")
        } else {
          stepsize <- max(stepsize.F, stepsize*0.9)
          if(stepsize > stepsize.F) converge <- 0
        }
      } else if(converge >= max.conv) {
      # GA converged, terminate main loop
        print("GA converged")
        converged <- 0
        cat("final ", "value", -best.F, "\n")
        cat("converged\n")
      }
    } else if(nEval > max.eval) { # max number of GA iterations exceeded
    # terminate main GA loop. Note: its ok to exceed nEval while scheme
    # is still converging
      cat("final ", "value", -best.F, "\n")
      cat("stopped after", max.eval, "iterations\n")
      break
    } else {           # keep iterating: reset counter convergent iterations...
#      elitism <- elitism.0
      converge <- 0
    }

    old.best.F <- best.F # keep best estimate of location of maximum fitness
    if(g %% 2 == 0) cat("iter ", nEval, "value", -best.F, "\n")

  } # end GA.loop
  cpu.time <- proc.time() - cpu.time
   # ---
  # End procedure here
  return(list(par = best.X, value = best.F, counts = nEval, cpu.time = cpu.time[2],
              elap.time = cpu.time[3], convergence = converged))
}
#-------------------------------------------------------------------------------
GA.evolve <- function(fit, u, fn1, prob.cross, prob.mut, crossover,
                      mutation, stepsize, lower, upper, integer.code) {
# Implements genetic operators selection, crossover and mutation.
# G. B. Kingston

# Start procedure here
# ---
  crossover <- pmatch(crossover, c("single point", "uniform ave", "child ave"))
  mutation <- pmatch(mutation, c("random", "stepsize"))
  mPop <- nrow(u)
  nDim <- ncol(u)
  parent <- array(dim = c(mPop, nDim))
  child <- array(dim = c(mPop, nDim))
#------------------------------------------------------------------------------
#[1] Tournament Selection:
#------------------------------------------------------------------------------
# Call subroutine to randomise the pairings in the pairing array

  pairing <- cbind(sample(1:mPop, mPop), sample(1:mPop, mPop))
# run tournaments and fill mating pool
  s <- pairing[,1]
  t <- pairing[,2]
  sel.s <- fit[s] > fit[t]
  parent[sel.s,] <- u[s[sel.s],]
  parent[!sel.s,] <- u[t[!sel.s],]
# ---
# End tournament selection

#------------------------------------------------------------------------------
#[2] Crossover:
#------------------------------------------------------------------------------
  r.cross <- runif(mPop/2)
  if(crossover == 1) {
    cross <- function(x, prob) {
      n <- length(x)/2
      x.1 <- c.1 <- x[1:n]
      x.2 <- c.2 <- x[-(1:n)]
      r.cross <- runif(1)
      if(r.cross < prob) {
        pt <- round(runif(1, 1, n - 1))
        c.1[(pt + 1):n] <- x.2[(pt + 1):n]
        c.2[(pt + 1):n] <- x.1[(pt + 1):n]
      }
      c <- c(c.1, c.2)
      return(c)
    }
    ids <- rep(1:(mPop/2), each = nDim*2)
    tmp <- matrix(unlist(tapply(c(t(parent)), ids, cross, prob = prob.cross)),
                  byrow = TRUE, ncol = nDim)
    child <- tmp
  } else if(crossover == 2) {
    r.tmp <- runif(mPop/2)
    for(i in seq(2, mPop, by = 2)) {
      if(r.cross[i/2] < prob.cross) {
        if(r.tmp[i/2] < 1/6) {
          child[i - 1,] <- (parent[i - 1,] + parent[i,]) / 2
          child[i,] <- parent[i,]
        } else if(r.tmp[i/2] < 2/6) {
          child[i - 1,] <- (parent[i - 1,] + parent[i,]) / 2
          child[i,] <- parent[i - 1,]
        } else if(r.tmp[i/2] < 3/6) {
          child[i - 1,] = parent[i,]
          child[i,] <- (parent[i - 1,] + parent[i,]) / 2
        } else if(r.tmp[i/2] < 4/6) {
          child[i - 1,] <- parent[i - 1,]
          child[i,] <- (parent[i - 1,] + parent[i,]) / 2
        } else if(r.tmp[i/2] < 5/6) {
          child[i - 1,] <- parent[i,]
          child[i,] <- parent[i - 1,]
        } else {
          child[i - 1,] <- parent[i - 1,]
          child[i,] <- parent[i,]
        }
      } else {
        child[i - 1,] <- parent[i - 1,]
        child[i,] <- parent[i,]
      }
    }
  } else if(crossover == 3) {
    i.tmp <- round(runif(mPop/2, 1, nDim))
    for(i in seq(2, mPop, by = 2)) {
      pt <- i.tmp[i/2]
      if(r.cross[i/2] < prob.cross) {
        child[i - 1,1:pt] <- (parent[i - 1,1:pt] + parent[i,1:pt]) / 2
        child[i,1:pt] <- parent[i,1:pt]
        if(pt < nDim) {
          child[i - 1,(pt + 1):nDim] <- parent[i - 1,(pt + 1):nDim]
          child[i,(pt + 1):nDim] <- (parent[i - 1,(pt + 1):nDim] +
                                     parent[i,(pt + 1):nDim]) / 2
        }
      } else {
        child[i - 1,] <- parent[i - 1,]
        child[i,] <- parent[i,]
      }
    }
  } else {
    print("Invalid crossover method")
  }
# ---
# End crossover

#------------------------------------------------------------------------------
#[3] Mutation
#------------------------------------------------------------------------------
  if(mutation == 1) {
  # Random Replacement Mutation
    mut <- function(x, prob, low, upp) {
      n <- length(x)
      r.mut <- runif(n)
      sel <- r.mut < prob
      if(integer.code) {
        x[sel] <- round(runif(sum(sel), low, upp))
      } else {
        x[sel] <- runif(sum(sel), low, upp)      
      }
      return(x)
    }
    tmp <- t(apply(child, 1, mut, prob = prob.mut, low = lower, upp = upper))
    child <- tmp
  } else if(mutation == 2) {
  # Stepsize Mutation
    mut <- function(x, prob, step, upp, low) {
      n <- length(x)
      r.mut <- runif(n)
      sel <- r.mut < prob
      mutate <- runif(sum(sel), -1, 1)
      mutate <- mutate * step
      rev <- (x[sel] + mutate) > upp | (x[sel] + mutate) < low
      mutate[rev] <- -mutate[rev]
      x[sel] <- x[sel] + mutate
#      x[sel] <- pmin(upp, x[sel] + mutate)
#      x[sel] <- pmax(low, x[sel])
      return(x)
    }
    tmp <- t(apply(child, 1, mut, prob = prob.mut, step = stepsize, upp = upper,
                   low = lower))
    child <- tmp
  } else {
    print("Invalid mutation method")
  }
# ---
# End mutation

#------------------------------------------------------------------------------
#[4] Evaluate fitness of each member in new population
#------------------------------------------------------------------------------

# Generate mPop samples (parents) in the search hyperspace
  fit <- apply(child, 1, my.fitness, fn1) # get fitness  
  return(list(fit = fit, child = child))
# ---
# End procedure here
}
#-------------------------------------------------------------------------------
my.fitness <- function(x, fn1) {
# Purpose: Allows more effective runtime diagnostics

  fx <- fn1(x)
  return(fx)
}
#-------------------------------------------------------------------------------



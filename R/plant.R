## Support functions for growing an individual to a given size, etc.  Used
## in dfalster/TraitGrowthTrajectories

##' Grow an individual up to particular sizes.
##'
##' @title Grow individual to given size
##' @param individual An \code{Individual} object.
##' @param sizes A vector of sizes to grow the plant to (increasing in
##' size).
##' @param size_name The name of the size variable within
##' \code{individual$rates} (e.g., height).
##' @param env An \code{Environment} object.
##' @param time_max Time to run the ODE out for -- only exists to
##' prevent an infinite loop (say, on an unreachable size).
##' @param warn Warn if requesting a plant that is too large?
##' @param filter Filter individuals that are too large?
##' @return A list with elements \code{time} (the time that a given
##' size was reached), \code{state} (the \emph{ode state} at these
##' times, as a matrix) and \code{plant} a list of individuals grown to the
##' appropriate size.  Note that if only a single size is given,
##' a list of length 1 is returned.
##' @export
grow_individual_to_size <- function(individual, sizes, size_name, env,
                               time_max=Inf, warn=TRUE, filter=FALSE) {
  obj <- grow_individual_bracket(individual, sizes, size_name, env, time_max, warn)

  polish <- function(i) {
    grow_individual_bisect(obj$runner, sizes[[i]], size_name,
                      obj$t0[[i]], obj$t1[[i]], obj$y0[i, ])
  }
  res <- lapply(seq_along(sizes), polish)

  state <- t(sapply(res, "[[", "state"))
  colnames(state) <- colnames(obj$state)

  ret <- list(time=vnapply(res, "[[", "time"),
              state=state,
              individual=lapply(res, "[[", "individual"),
              trajectory=cbind(time=obj$time, state=obj$state),
              env=env)
  if (filter) {
    i <- !vlapply(ret$individual, is.null)
    if (!all(i)) {
      ret$time  <- ret$time[i]
      ret$state <- ret$state[i, , drop=FALSE]
      ret$individual <- ret$individual[i]
    }
  }
  ret
}

##' @export
##' @rdname grow_individual_to_size
##' @param ... Additional parameters passed to
##' \code{grow_individual_to_size}.
##' @param heights Heights (when using \code{grow_individual_to_height})
grow_individual_to_height <- function(individual, heights, env, ...) {
  grow_individual_to_size(individual, heights, "height", env, ...)
}

##' Grow a plant up for particular time lengths
##'
##' @title Grow a plant
##' @param individual An \code{Individual} object
##' @param times A vector of times
##' @param env An \code{Environment} object
##' @export
grow_individual_to_time <- function(individual, times, env) {
  if (any(times < 0.0)) {
    stop("Times must be positive")
  }
  n <- length(times)
  if (n == 0L) {
    stop("At least one time must be given")
  }

  y0 <- individual$ode_state
  t0 <- 0.0
  i <- 1L
  t_next <- times[[i]]
  strategy_name <- individual$strategy_name

  ir1 <- IndividualRunner(strategy_name, environment_type(strategy_name))(individual, env)
  ir2 <- IndividualRunner(strategy_name, environment_type(strategy_name))(individual, env)

  runner <- OdeRunner(strategy_name)(ir1)
  runner_detail <- OdeRunner(strategy_name)(ir2)

  ## TODO: This could also be done by better configuring the
  ## underlying ODE runner, but this seems a reasonable way of getting
  ## things run for now.
  state <- matrix(NA, n, length(y0))
  colnames(state) <- individual$ode_names

  individual <- vector("list", n)
  while (i <= n) {
    runner$step()
    t1 <- runner$time
    y1 <- runner$state
    while (t_next < t1 && i <= n) {
      runner_detail$set_state(y0, t0)
      runner_detail$step_to(t_next)
      state[i, ] <- runner_detail$state
      individual[[i]] <- runner_detail$object$individual
      i <- i + 1L
      t_next <- times[i] # allows out-of-bounds extraction
    }
    t0 <- t1
    y0 <- y1
  }

  list(time=times, state=state, individual=individual, env=env)
}

## internal funciton to grab the internal state of the ode runner:
get_individual_internals_fun <- function (individual) {
  get(paste0(individual$strategy_name, '_oderunner_individual_internals')) #!
}

grow_individual_bracket <- function(individual, sizes, size_name, env,
                               time_max=Inf, warn=TRUE) {
  if (length(sizes) == 0L || is.unsorted(sizes)) {
    stop("sizes must be non-empty and sorted")
  }
  if (individual$state(size_name) > sizes[[1]]) {
    stop("Individual already bigger than smallest target size")
  }
  strategy_name <- individual$strategy_name

  # TODO: size index uses index from 0
  # can we clarify?
  size_index <- (which(individual$ode_names == size_name) - 1)

  runner <- OdeRunner(strategy_name)(IndividualRunner(strategy_name, environment_type(strategy_name))(individual, env))
  internals <- get_individual_internals_fun(runner$object$individual)
  i <- 1L
  n <- length(sizes)
  j <- rep_len(NA_integer_, n)
  state <- list(list(time=runner$time, state=runner$state))

  while (i <= n & runner$time < time_max) {
    ok <- tryCatch({
      runner$step()
      TRUE
    },
    error=function(e) {
      msg <- c("Stopping early as integration failed with error: ", e$message,
               sprintf("  %d larger sizes dropped", sum(is.na(j))))
      if (warn) {
        warning(paste(msg, collapse="\n"), immediate.=TRUE)
      }
      ## TODO: Consider making this an error, or making the test a bit better.
      if (runner$object$individual$ode_rates[[2]] < 1e-10) {
        warning("Integration may have failed for reasons other than mortality",
                immediate.=TRUE)
      }
      FALSE
    })
    if (!ok) {
      break
    }
    state <- c(state, list(list(time=runner$time, state=runner$state)))


    while (i <= n && internals(runner)$state(size_index) > sizes[[i]]) {
      j[[i]] <- length(state) - 1L
      i <- i + 1L
    }
    if (runner$time >= time_max) {
      if (warn) {
        warning(sprintf("Time exceeded time_max, %d larger sizes dropped",
                        sum(is.na(j))), immediate.=TRUE)
      }
      break
    }
  }

  t <- vnapply(state, "[[", "time")
  m <- t(sapply(state, "[[", "state"))
  k <- j + 1L
  colnames(m) <- runner$object$individual$ode_names
  list(t0=t[j],
       t1=t[k],
       y0=m[j, , drop=FALSE],
       y1=m[k, , drop=FALSE],
       time=t,
       state=m,
       index=j,
       runner=runner)
}

##' @noRd
##' @importFrom stats uniroot
grow_individual_bisect <- function(runner, size, size_name, t0, t1, y0) {

  
  # TODO: size index uses index from 0
  # can we clarify?
  size_index <- (which(runner$object$individual$ode_names == size_name) - 1)

  internals <- get_individual_internals_fun(runner$object$individual)
  f <- function(t1) {
    runner$set_state(y0, t0)
    runner$step_to(t1)
    internals(runner)$state(size_index) - size
  }

  if (is.na(t0) || is.na(t1) || any(is.na(y0))) {
    y0[] <- NA_real_
    list(time=NA_real_, state=y0, individual=NULL)
  } else {
    root <- uniroot(f, lower=t0, upper=t1)
    list(time=root$root, state=runner$state, individual=runner$object$individual)
  }
}

#!
#' Compute the whole plant light compensation point for a single
#' plant.
#' @title Whole plant light compensation point
#' @param p A \code{PlantPlus}, with strategy, height, etc set.
#' @param ... Additional arguments that are ignored
#' @export
#' @author Rich FitzJohn
resource_compensation_point <- function(p, ...) {
  UseMethod("resource_compensation_point")
}

##' @export
resource_compensation_point.Plant <- function(p, ...) {
  resource_compensation_point(p, ...)
}

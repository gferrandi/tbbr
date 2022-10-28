# norm2 --------
#' L2 norm of a vector
#'
#' @param x vector
#'
#' @return L2 norm of a vector
norm2 <- function(x) {sqrt(sum(x^2))}

# angle -------------------------------------------------------------------
#' Angle between two vectors
#'
#' @param a first vector
#' @param b second vector
#'
#' @return angle between the two vectors in radiants
angle <- function(a,b){
  a <- a/norm2(a); b <- b/norm2(b)
  asin(min(1, norm2(a- sum(a*b)*b)))
}

#' Gallery of targets
#'
#' @param target_option string corresponding to one of the available targets (add which)
#'
#' @details e.g.,
#' \tabular{ll}{
#' Target \tab Description \cr 
#' cotan1_1 \tab cotangent step with $q=r=1$ \cr
#' cotan1_2 \tab cotangent step with $q=1$ and $r=2$ \cr
#' cotan2_1 \tab cotangent step with $q=2$ and $r=1$ \cr
#' cotan1_H \tab cotangent step with $q=1$ and $r=0.5$ \cr
#' cotanH_1 \tab cotangent step with $q=0.5$ and $r=1$ \cr
#' inv_bb2_2_01 \tab positive target $2.01$*INV BB2 \cr
#' inv_bb2_100 \tab positive target $100$*INV BB2 \cr
#' iter \tab positive target $k$*INV BB2 k \cr
#' }
#' 
#' @md
#'
#' @return function with three arguments: angle between $s_k$ and $y_k$, inverse BB2 $\|y_k\|^2/s_k^Ty_k$, number of iterations
#' @export
#'
#' @examples
#' target_gallery("cotan1_1")
target_gallery <- function(target_option){
  switch(target_option,
         "bb1"           = function (gy, y, sng, ng, alpha, iter) -alpha*sng/gy,
         "bb2"           = function (gy, y, sng, ng, alpha, iter) {yy <- sum(y*y); -alpha*gy/yy},
         "homo_abb"       = function (gy, y, sng, ng, alpha, iter) {
           yy <- sum(y*y); cosa <- abs(gy)/(ng*sqrt(yy))
           if(cosa^2 < 0.8) -0.5*alpha*(sng - yy/alpha^2 + sqrt((sng - yy/alpha^2)^2 + 4*(gy/alpha)^2))/gy else -alpha*sng/gy
         },
         "homo_bb"       = function (gy, y, sng, ng, alpha, iter) {
           yy <- sum(y*y)
           -0.5*alpha*(sng - yy/alpha^2 + sqrt((sng - yy/alpha^2)^2 + 4*(gy/alpha)^2))/gy
         },
         "cotan1_1"      = function (gy, y, sng, ng, alpha, iter) {
           sy <- -alpha*gy; yy <- sum(y*y); cosa <- abs(gy)/(ng*sqrt(yy)); sina <- sqrt(1-min(cosa^2,1))
           tau <- (cosa); tau0 <- (sina); (tau0*sy+tau*alpha^2*sng)/(tau0*yy+tau*sy)
         },
         "cotan2_1"      = function (gy, y, sng, ng, alpha, iter) {
           sy <- -alpha*gy; yy <- sum(y*y); cosa <- abs(gy)/(ng*sqrt(yy)); sina <- sqrt(1-min(cosa^2,1))
           tau <- (cosa)^2; tau0 <- (sina); (tau0*sy+tau*alpha^2*sng)/(tau0*yy+tau*sy)
         },
         "cotan1_2"      = function (gy, y, sng, ng, alpha, iter) {
           sy <- -alpha*gy; yy <- sum(y*y); cosa <- abs(gy)/(ng*sqrt(yy)); sina <- sqrt(1-min(cosa^2,1))
           tau <- (cosa); tau0 <- (sina)^2; (tau0*sy+tau*alpha^2*sng)/(tau0*yy+tau*sy)
         },
         "cotan1_H"      = function (gy, y, sng, ng, alpha, iter) {
           sy <- -alpha*gy; yy <- sum(y*y); cosa <- abs(gy)/(ng*sqrt(yy)); sina <- sqrt(1-min(cosa^2,1))
           tau <- (cosa); tau0 <- sqrt(sina); (tau0*sy+tau*alpha^2*sng)/(tau0*yy+tau*sy)
         },
         "cotanH_1"      = function (gy, y, sng, ng, alpha, iter) {
           sy <- -alpha*gy; yy <- sum(y*y); cosa <- abs(gy)/(ng*sqrt(yy)); sina <- sqrt(1-min(cosa^2,1))
           tau <- sqrt(cosa); tau0 <- (sina); (tau0*sy+tau*alpha^2*sng)/(tau0*yy+tau*sy)
         },
         "inv_bb2_2_01"  = function (gy, y, sng, ng, alpha, iter) {
           yy <- sum(y*y); k <- 2.01; alpha*(gy/yy-k*sng/gy)/(k-1)
         },
         "inv_bb2_100"   = function (gy, y, sng, ng, alpha, iter) {
           yy <- sum(y*y); k <- 100; alpha*(gy/yy-k*sng/gy)/(k-1)
         },
         "iter"          = function (gy, y, sng, ng, alpha, iter) {
           yy <- sum(y*y); 
           if(iter == 0) {-alpha*gy/yy} else {k <- iter+1; alpha*(gy/yy-k*sng/gy)/(k-1)}
         }
  )
}

# parameters of sg --------------------------------------------------------
#' Parameters of global TBB
#'
#' @param alpha_min min value for the stepsize, default = 1e-30
#' @param alpha_max max value for the stepsize, default = 1e30
#' @param alpha initial value for the stepsize, default = 1
#' @param gamma line search param, default = 1e-4
#' @param M number of stored previous function evaluations, default = 10
#' @param tol tolerance for the stopping criterion, default = 1e-6
#' @param iter_max maximum number of iterations, default = 1e4
#' @param iter_int_max maximum number of iterations for the backtracking, default = 100
#' @param target_option string or function. See target_gallery for more details; default = "bb1"
#' @param ... additional arguments (not needed at the moment)
#'
#' @details Other stepsizes:
#' \tabular{ll}{
#' Target \tab Description \cr 
#' bb1 \tab Barzilai -- Borwein 1 \cr
#' bb2 \tab Barzilai -- Borwein 2 \cr
#' abb \tab Adaptive BB with threshold 0.8 \cr
#' homo_bb \tab Homogeneous BB \cr
#' homo_abb \tab Adaptive homo BB with threshold 0.8 \cr
#' }
#' 
#' @md
#'
#' @return list with parameters for the 
#' @export
#'
#' @examples
#' gtbb_control()
gtbb_control <- function(alpha_min = 1e-30, alpha_max = 1e30, alpha = NULL, 
                               gamma = 1e-4, M = 10, M2 = 5, tol = 1e-6, iter_max = 5e4, iter_int_max = 100, 
                               target_option = "bb1", negstep = "normg", t = 0.8, ...){
  args <- list(...)
  c(list(alpha_min = alpha_min, alpha_max = alpha_max, alpha = alpha, 
         gamma = gamma, M = M, M2 = M2, tol = tol, iter_max = iter_max, iter_int_max = iter_int_max, 
         target_option = target_option, negstep = negstep, t = t), args)
}

# unconstrained spectral gradient + nonmonotone ls ------------------------
#' Global TBB 
#' 
#' Gradient method with nonmonotone line search for continuously differentiable functions. 
#' 
#' It is endowed with the TBB stepsize, a Barzilai--Borwein like stepsize with an additional parameter, called target.
#'
#' @param x starting vector
#' @param obj objective function. Must be of the form obj(x)
#' @param grad gradient function. Must be of the form grad(x)
#' @param p control parameters. See gtbb_control for more details
#' @param quadratic whether the objective function is quadratic or not, default = FALSE
#' @param reltol whether the tolerance is multiplied by the norm of the first gradient (TRUE) or not. default = TRUE
#' @param info save information of the whole method, default = FALSE
#'
#' @return minimizer and some information: number of function evaluations, number of iterations, last gradient norm, tolerance, final objective value.
#' It optionally returns the aforementioned information about all the iterations.
#' @export
#'
#' @examples
#' A <- c(1:100)
#' x <- rep(1:100)
#' 
#' obj <- function(x) sum(A * x^2)/2
#' grad <- function(x) A * x
#' 
#' gtbb(x, obj, grad, quadratic = TRUE)
gtbb <- function(x, obj, grad, p = gtbb_control(), quadratic = FALSE, reltol = TRUE, info = FALSE){
  
  ctime <- Sys.time()
  
  # get params
  p <- do.call(gtbb_control, as.list(p))
  # p$M <- 10 # fix M
  # p$alpha <- 1 # fix alpha
  # p$iter_max <- 5e4 # fix iter max
  
  # define target function
  if(!(p$target_option %in% c("bon", "fra1", "abb"))){
    if(!is.function(p$target_option)){
      tau_fun <- target_gallery(p$target_option)
    } else {
      tau_fun <- p$target_option
      p$target_option <- "new target"
    }
  }
  
  # init params for ABB 
  if(p$target_option == "bon") p$t <- 0.5
  if(p$target_option == "abb") p$M2 <- 0
  mem <- rep(Inf, p$M2) 
  
  # init f, gradient, alpha
  f0 <- f <- obj(x); g <- grad(x)
  alpha <- if(is.null(p$alpha)) max(p$alpha_min, min(1/norm2(g), p$alpha_max)) else p$alpha
  x0norm <- sum(abs(x))
  alpha0 <- alpha
  
  # init objects to store previous iterations
  if(!quadratic) f_hist <- rep(f, p$M)
  
  # save for reproducibility
  x0 <- x
  
  # tracking parameters
  if(info) alg_log <- NULL; 
  nfe <- 1; const_step <- 0; gy <- NA
  
  # relative tolerance
  if(reltol) p$tol <- p$tol*norm2(g)
  
  tryCatch({
    for(iter in 0:(p$iter_max-1)){
      
      # stopping criterion
      sng <- sum(g^2)
      ng <- sqrt(sng)
      
      # save info 
      if(info){
        iter_info <- c(iter = iter, stop_criterion = ng, obj = f, alpha = alpha, nfe = nfe, gy = gy)
        alg_log <- rbind(alg_log, iter_info)
      }
      
      # is current point stationary? 
      if(ng < p$tol) break
      
      # BACKTRACKING 
      if(!quadratic){
        xnew <- x - alpha*g; f <- obj(xnew); nfe <- nfe + 1
        fref <- max(f_hist, na.rm = TRUE)
        while(f > fref - 1e-4*alpha*sng){
          alpha <- alpha/2
          stopifnot(alpha > p$alpha_min)
          xnew <- x - alpha*g; f <- obj(xnew); nfe <- nfe + 1; 
        }
        x <- xnew; f_hist <- c(f_hist, f); f_hist <- f_hist[-1] 
      } else {
        x <- x -alpha*g; f <- NA
      }
      
      # update gradients
      g0 <- g; g <- grad(x); y <- g - g0
      
      # update stepsize 
      gy <- sum(g0*y) 
      
      if(gy < 0){ # -gy > 0
        if(p$target_option %in% c("bon", "fra1", "abb")){
          yy <- sum(y*y); tmp <- -alpha*gy/yy
          #if(!quadratic & tmp < 0) tmp <- 1/max(1e-5, min(norm2(g), 1)) 
          mem <- c(mem, tmp)
          cosa <- abs(gy)/(ng*sqrt(yy))
          if(cosa^2 < p$t) {
            alpha <- min(mem)
            if(p$target_option == "bon") p$t <- 0.9*p$t
          } else {
            alpha <- -alpha*sng/gy
            if(p$target_option == "bon") p$t <- 1.1*p$t
          }
          mem <- mem[-1]
        } else {
          alpha <- tau_fun(gy, y, sng, ng, alpha, iter)
        }
        #alpha0 <- alpha
      } else {
        alpha <- 1/max(1e-5, min(norm2(g), 1)); const_step <- const_step + 1
        #alpha <-alpha0; const_step <- const_step + 1
      }
      
      
      if(!quadratic){
        alpha <- max(p$alpha_min, min(alpha, p$alpha_max)) 
      }
    }
  },
  error = function(e){
    cat("ERROR:", conditionMessage(e), "\n")
    return(NA)
  }
  )
  
  # prepare output
  last_iter <- data.frame(id_algo = p$target_option, iter = iter, nfe = nfe, stop_criterion = ng, fstar = obj(x))
  last_iter$const_step <- const_step
  last_iter$ctime <- difftime(Sys.time(), ctime, units = "secs") 
  last_iter$tol <- p$tol
  last_iter$f0 <- f0
  last_iter$x0 <- x0norm
  last_iter$xstar <- sum(abs(x))
  
  if(info){
    alg_log <- as.data.frame(alg_log); rownames(alg_log) <- NULL
    return(list(alg_log = alg_log, last_iter = last_iter, x = x))
  } else {
    return(list(last_iter = last_iter, x = x))
  }
}

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
  
  implemented <- c("cotan1_1", "cotan2_1", "cotan1_2", "cotan1_H", "cotanH_1", 
                   "inv_bb2_2_01", "inv_bb2_100", "iter", "bb1", "bb2", "abb", "homo_bb", "homo_abb")
  
  if(!(target_option %in% implemented)) {
    stop(cat("Target option not implemented. Try\n", paste0(implemented, collapse = "\n")))
  }
  
  switch(target_option,
         "cotan1_1"      = function (a,b,k) -cos(a)/sin(a),
         "cotan2_1"      = function (a,b,k) -(cos(a))^2/sin(a),
         "cotan1_2"      = function (a,b,k) -cos(a)/(sin(a))^2,
         "cotan1_H"      = function (a,b,k) -cos(a)/sqrt(sin(a)),
         "cotanH_1"      = function (a,b,k) -sqrt(cos(a))/sin(a),
         "inv_bb2_2_01"  = function (a,b,k) 2.01*b,
         "inv_bb2_100"   = function (a,b,k) 100*b,
         "iter"          = function (a,b,k) {if(k == 0) 0 else (k+1)*b})
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
gtbb_control <- function(alpha_min = 1e-8, alpha_max = 1e30, alpha = NULL, 
                               gamma = 1e-4, M = 5, tol = 1e-6, iter_max = 1e4, iter_int_max = 30, 
                               target_option = "bb1", ...){
  args <- list(...)
  c(list(alpha_min = alpha_min, alpha_max = alpha_max, alpha = alpha, 
         gamma = gamma, M = M, tol = tol, iter_max = iter_max, iter_int_max = iter_int_max, 
         target_option = target_option), args)
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
  
  # get params
  p <- do.call(gtbb_control, as.list(p))
  
  # define target function
  if(!is.function(p$target_option)){
    tau_fun <- target_gallery(p$target_option)
  } else {
    tau_fun <- p$target_option
    p$target_option <- "new target"
  }
  
  # init f, gradient, alpha
  f <- obj(x); g <- grad(x)
  alpha <- if(is.null(p$alpha)) max(p$alpha_min, min(1/norm2(g), p$alpha_max)) else p$alpha
  alpha0 <- alpha
  
  # init objects to store previous iterations
  if(!quadratic){f_hist <- rep(f, p$M)}
  
  # save for reproducibility
  x0 <- x
  
  # tracking parameters
  if(info) alg_log <- NULL; 
  nfe <- 1; const_step <- 0
  
  # relative tolerance
  if(reltol) p$tol <- p$tol*norm2(g)
  
  tryCatch({
    for(iter in 0:(p$iter_max-1)){
      
      # stopping criterion
      sng <- sum(g^2)
      ng <- sqrt(sng)
      
      # save info 
      if(info){
        iter_info <- c(iter = iter, stop_criterion = ng, obj = f, alpha = alpha, nfe = nfe)
        alg_log <- rbind(alg_log, iter_info)
      }
      
      # is current point stationary? 
      if(ng < p$tol) break
      
      # START BACKTRACKING 
      if(!quadratic){
        for (k in 1:p$iter_int_max) {
          xnew <- x - alpha*g; f <- obj(xnew); nfe <- nfe + 1
          if (f <= max(f_hist, na.rm = TRUE)-p$gamma*alpha*sng) break else alpha <- alpha/2
        }
        x <- xnew
        f_hist<- if(p$M == 1) f else c(f_hist[2:p$M], f) # save previous function value
        # END BACKTRACKING 
      } else {
        x <- x -alpha*g; f <- obj(x); nfe <- nfe + 1
      }
      
      # update gradients
      g0 <- g; g <- grad(x); y <- g - g0
      
      # update stepsize 
      p0y <- sum(g0*y); yy <- sum(y*y); sy <- -alpha*p0y
      
      if(p$target_option == "bb1"){
        alpha <- -alpha*sng/p0y
      } else if(p$target_option == "bb2"){
        alpha <- sy/yy
      } else if(p$target_option == "abb"){
        if(cos(angle(g0, y))^2 < 0.8) alpha <- sy/yy else alpha <- -alpha*sng/p0y
      } else if(p$target_option == "homo_bb"){
        mu <- 0.5*(sng - yy/alpha^2 + sqrt((sng - yy/alpha^2)^2 + 4*(p0y/alpha)^2))
        alpha <- -alpha*mu/p0y
      } else if(p$target_option == "homo_abb"){
        if(cos(angle(g0, y))^2 < 0.8){
          mu <- 0.5*(sng - yy/alpha^2 + sqrt((sng - yy/alpha^2)^2 + 4*(p0y/alpha)^2))
          alpha <- -alpha*mu/p0y
        } else {
          alpha <- -alpha*sng/p0y
        }
      } else {
        tauk <- tau_fun(angle(g0, y),yy/sy,iter)
        alpha <- (sy-tauk*alpha^2*sng)/(yy-tauk*sy)  
      }
      
      if(!quadratic){
        if(alpha > 0) {alpha0 <- alpha} else {alpha <- alpha0; const_step <- const_step + 1}
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
  last_iter <- data.frame(id_algo = p$target_option, iter = iter, nfe = nfe, stop_criterion = ng, obj = f)
  last_iter$const_step <- const_step
  last_iter$tol <- p$tol
  
  if(info){
    alg_log <- as.data.frame(alg_log); rownames(alg_log) <- NULL
    return(list(alg_log = alg_log, last_iter = last_iter, x = x))
  } else {
    return(list(last_iter = last_iter, x = x))
  }
}
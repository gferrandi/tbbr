#' Unconstrained problems
#'
#' Gallery of unconstrained optimization problems from different sources. They are all based on sum of functions, where the number of terms is related to the size of the problem.
#'
#' @param exp string that identifies the experiment in the gallery
#' @param n dimension of the problem 
#'
#' @details e.g.,
#' \tabular{lll}{
#' Name \tab Function \tab Reference \cr
#' 1 \tab Strictly Convex 1 \tab Raydan 1997 \cr
#' 2 \tab Strictly convex 2 \tab Raydan 1997 \cr
#' 3 \tab Brown almost linear \tab Raydan 1997 \cr
#' 4 \tab Trigonometric \tab Raydan 1997 \cr
#' 5 \tab Broyden tridiagonal \tab Raydan 1997 \cr
#' 7 \tab Extended Rosenbrock \tab Raydan 1997 \cr
#' 8 \tab Penalty 1 \tab Raydan 1997 \cr
#' 10 \tab Variably dimensioned \tab Raydan 1997 \cr
#' 11 \tab Extended Powell \tab Raydan 1997 \cr
#' 12 \tab Generalized Rosenbrock \tab Raydan 1997 \cr
#' hager \tab Hager \tab Andrei 2008 \cr
#' gentridiag1 \tab Generalized tridiagonal 1 \tab Andrei 2008 \cr
#' genWH \tab Generalized White and Holst \tab Andrei 2008 \cr
#' diag5 \tab Diagonal 5 \tab Andrei 2008 \cr
#' }
#' 
#' @md
#'
#' @return a list with suggested starting point, objective and gradient functions
#' @export
#'
#' @examples
#' function_gallery("genWH", 10)
function_gallery <- function(exp, n){
  define_f <- switch(exp,
                     "3" = function(x, n){
                       f          <- rep(0, n)
                       f[1:(n-1)] <- x[1:(n-1)] - (n+1) + sum(x)
                       f[n]       <- prod(x) - 1
                       f
                     },
                     "4" = function(x, n){
                       n - sum(cos(x)) + c(1:n)*(1-cos(x)) - sin(x)
                     },
                     "5" = function(x, n){
                       (3-2*x)*x - c(0, x[1:(n-1)]) - 2*c(x[2:n], 0) + 1
                     },
                     "7" = function(x, n){ # extended rosenbrock (easy version)
                       odd <- as.logical(1:n %% 2)
                       f <- rep(0, n)
                       f[odd] <- 10*(x[!odd] - x[odd]^2)
                       f[!odd] <- 1 - x[odd]
                       f
                     },
                     "8" = function(x, n){
                       f <- sqrt(1e-5)*(x - 1)
                       c(f, sum(x^2) - 1/4)
                     },
                     "10" = function(x, n){
                       f      <- rep(0, n)
                       f[1:n] <- x[1:n]-1
                       f[n+1] <- sum(c(1:n)*f[1:n])
                       f[n+2] <- f[n+1]^2
                       f
                     },
                     "11" = function(x, n){
                       f <- rep(0,n); i <- 1:n %% 4
                       f[i == 0] <- sqrt(10)*(x[i == 1] - x[i == 0])^2
                       f[i == 3] <- (x[i == 2] - 2*x[i ==3])^2
                       f[i == 2] <- sqrt(5)*(x[i == 3] - x[i == 0])
                       f[i == 1] <- x[i == 1] + 10*x[i == 2]
                       f
                     }
  )
  
  exp_ls <- switch(exp,
                   "1" = list(x0    = c(1:n)/n,
                              obj   = function(x) {sum(exp(x) - x)},
                              grad  = function(x) {exp(x) - 1}),
                   "2" = list(x0    = rep(1,n),
                              obj   = function(x) {sum(c(1:n)*(exp(x) - x))/10},
                              grad  = function(x) {c(1:n)*(exp(x) - 1)/10}),
                   "3" = list(x0    = rep(0.5, n),
                              obj   = function(x) {sum(define_f(x, n)^2)/2},
                              grad  = function(x) {f <- define_f(x, n)
                              out <- c(f[1:(n-1)],0) + sum(f[1:(n-1)])
                              tmp <- numeric(); for(i in 1:n) {tmp[i] <- prod(x[-i])}
                              out + f[n]*tmp
                              }),
                   "4" = list(x0    = rep(1,n)/n,
                              obj   = function(x) {sum(define_f(x, n)^2)/2},
                              grad  = function(x) {f <- define_f(x, n); sum(f)*sin(x) + f*(c(1:n)*sin(x) - cos(x))}),
                   "5" = list(x0    = rep(-1,n),
                              obj   = function(x) {sum(define_f(x, n)^2)/2},
                              grad  = function(x) {f <- define_f(x, n); f*(3-4*x) - c(f[2:n], 0) - 2*c(0, f[1:(n-1)])}),
                   "7" = list(x0    = rep(c(-1.2, 1), floor(n/2)),
                              obj   = function(x) {sum(define_f(x, n)^2)/2},
                              grad  = function(x) {f <- define_f(x, n); odd <- as.logical(1:n %% 2); g <- rep(0,n)
                              g[odd] <- -20*f[odd]*x[odd] - f[!odd]; g[!odd] <- 10*f[odd]; g
                              }),
                   "8" = list(x0    = c(1:n),
                              obj   = function(x) {sum(define_f(x, n)^2)/2},
                              grad  = function(x) {f <- define_f(x, n); sqrt(1e-5)*f[1:n] + 2*f[n+1]*x}),
                   "10" = list(x0   = 1 - c(1:n)/n,
                               obj  = function(x) {sum(define_f(x, n)^2)/2},
                               grad = function(x) {f <- define_f(x, n); f[1:n] + c(1:n)*f[n+1]*(1 + 2*f[n+2])}),
                   "11" = list(x0   = rep(c(3,-1,0,1), floor(n/4)),
                               obj  = function(x) {sum(define_f(x, n)^2)/2},
                               grad = function(x) {f <- define_f(x, n); i <- 1:n %% 4; g <- rep(0,n)
                               g[i == 0] <- -2*sqrt(10)*f[i == 0]*(x[i == 1] - x[i == 0]) - sqrt(5)*f[i == 2]
                               g[i == 3] <- -4*f[i == 3]*(x[i == 2] - 2*x[i == 3]) + sqrt(5)*f[i == 2]
                               g[i == 2] <- 2*f[i == 3]*(x[i == 2] - 2*x[i == 3]) + 10*f[i == 1]
                               g[i == 1] <- 2*sqrt(10)*f[i == 0]*(x[i == 1] - x[i == 0]) + f[i == 1]
                               g
                               }),
                   "12" = list(x0    = rep(c(-1.2, 1), floor(n/2)), # n must be even (?)
                               obj   = function(x) {sum(100*(x[2:n] - x[1:(n-1)]^2)^2 + (1-x[1:(n-1)])^2)/2},
                               grad  = function(x) {g <- c(x[1:(n-1)]-1,0) 
                               g[1] <- g[1] -200*(x[2] - x[1]^2)*x[1]
                               g[2:(n-1)] <- g[2:(n-1)] - 200*(x[3:n] - x[2:(n-1)]^2)*x[2:(n-1)] +100*(x[2:(n-1)] - x[1:(n-2)]^2)
                               g[n] <- g[n] + 100*(x[n] - x[n-1]^2)
                               g
                               }),
                   "diag5" = list(x0 = rep(1,n), # too easy
                                  obj = function(x){sum(log(exp(x)+exp(-x)))},
                                  grad = function(x){(exp(x)-exp(-x))/(exp(x)+exp(-x))}),
                   "hager" = list(x0 = rep(1,n),
                                  obj = function(x){sum(exp(x) - sqrt(c(1:n))*x)},
                                  grad = function(x){exp(x) - sqrt(c(1:n))}),
                   "gentridiag1" = list(x0 = rep(2, n),
                                        obj = function(x){sum((x[2:n] + x[1:(n-1)] - 3)^2 + (x[1:(n-1)] -x[2:n] + 1)^4)/2},
                                        grad = function(x){c(x[2:n] + x[1:(n-1)] - 3, 0) + c(0, x[2:n] + x[1:(n-1)] - 3) +
                                            2*c(-x[2:n] + x[1:(n-1)] + 1, 0)^3 - 2*c(0, -x[2:n] + x[1:(n-1)] + 1)^3}),
                   "genWH" = list(x0    = rep(c(-1.2, 1), floor(n/2)), # n must be even (?)
                                  obj   = function(x) {sum(100*(x[2:n] - x[1:(n-1)]^3)^2 + (1-x[1:(n-1)])^2)/2},
                                  grad  = function(x) {g <- c(x[1:(n-1)]-1,0) 
                                  g[1] <- g[1] -300*(x[2] - x[1]^3)*(x[1]^2)
                                  g[2:(n-1)] <- g[2:(n-1)] - 300*(x[3:n] - x[2:(n-1)]^3)*(x[2:(n-1)]^2) +100*(x[2:(n-1)] - x[1:(n-2)]^3)
                                  g[n] <- g[n] + 100*(x[n] - x[n-1]^3)
                                  g
                                  })
  )
  exp_ls
}
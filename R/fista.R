#' Solve SLOPE with FISTA
#'
#' @param x design matrix
#' @param y response vector
#' @param lambda regularization
#' @param maxit max iterations
#' @param opt optimal value
#' @param opt_tol relative tolerance for optimal value
#'
#' @return A list of time, beta, and loss
#' @export
fista <- function(x, y, lambda, maxit = 100, opt = -Inf, opt_tol = 1e-3) {

  p <- ncol(x)

  learning_rate <- 1
  eta <- 0.5

  beta_old <- beta_tilde <- beta <- double(p)

  time <- double(maxit)
  loss <- double(maxit)

  start_time <- proc.time()[3]

  t <- 0
  t_old <- 0

  for (it in 1:maxit) {
    lin_pred <- x %*% beta
    gradient <- t(x) %*% (lin_pred - y)

    g <- 0.5*norm(y - lin_pred, "2")^2
    h <- sum(sort(abs(beta), decreasing = TRUE)*lambda)
    f <- g + h

    loss[it] <- f
    time[it] <- proc.time()[3] - start_time

    if (abs(f - opt)/f <= opt_tol)
      break

    g_old <- g

    repeat {

      beta_tilde <- sorted_l1_prox(beta - learning_rate*gradient,
                                   lambda*learning_rate)
      d <- beta_tilde - beta
      g <- 0.5*norm(y - x %*% beta_tilde, "2")^2
      q <- g_old + sum(d*gradient) + sum(d^2)/(2*learning_rate)

      if (q >= g*(1 - 1e-12))
        break

      learning_rate <- learning_rate*eta;
    }

    t <- 0.5*(1 + sqrt(1 + 4*t_old*t_old))
    beta <- beta_tilde + (t_old - 1)/t * (beta_tilde - beta_old)

    beta_old <- beta
    t_old <- t
  }

  list(loss = loss,
       time = time,
       beta = beta)

}
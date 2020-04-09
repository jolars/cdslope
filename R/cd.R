#' Coordinate Descent for SLOPE
#'
#' @inheritParams fista
#'
#' @export
cd <- function(x, y, lambda, maxit = 100, opt = -Inf, opt_tol = 1e-3) {
  p <- ncol(x)

  time <- double(maxit)
  loss <- double(maxit)

  beta <- double(p)

  start_time <- proc.time()[3]

  xTx <- colSums(x^2)

  for (it in 1:maxit) {

    lin_pred <- x %*% beta
    residual <- y - lin_pred

    g <- 0.5*norm(residual, "2")^2
    h <- sum(sort(abs(beta), decreasing = TRUE)*lambda)
    f <- g + h

    loss[it] <- f
    time[it] <- proc.time()[3] - start_time

    if (abs(f - opt)/f <= opt_tol)
      break

    for (j in 1:p) {
      residual <- residual + x[, j]*beta[j]
      beta[j] <- sum(residual*x[, j])

      beta_new <- sorted_l1_prox(beta/xTx, lambda/xTx)
      beta[j] <- beta_new[j]

      residual <- residual - x[, j]*beta[j]
    }
  }

  list(loss = loss,
       time = time,
       beta = beta)
}
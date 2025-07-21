#' Fast forecast for ARDL model using data.table
#' @export
forecast.model_ardl <- function(object, new_data, specials = NULL, ...) {
  requireNamespace("data.table")
  requireNamespace("tsibble")
  requireNamespace("distributional")

  # Extract key/index info
  key_cols <- tsibble::key_vars(new_data)
  idx_col <- tsibble::index_var(new_data)

  old_xreg <- object$xreg
  xreg <- specials$xreg[[1]]$xreg
  yvar <- object$y_name
  order <- object$order
  coef <- object$coef
  lag_starts <- object$lag_starts
  sigma2 <- object$sigma2

  # Combine old and new xreg for lagging
  xreg_full <- data.table::rbindlist(
    list(
      data.table::as.data.table(old_xreg),
      data.table::as.data.table(
        cbind(new_data[, key_cols, drop = FALSE], xreg)
      )
    ),
    use.names = TRUE, fill = TRUE
  )

  # Compute lags for all exogenous variables efficiently
  for (k_name in names(order)[names(order) != yvar]) {
    for (lag_val in lag_starts[[k_name]]:order[[k_name]]) {
      colname <- if (lag_val == 0) k_name else paste0("lag(", k_name, ", ", lag_val, ")")
      xreg_full[, (colname) := data.table::shift(.SD[[k_name]], lag_val, type = "lag"), by = key_cols]
    }
  }

  # Select only the columns needed for the forecast
  needed_cols <- c(key_cols, gsub("`", "", names(coef)[!grepl(paste0("^\\(Intercept\\)$|", yvar), names(coef))]))
  xreg_lag <- xreg_full[(.N - nrow(new_data) + 1):.N, ..needed_cols]

  # Prepare matrix for regression
  xreg_lag_matrix <- as.matrix(xreg_lag[, setdiff(names(xreg_lag), key_cols), with = FALSE])

  # Compute regression part
  alpha <- rep(unname(coef["(Intercept)"]), nrow(new_data))
  beta <- unname(coef[(order[yvar] + 1 + 1):length(coef)])
  phi <- coef[grepl(paste0("^`lag\\(", yvar), names(coef))]
  xresult <- as.numeric(alpha + xreg_lag_matrix %*% beta)

  # Prepare output table
  out <- data.table::as.data.table(new_data[, c(idx_col, key_cols), drop = FALSE])
  out[, yx := xresult]
  out[, h := seq_len(.N), by = key_cols]

  # Recursive AR filter (vectorized per group)
  ar_filter <- function(yx, phi, init) {
    stats::filter(x = yx, filter = phi, method = "recursive", init = init)
  }
  # Get initial values for AR recursion
  init_vals <- rev(unclass(object$mod_obj$model[[yvar]]))[lag_starts[[yvar]]:order[[yvar]]]
  out[, .sim := as.numeric(ar_filter(yx, coef[2:(order[yvar] - lag_starts[[yvar]] + 1 + 1)], init_vals)), by = key_cols]

  # Compute psi and variance for each forecast horizon
  p <- length(phi)
  H <- max(out$h)
  psi <- numeric(H)
  psi[1] <- 1
  for (step in 2:H) {
    s <- 0
    for (i in 1:min(p, step - 1)) {
      s <- s + phi[i] * psi[step - i]
    }
    psi[step] <- s
  }
  out[, .var := sigma2 * cumsum(psi^2)[h]]
  out[, .sd := sqrt(.var)]
  out[, .dist := ifelse(is.na(.sim), NA, distributional::dist_normal(.sim, .sd))]

  # Return only the distribution column as in the original
  out$.dist
}
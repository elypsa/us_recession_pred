# ==============================================================================
# AR(1) Probit Model — Marginal Composite Likelihood (MCL) Implementation
# ==============================================================================

#' Compute Conditional Mean μ_t Recursively
#'
#' @description
#' Computes the conditional mean of the latent variable y*_t for all time periods
#' using the recursive formula.
#'
#' @param alpha Intercept parameter
#' @param rho AR(1) coefficient (must satisfy |rho| < 1)
#' @param beta Vector of coefficients for predictors
#' @param X Matrix of predictors (T x p), already lagged by m periods
#'
#' @return Vector of conditional means μ_t for t = 1, ..., T
#'
#' @details
#' **Theory Ref:** Section 4: "Conditional Mean μ_t"
#'
#' **Equation (Recursive):**
#' $$\mu_t = \alpha + \rho \, \mu_{t-1} + x_{t-m}' \beta, \qquad \mu_0 = \frac{\alpha}{1-\rho}$$
#'
#' **Moments:**
#' $$y_t^* \mid x \sim \mathcal{N}(\mu_t, \, 1)$$
compute_mu_recursive <- function(alpha, rho, beta, X) {
  if (abs(rho) > 1) {
    warning('Abs. value of AR(1) coefficient larger than 1')
  }
  T <- nrow(X)
  mu <- numeric(T)

  # Initial condition from stationary distribution
  mu_0 <- alpha / (1 - rho)

  # Recursive computation
  for (t in 1:T) {
    if (t == 1) {
      mu[t] <- alpha + rho * mu_0 + as.numeric(X[t, ] %*% beta)
    } else {
      mu[t] <- alpha + rho * mu[t - 1] + as.numeric(X[t, ] %*% beta)
    }
  }

  return(mu)
}


#' Compute Gradient of μ_t with Respect to θ
#'
#' @description
#' Recursively computes ∂μ_t/∂θ for all time periods, where θ = (ρ, α, β')'.
#'
#' @param alpha Intercept parameter
#' @param rho AR(1) coefficient
#' @param beta Vector of coefficients for predictors
#' @param X Matrix of predictors (T x p)
#' @param mu Vector of conditional means (from compute_mu_recursive)
#'
#' @return Matrix (T x k) where k = 1 + 1 + p = length(theta), containing
#'         ∂μ_t/∂θ for each time period
#'
#' @details
#' **Theory Ref:** Section 6: "Gradient (Score)" - Derivatives of μ_t
#'
#' **Equation:**
#' $$\frac{\partial \mu_t}{\partial \theta} = \rho \, \frac{\partial \mu_{t-1}}{\partial \theta}
#' + \begin{pmatrix} \mu_{t-1} \\ 1 \\ x_{t-m} \end{pmatrix}$$
#'
#' **Initial condition:**
#' $$\frac{\partial \mu_0}{\partial \theta} = \begin{pmatrix}
#' \frac{\alpha}{(1-\rho)^2} \\ \frac{1}{1-\rho} \\ 0
#' \end{pmatrix}$$
#'
#' where blocks correspond to (∂/∂ρ, ∂/∂α, ∂/∂β).
compute_dmu_dtheta <- function(alpha, rho, beta, X, mu) {
  T <- nrow(X)
  p <- length(beta)
  k <- 1 + 1 + p # rho, alpha, beta

  dmu <- matrix(0, nrow = T, ncol = k)

  # Initial condition at t=0
  mu_0 <- alpha / (1 - rho)
  dmu_0 <- c(
    alpha / (1 - rho)^2, # ∂μ_0/∂ρ
    1 / (1 - rho), # ∂μ_0/∂α
    rep(0, p) # ∂μ_0/∂β = 0
  )

  # Recursive computation
  for (t in 1:T) {
    if (t == 1) {
      # Derivative components at t=1
      dmu_prev <- dmu_0
      mu_prev <- mu_0
    } else {
      dmu_prev <- dmu[t - 1, ]
      mu_prev <- mu[t - 1]
    }

    # Build the added term: (μ_{t-1}, 1, x_{t-m})'
    added_term <- c(mu_prev, 1, X[t, ])

    # Recursive formula
    dmu[t, ] <- rho * dmu_prev + added_term
  }

  return(dmu)
}


#' Compute Second Derivative of μ_t with Respect to θ
#'
#' @description
#' Recursively computes ∂²μ_t/∂θ∂θ' for all time periods.
#'
#' @param alpha Intercept parameter
#' @param rho AR(1) coefficient
#' @param beta Vector of coefficients for predictors
#' @param X Matrix of predictors (T x p)
#' @param dmu Matrix of first derivatives (from compute_dmu_dtheta)
#'
#' @return Array (T x k x k) containing the Hessian of μ_t for each time period
#'
#' @details
#' **Theory Ref:** Section 7: "Hessian" - Second derivative of μ_t
#'
#' **Equation:**
#' $$\frac{\partial^2 \mu_t}{\partial \theta \, \partial \theta'} =
#' \rho \, \frac{\partial^2 \mu_{t-1}}{\partial \theta \, \partial \theta'} +
#' e_1 \frac{\partial \mu_{t-1}}{\partial \theta'} +
#' \frac{\partial \mu_{t-1}}{\partial \theta} e_1'$$
#'
#' where $e_1 = (1, 0, \ldots, 0)'$ selects the ρ component.
compute_d2mu_dtheta2 <- function(alpha, rho, beta, X, dmu) {
  T <- nrow(X)
  p <- length(beta)
  k <- 1 + 1 + p

  d2mu <- array(0, dim = c(T, k, k))

  # e_1 vector (selects rho component)
  e1 <- c(1, rep(0, k - 1))

  # Initial condition at t=0
  mu_0 <- alpha / (1 - rho)
  dmu_0 <- c(alpha / (1 - rho)^2, 1 / (1 - rho), rep(0, p))

  # Second derivative of mu_0
  d2mu_0 <- matrix(0, k, k)
  d2mu_0[1, 1] <- 2 * alpha / (1 - rho)^3 # ∂²μ_0/∂ρ²
  d2mu_0[1, 2] <- d2mu_0[2, 1] <- 1 / (1 - rho)^2 # ∂²μ_0/∂ρ∂α

  # Recursive computation
  for (t in 1:T) {
    if (t == 1) {
      d2mu_prev <- d2mu_0
      dmu_prev <- dmu_0
    } else {
      d2mu_prev <- d2mu[t - 1, , ]
      dmu_prev <- dmu[t - 1, ]
    }

    # Recursive formula
    d2mu[t, , ] <- rho * d2mu_prev + outer(e1, dmu_prev) + outer(dmu_prev, e1)
  }

  return(d2mu)
}


#' Marginal Composite Log-Likelihood (MCL)
#'
#' @description
#' Computes the negative MCL objective (for minimization).
#'
#' @param theta Parameter vector (ρ, α, β')
#' @param y Binary outcome vector (0/1)
#' @param X Matrix of predictors (already lagged)
#'
#' @return Negative average log-likelihood (scalar)
#'
#' @details
#' **Theory Ref:** Section 5: "Marginal Composite Log-Likelihood (MCL)"
#'
#' **Equation:**
#' $$\mathcal{L}_{\mathrm{MCL}}(\theta) = \frac{1}{T} \sum_{t=1}^{T}
#' \left[ y_t \ln \Phi(\mu_t) + (1-y_t) \ln \Phi(-\mu_t) \right]$$
#'
#' This function returns the **negative** for use with minimization routines.
mcl_objective <- function(theta, y, X) {
  p <- ncol(X)
  rho <- theta[1]
  alpha <- theta[2]
  beta <- theta[3:(2 + p)]

  # Compute conditional means
  mu <- compute_mu_recursive(alpha, rho, beta, X)

  # Compute log-likelihood contributions
  ll <- y * pnorm(mu, log.p = TRUE) + (1 - y) * pnorm(-mu, log.p = TRUE)

  # Return negative average log-likelihood
  return(-mean(ll))
}


#' Gradient of MCL (Score Function)
#'
#' @description
#' Computes the gradient of the negative MCL objective.
#'
#' @param theta Parameter vector (ρ, α, β')
#' @param y Binary outcome vector (0/1)
#' @param X Matrix of predictors
#'
#' @return Gradient vector (length k)
#'
#' @details
#' **Theory Ref:** Section 6: "Gradient (Score)"
#'
#' **Equation:**
#' $$s_t(\theta) = \frac{\partial \ln f(y_t \mid x; \theta)}{\partial \theta} =
#' \left[ y_t \frac{\phi(\mu_t)}{\Phi(\mu_t)} - (1-y_t) \frac{\phi(\mu_t)}{\Phi(-\mu_t)} \right]
#' \frac{\partial \mu_t}{\partial \theta}$$
#'
#' The aggregate score is $s(\theta) = \frac{1}{T} \sum_t s_t(\theta)$.
mcl_gradient <- function(theta, y, X) {
  p <- ncol(X)
  rho <- theta[1]
  alpha <- theta[2]
  beta <- theta[3:(2 + p)]

  T <- length(y)

  # Compute conditional means and derivatives
  mu <- compute_mu_recursive(alpha, rho, beta, X)
  dmu <- compute_dmu_dtheta(alpha, rho, beta, X, mu)

  # Compute score contributions for each t
  phi_mu <- dnorm(mu)
  Phi_mu <- pnorm(mu)
  Phi_neg_mu <- pnorm(-mu)

  # Weight for each observation
  weight <- y * phi_mu / Phi_mu - (1 - y) * phi_mu / Phi_neg_mu

  # Score vector: average of weighted derivatives
  score <- colMeans(weight * dmu)

  # Return negative gradient (for minimization)
  return(-score)
}


#' Hessian of MCL
#'
#' @description
#' Computes the Hessian matrix of the negative MCL objective.
#'
#' @param theta Parameter vector (ρ, α, β')
#' @param y Binary outcome vector (0/1)
#' @param X Matrix of predictors
#'
#' @return Hessian matrix (k x k)
#'
#' @details
#' **Theory Ref:** Section 7: "Hessian"
#'
#' **Equation:**
#' $$h_t(\theta) = \frac{\partial^2 \ln f}{\partial \theta \, \partial \theta'} =
#' A_t(\theta) \, \frac{\partial^2 \mu_t}{\partial \theta \, \partial \theta'} -
#' B_t(\theta) \, \frac{\partial \mu_t}{\partial \theta} \frac{\partial \mu_t}{\partial \theta'}$$
#'
#' where:
#' $$A_t = y_t \frac{\phi(\mu_t)}{\Phi(\mu_t)} - (1-y_t) \frac{\phi(\mu_t)}{\Phi(-\mu_t)}$$
#'
#' $$B_t = y_t \left[ \frac{\mu_t \phi(\mu_t)}{\Phi(\mu_t)} + \frac{\phi(\mu_t)^2}{\Phi(\mu_t)^2} \right] +
#' (1-y_t) \left[ \frac{-\mu_t \phi(\mu_t)}{\Phi(-\mu_t)} + \frac{\phi(\mu_t)^2}{\Phi(-\mu_t)^2} \right]$$
mcl_hessian <- function(theta, y, X) {
  p <- ncol(X)
  rho <- theta[1]
  alpha <- theta[2]
  beta <- theta[3:(2 + p)]

  T <- length(y)
  k <- length(theta)

  # Compute conditional means and derivatives
  mu <- compute_mu_recursive(alpha, rho, beta, X)
  dmu <- compute_dmu_dtheta(alpha, rho, beta, X, mu)
  d2mu <- compute_d2mu_dtheta2(alpha, rho, beta, X, dmu)

  # Compute quantities needed for A_t and B_t
  phi_mu <- dnorm(mu)
  Phi_mu <- pnorm(mu)
  Phi_neg_mu <- pnorm(-mu)

  # A_t weights
  A_t <- y * phi_mu / Phi_mu - (1 - y) * phi_mu / Phi_neg_mu

  # B_t weights
  B_t <- y *
    (mu * phi_mu / Phi_mu + phi_mu^2 / Phi_mu^2) +
    (1 - y) * (-mu * phi_mu / Phi_neg_mu + phi_mu^2 / Phi_neg_mu^2)

  # Accumulate Hessian
  H <- matrix(0, k, k)
  for (t in 1:T) {
    H <- H + A_t[t] * d2mu[t, , ] - B_t[t] * outer(dmu[t, ], dmu[t, ])
  }

  # Average and return negative (for minimization)
  return(-H / T)
}


#' Compute Sandwich Variance Estimator (HAC)
#'
#' @description
#' Computes robust asymptotic variance using Newey-West HAC estimator with
#' Parzen kernel.
#'
#' @param theta Parameter estimates
#' @param y Binary outcome vector
#' @param X Matrix of predictors
#'
#' @return List containing:
#'   - vcov: Variance-covariance matrix
#'   - se: Standard errors
#'
#' @details
#' **Theory Ref:** Section 8: "Asymptotic Variance (Sandwich Estimator)"
#'
#' **Equation:**
#' $$\widehat{\mathrm{Var}}(\hat{\theta}) = \hat{H}^{-1} \, \hat{\Omega} \, \hat{H}^{-1}$$
#'
#' - $\hat{H} = \frac{1}{T} \sum_t h_t(\hat{\theta})$ — estimated Hessian
#' - $\hat{\Omega}$ — HAC long-run variance of scores
#' - Bandwidth: $M_T = \lfloor 12(T/100)^{4/25} \rfloor$, Parzen kernel
compute_sandwich_variance <- function(theta, y, X) {
  p <- ncol(X)
  rho <- theta[1]
  alpha <- theta[2]
  beta <- theta[3:(2 + p)]

  T <- length(y)
  k <- length(theta)

  # Compute conditional means and derivatives
  mu <- compute_mu_recursive(alpha, rho, beta, X)
  dmu <- compute_dmu_dtheta(alpha, rho, beta, X, mu)

  # Compute individual score contributions
  phi_mu <- dnorm(mu)
  Phi_mu <- pnorm(mu)
  Phi_neg_mu <- pnorm(-mu)

  weight <- y * phi_mu / Phi_mu - (1 - y) * phi_mu / Phi_neg_mu
  scores <- weight * dmu # T x k matrix

  # Compute Hessian
  H <- -mcl_hessian(theta, y, X) # Remove negative sign

  # HAC variance of scores (Newey-West with Parzen kernel)
  M_T <- floor(12 * (T / 100)^(4 / 25))

  # Parzen kernel weights
  parzen_weight <- function(j, M) {
    x <- j / M
    if (abs(x) <= 0.5) {
      return(1 - 6 * x^2 + 6 * abs(x)^3)
    } else if (abs(x) <= 1) {
      return(2 * (1 - abs(x))^3)
    } else {
      return(0)
    }
  }

  # Compute Omega (long-run variance)
  Omega <- crossprod(scores) / T # Variance at lag 0

  for (j in 1:M_T) {
    w <- parzen_weight(j, M_T)
    Gamma_j <- matrix(0, k, k)

    for (t in (j + 1):T) {
      Gamma_j <- Gamma_j + outer(scores[t, ], scores[t - j, ])
    }
    Gamma_j <- Gamma_j / T

    Omega <- Omega + w * (Gamma_j + t(Gamma_j))
  }

  # Sandwich variance
  H_inv <- solve(H)
  vcov <- H_inv %*% Omega %*% H_inv / T
  se <- sqrt(diag(vcov))

  return(list(vcov = vcov, se = se))
}


#' Estimate AR(1) Probit Model via MCL
#'
#' @description
#' Main estimation function for the AR(1) Probit model using Marginal Composite
#' Likelihood. Returns an S3 object of class 'ar_probit_mcl'.
#'
#' @param formula Formula object specifying the model (e.g., y ~ x1 + x2)
#' @param data Data frame containing the variables
#' @param lag_m Number of periods to lag the predictors (default: 1)
#' @param start Optional starting values for (rho, alpha, beta). If NULL, uses
#'              default values: rho=0, alpha=qnorm(mean(y)), beta=0
#' @param ... Additional arguments passed to optim()
#'
#' @return Object of class 'ar_probit_mcl' containing:
#'   - coefficients: Named vector of parameter estimates
#'   - vcov: Variance-covariance matrix
#'   - se: Standard errors
#'   - convergence: Convergence code from optim
#'   - loglik: Final MCL log-likelihood value
#'   - formula: Model formula
#'   - data: Data used for estimation
#'   - lag_m: Lag parameter
#'   - mu: Fitted conditional means
#'   - call: The matched call
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- estimate_ar_probit_mcl(
#'   recession_flag ~ yield_spread + NFCI,
#'   data = df_quarterly,
#'   lag_m = 4
#' )
#' summary(fit)
#' }
estimate_ar_probit_mcl <- function(
  formula,
  data,
  lag_m = 1,
  start = NULL,
  ...
) {
  # Store the call
  cl <- match.call()

  # Extract model frame
  mf <- model.frame(formula, data, na.action = na.omit)
  y <- model.response(mf)
  X_orig <- model.matrix(formula, mf)[, -1, drop = FALSE] # Remove intercept

  # Create lagged predictors
  T <- length(y)
  p <- ncol(X_orig)

  # Apply lag
  if (lag_m > 0) {
    # Shift predictors back by lag_m periods
    X <- X_orig[1:(T - lag_m), , drop = FALSE]
    y <- y[(lag_m + 1):T]
    T_eff <- length(y)
  } else {
    X <- X_orig
    T_eff <- T
  }

  # Set starting values
  if (is.null(start)) {
    mean_y <- mean(y)
    mean_y <- pmax(pmin(mean_y, 0.99), 0.01) # Bound away from 0/1

    start <- c(
      rho = 0,
      alpha = qnorm(mean_y),
      beta = rep(0, p)
    )
  }

  # Parameter names
  var_names <- colnames(X)
  param_names <- c("rho", "alpha", var_names)
  names(start) <- param_names

  # Optimization with constraint |rho| < 1
  opt_result <- optim(
    par = start,
    fn = mcl_objective,
    gr = mcl_gradient,
    method = "L-BFGS-B",
    lower = c(-0.99, -Inf, rep(-Inf, p)),
    upper = c(0.99, Inf, rep(Inf, p)),
    y = y,
    X = X,
    hessian = FALSE, # We compute analytical Hessian
    ...
  )

  # Extract estimates
  theta_hat <- opt_result$par
  names(theta_hat) <- param_names

  # Compute sandwich variance
  var_est <- compute_sandwich_variance(theta_hat, y, X)
  names(var_est$se) <- param_names
  rownames(var_est$vcov) <- colnames(var_est$vcov) <- param_names

  # Compute fitted mu values
  rho_hat <- theta_hat[1]
  alpha_hat <- theta_hat[2]
  beta_hat <- theta_hat[3:(2 + p)]
  mu_fitted <- compute_mu_recursive(alpha_hat, rho_hat, beta_hat, X)

  # Create return object
  result <- list(
    coefficients = theta_hat,
    vcov = var_est$vcov,
    se = var_est$se,
    convergence = opt_result$convergence,
    loglik = -opt_result$value * T_eff, # Convert back to log-likelihood
    n_obs = T_eff,
    formula = formula,
    data = data,
    lag_m = lag_m,
    mu = mu_fitted,
    y = y,
    X = X,
    var_names = var_names,
    call = cl
  )

  class(result) <- "ar_probit_mcl"
  return(result)
}


#' Fitted Values for AR(1) Probit MCL Model
#'
#' @description
#' Extracts in-sample fitted values (predicted probabilities or classes).
#'
#' @param object Object of class 'ar_probit_mcl'
#' @param type Either "prob" for probabilities or "class" for binary predictions
#' @param threshold Threshold for classification (default: 0.5)
#' @param ... Additional arguments (not used)
#'
#' @return Vector of fitted values
#'
#' @export
fitted.ar_probit_mcl <- function(
  object,
  type = c("prob", "class"),
  threshold = 0.5,
  ...
) {
  type <- match.arg(type)

  # Compute probabilities
  prob <- pnorm(object$mu)

  if (type == "prob") {
    return(prob)
  } else {
    return(as.numeric(prob >= threshold))
  }
}


#' Predictions for AR(1) Probit MCL Model (Time Series Forecasting)
#'
#' @description
#' Generates h-step-ahead forecasts from time T (end of training period).
#' This is a TRUE time series forecast: uses ONLY information available at
#' time T, does NOT use future data.
#'
#' @param object Object of class 'ar_probit_mcl' fitted on data up to time T
#' @param newdata Not used. Kept for S3 method compatibility.
#' @param h Forecast horizon(s). Can be a single integer or vector of horizons
#'   (e.g., 1:4 for 1-4 quarters ahead). Default is 1.
#' @param type Either "prob" for probabilities or "class" for binary predictions
#' @param threshold Threshold for classification (default: 0.5)
#' @param ... Additional arguments (not used)
#'
#' @return Named vector of forecasts for horizons h. Names are "h=1", "h=2", etc.
#'   If h is a single value, returns a single unnamed value.
#'
#' @details
#' **Theory Ref:** Section 9: "Forecasting"
#'
#' **TIME SERIES FORECASTING LOGIC:**
#'
#' The model is fitted on data up to time T. At time T, we know:
#' - μ_T (last conditional mean from training)
#' - x_T, x_{T-1}, ..., x_{T-m+1} (recent predictors from training data)
#'
#' To forecast h steps ahead, we compute:
#'
#' **Equation:**
#' $$P_T(y_{T+h}=1) = \Phi\left( \frac{\mu_{T+h|T}}{\sqrt{1-\rho^{2h}}} \right)$$
#'
#' where:
#' $$\mu_{T+h|T} = \frac{1-\rho^h}{1-\rho} \alpha + \rho^h \mu_T +
#' \sum_{k=0}^{h-1} \rho^k x_{T+h-m-k}' \beta$$
#'
#' **Key constraint:** We need x_{T+h-m-k} for k=0,...,h-1. These are:
#' x_{T+h-m}, x_{T+h-m-1}, ..., x_{T+1-m}.
#'
#' For these to be known at time T, we need T+h-m ≤ T, which gives h ≤ m.
#'
#' **Example with m=4 (4-quarter lag):**
#' - h=1: need x_{T-3} (known at T) ✓
#' - h=2: need x_{T-2}, x_{T-3} (known at T) ✓
#' - h=4: need x_T, x_{T-1}, x_{T-2}, x_{T-3} (all known at T) ✓
#'
#' **NO FUTURE DATA IS USED.** All forecasts use only information up to time T.
#'
#' @export
predict.ar_probit_mcl <- function(
  object,
  newdata = NULL,
  h = 1,
  type = c("prob", "class"),
  threshold = 0.5,
  ...
) {
  type <- match.arg(type)

  # Extract parameters
  rho <- object$coefficients["rho"]
  alpha <- object$coefficients["alpha"]
  beta <- object$coefficients[object$var_names]
  lag_m <- object$lag_m

  # Get the last fitted mu from training data (this is μ_T)
  mu_T <- object$mu[length(object$mu)]

  # Validate horizons don't exceed lag
  if (max(h) > lag_m) {
    warning(
      "Maximum forecast horizon (", max(h), ") exceeds lag_m (", lag_m, "). ",
      "Forecasts beyond h=", lag_m, " require unobserved future predictors."
    )
  }

  # Get the last lag_m observations of predictors from TRAINING data
  # These are x_T, x_{T-1}, ..., x_{T-m+1}
  T_train <- nrow(object$X)
  X_recent <- object$X[max(1, T_train - lag_m + 1):T_train, , drop = FALSE]

  # Compute forecasts for each horizon
  forecasts <- numeric(length(h))
  names(forecasts) <- paste0("h=", h)

  for (i in seq_along(h)) {
    h_curr <- h[i]

    # Compute forecast mean μ_{T+h|T}
    # Term 1: (1-ρ^h)/(1-ρ) * α
    term1 <- (1 - rho^h_curr) / (1 - rho) * alpha

    # Term 2: ρ^h * μ_T
    term2 <- rho^h_curr * mu_T

    # Term 3: Σ_{k=0}^{h-1} ρ^k x_{T+h-m-k}' β
    # We need x_{T+h-m}, x_{T+h-m-1}, ..., x_{T+1-m}
    term3 <- 0

    for (k in 0:(h_curr - 1)) {
      # Time index we need: T + h - m - k
      # Relative to T: h - m - k (can be negative, zero, or positive)
      # In X_recent (which starts at T-m+1 and ends at T):
      #   T is at index nrow(X_recent)
      #   T-1 is at index nrow(X_recent)-1, etc.

      time_offset <- h_curr - lag_m - k # Relative to T

      # Convert to index in X_recent (1-indexed, where last row is T)
      idx_in_recent <- nrow(X_recent) + time_offset

      if (idx_in_recent >= 1 && idx_in_recent <= nrow(X_recent)) {
        term3 <- term3 + rho^k * sum(X_recent[idx_in_recent, ] * beta)
      } else if (time_offset > 0) {
        # Would need future data - this shouldn't happen if h ≤ m
        warning(
          "Horizon h=", h_curr, " requires predictor at T+", time_offset,
          " which is beyond the training period. Setting contribution to 0."
        )
      }
    }

    # Compute conditional mean
    mu_forecast <- term1 + term2 + term3

    # Compute forecast probability with proper scaling
    scaling <- sqrt(1 - rho^(2 * h_curr))
    prob_forecast <- pnorm(mu_forecast / scaling)

    if (type == "prob") {
      forecasts[i] <- prob_forecast
    } else {
      forecasts[i] <- as.numeric(prob_forecast >= threshold)
    }
  }

  # Return forecasts
  # If single horizon, return unnamed scalar for simplicity
  if (length(h) == 1) {
    return(unname(forecasts))
  } else {
    return(forecasts)
  }
}


#' Print Method for AR(1) Probit MCL
#'
#' @param x Object of class 'ar_probit_mcl'
#' @param ... Additional arguments (not used)
#'
#' @export
print.ar_probit_mcl <- function(x, ...) {
  cat("\nAR(1) Probit Model (MCL Estimation)\n")
  cat("=====================================\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Coefficients:\n")
  print(round(x$coefficients, 4))
  cat("\n")
  cat("Log-likelihood:", round(x$loglik, 2), "\n")
  cat("Observations:", x$n_obs, "\n")
  cat("Convergence code:", x$convergence, "\n")
  if (x$convergence != 0) {
    cat("Warning: Optimization may not have converged!\n")
  }
  invisible(x)
}


#' Summary Method for AR(1) Probit MCL
#'
#' @param object Object of class 'ar_probit_mcl'
#' @param ... Additional arguments (not used)
#'
#' @export
summary.ar_probit_mcl <- function(object, ...) {
  # Compute z-statistics and p-values
  z_stat <- object$coefficients / object$se
  p_val <- 2 * pnorm(-abs(z_stat))

  # Create coefficient table
  coef_table <- cbind(
    Estimate = object$coefficients,
    `Std. Error` = object$se,
    `z value` = z_stat,
    `Pr(>|z|)` = p_val
  )

  # Create summary object
  summary_obj <- list(
    call = object$call,
    coefficients = coef_table,
    loglik = object$loglik,
    n_obs = object$n_obs,
    convergence = object$convergence,
    lag_m = object$lag_m
  )

  class(summary_obj) <- "summary.ar_probit_mcl"
  return(summary_obj)
}


#' Print Summary Method for AR(1) Probit MCL
#'
#' @param x Object of class 'summary.ar_probit_mcl'
#' @param ... Additional arguments (not used)
#'
#' @export
print.summary.ar_probit_mcl <- function(x, ...) {
  cat("\nAR(1) Probit Model (MCL Estimation)\n")
  cat("=====================================\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Coefficients:\n")
  printCoefmat(x$coefficients, digits = 4, signif.stars = TRUE)
  cat("\n")
  cat("---\n")
  cat("Log-likelihood:", round(x$loglik, 2), "\n")
  cat("Number of observations:", x$n_obs, "\n")
  cat("Predictor lag (m):", x$lag_m, "period(s)\n")
  cat("Convergence code:", x$convergence, "\n")
  if (x$convergence != 0) {
    cat("\nWarning: Optimization may not have converged!\n")
  }
  invisible(x)
}

#' Penalized synthetic control estimator
#'
#' For a given set of variable weights (v) this function estimates
#' the unit weights for a synthetic control with penalization
#' according to Abadie & L'Hour (2021). This function deals with only a
#' single treated unit.
#'
#' @param X1 `N_covars by 1 matrix` of treated unit covariates
#' @param X0 `N_covars by N_donors matrix` of donor unit covariates
#' @param v `N_covars vector` of variable weights
#' @param lambda `numeric` penalization parameter
#' @param opt_pars `clarabel` settings using [clarabel::clarabel_control()]
#' @param standardize `boolean` whether to standardize the input matrices (default TRUE)
#'
#' @details This routine uses the same notation of the original [Synth::synth()] implementation
#' but uses a different, faster quadratic program solver (namely, [osqp::osqp()]). Additionally, it
#' implements the penalization procedure described in Abadie & L'Hour (2021), such that the loss
#' function is as in equation 5 of that paper (but for a single treated unit).
#'
#' Variable weights are not optimized by this function, meaning they need to be pre-specified.
#' This is by design.
#'
#' The original synthetic control method can be recovered by setting lambda = 0. For determining
#' lambda based on data, see [cv_pensynth()].
#'
#' @references Abadie, A., & L’Hour, J. (2021).
#' A penalized synthetic control estimator for disaggregated data.
#' _Journal of the American Statistical Association, 116_(536), 1817-1834.
#'
#' @return A list with two values: `w`, the estimated weights; and
#' `solution`, the result of the optimization.
#'
#' @importFrom utils capture.output
#'
#' @examples
#' # generate some data
#' X0 <- matrix(
#'   c(1, 1.3,
#'     0.5, 1.8,
#'     1.1, 2.4,
#'     1.8, 1.8,
#'     1.3, 1.8), 2)
#' X1 <- matrix(c(0.8, 1.65), 2)
#' v <- rep(1, 2)
#'
#' # run classic synthetic control (no penalization)
#' res <- pensynth(X1, X0, v)
#' plot(t(X0))
#' points(t(X1), pch = 2)
#' points(t(X0%*%res$w), pch = 3)
#'
#' # run synthetic control with penalty
#' res <- pensynth(X1, X0, v, lambda = 0.5)
#' points(t(X0 %*% res$w), pch = 4)
#'
#' @seealso [cv_pensynth()] [Synth::synth()]
#'
#' @export
pensynth <- function(X1, X0, v, lambda = 0, return_solver_info= TRUE, opt_pars = clarabel::clarabel_control(), standardize = TRUE) {
  if (standardize) {
    st <- standardize_X(X1, X0)
    X0 <- st$X0
    X1 <- st$X1
  }
  N_donors <- ncol(X0)
  X0v <- X0*sqrt(v)
  X1v <- X1*sqrt(v)

  # components for quadratic program
  # see https://github.com/jeremylhour/pensynth/blob/master/functions/wsoll1.R
  X0VX0 <- crossprod(X0v)
  X1VX0 <- crossprod(X1v, X0v)
  Delta <- apply(X0v - c(X1v), 2, crossprod)

  # Constraint matrices
  Amat <- rbind(
    rep(1, N_donors), # Sum to 1 constraint
    -diag(N_donors) # Individ. weights gte 0 constraint
  )
  B <- c(
    1, # Sum to 1 constraint
    rep(0, N_donors) # Individ. weights gte 0 constraint
  )

  # Define function for solving qp for a given lambda
  solve_qp <- function(lambda) {
    delta <- 1e-6
    P_reg <- X0VX0 + Matrix::Diagonal(n = N_donors, x = delta)

    # run the quadratic program solver
    result <- clarabel::clarabel(
      P = P_reg,
      q = -X1VX0 + lambda*Delta,
      A = Amat,
      b = B,
      cones = list(
        z = 1L, # There is 1 equality
        l = N_donors # There are N_donors inequalities
      ),
      control = opt_pars
    )

    # clarabel only returns a numeric status code, so we'll add a
    # human-readable status column here (plus a description)
    result$status_description <- clarabel::solver_status_descriptions()[result$status][[1]]
    result$status <- names(clarabel::solver_status_descriptions()[result$status])

    # Return result
    return(result)
  }

  solver_output <- solve_qp(lambda)


  # Construct a list of outputs
  out_obj <- list(
      w_opt    = as.matrix(solver_output[["x"]]),
      l_opt    = lambda,
      lseq     = lambda,
      w_path   = as.matrix(solver_output[["x"]])
  )

  # If we've been requested to return info about the solving process, do so
  if (return_solver_info) {
    # Remove unneeded columns from the solver output matrix
    rows_to_drop <- c("x", "y", "s", "z")
    solver_output <- as.matrix(solver_output[!names(solver_output) %in% rows_to_drop])

    # Add each row from the solver output matrix to .Data
    for (i in seq_len(nrow(solver_output))) {
      row_name <- rownames(solver_output)[i]
      out_obj[[row_name]] <- unlist(solver_output[i, ])
    }
  }

  # Convert the list to a cvpensynth object
  out_obj <- structure(
    .Data = out_obj,
    class = "pensynth"
  )

  return(out_obj)
}

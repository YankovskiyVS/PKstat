#' Parameters of the pharmacokinetics
#'
#' Use the model of 2 compartment
#'
#' @param data  data frame, which include time and concentrations points
#' @param time  time row from data frame
#' @param conc  concentration row from data frame
#' @examples
#' \dontrun{
#' long_running_function()
#' }
#' @return PK parameters%: "A", "alpha", "B", "beta", "R^2", "AUC", "K21", "K10",
#' "K12", "V1", "V2", "Vdss", "CL1", "CL2", "T1/2 alpha", "T1/2 beta"
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs
#' @importFrom stats lm coef predict
#' @importFrom utils head tail
#' @importFrom minpack.lm nlsLM
#' @export
par_2c <- function(data, time, conc, dose) {
  data$conc <- data[[conc]]
  fit <- minpack.lm::nlsLM(conc ~ A*exp(-alpha*data[[time]]) + B*exp(-beta*data[[time]]),
                           data = data,
                           start = list(A = 9, B = 3, alpha = 0.6, beta = 0.2),
                           control = list(maxiter = 200))
  A <- coef(fit)[1]
  B <- coef(fit)[2]
  alpha <- coef(fit)[3]
  beta <- coef(fit)[4]

  # Calculate predicted concentrations
  data$predicted_conc <- A * exp(-alpha * data[[time]]) + B * exp(-beta*data[[time]])

  #R squared calculation
  data$residuals <- data[[conc]] - data$predicted_conc
  SS_residuals <- sum(data$residuals^2)
  SS_total <- sum((data[[conc]] - mean(data[[conc]]))^2)
  R_squared <- 1 - (SS_residuals / SS_total)

  # Calculate AUC using the trapezoidal rule
  auc <- A/alpha + B/beta

  AUMC <- (A/alpha^2) + (B/beta^2)
  MRT <- AUMC/auc

  # Calculate pharmacokinetic parameters
  K21 <- (A*alpha+B*beta)/(A+B)
  K10 <- (alpha*beta)/K21
  K12 <- alpha + beta - K21 - K10
  V1 <- 1/(A+B)
  V2 <- (V1/K21)*K12
  CL1 <- V1*K10
  CL2 <- V2*K12
  T_halph_alpha <- log(2)/alpha
  T_halph_beta <- log(2)/beta
  T_halph <- MRT * log(2)
  Vc <- dose/(A+B)
  CL <- dose/auc
  Vdss <- Vc*(1 + K12/K21)
  

  # Compile results into a data frame
  parameters <- data.frame(
    parameter = c("A", "alpha", "B", "beta", "R^2", "AUC", "K21", "K10",
                  "K12", "V1", "V2", "CL1", "CL2", "T1/2 alpha", "T1/2 beta", "Vc", "CL", "Vdss",
                 "T1/2"),
    value = c(
      A,
      alpha,
      B,
      beta,
      R_squared,
      auc,
      K21,
      K10,
      K12,
      V1,
      V2,
      CL1,
      CL2,
      T_halph_alpha,
      T_halph_beta,
      Vc,
      CL,
      Vdss,
      T_halph
    )
  )

  return(parameters)
}

#' Regression model of the pharmacokinetics
#' Use the model of 2 compartment
#' @param data data frame, which include time and concentrations points
#' @param time time row from data frame
#' @param conc concentration row from data frame
#' @return an exponential regression model of experimental data
#' @examples
#' \dontrun{
#' long_running_function()
#' }
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs
#' @importFrom stats lm coef predict
#' @importFrom utils head tail
#' @importFrom minpack.lm nlsLM
#' @export
graph_2c <- function(data, time, conc){
  data$conc <- data[[conc]]
  fit <- minpack.lm::nlsLM(conc ~ A*exp(-alpha*data[[time]]) + B*exp(-beta*data[[time]]),
                           data = data,
                           start = list(A = 9, B = 3, alpha = 0.6, beta = 0.2),
                           control = list(maxiter = 200))

  # Extract coefficients from the fitted model
  coef_fit <- coef(fit)
  A <- coef_fit["A"]
  B <- coef_fit["B"]
  alpha <- coef_fit["alpha"]
  beta <- coef_fit["beta"]

  # Calculate predicted concentrations
  data$predicted_conc <- A * exp(-alpha * data[[time]]) + B * exp(-beta * data[[time]])

  graph <- ggplot2::ggplot(data = data, aes(x = .data[[time]], y = .data[[conc]])) +
    ggplot2::geom_point(color = "black", size = 2) +
    ggplot2::geom_line(aes(y = .data$predicted_conc), color = "red") +
    ggplot2::labs(title = "Exponential Regression",
                  x = colnames(data)[1],
                  y = colnames(data)[2])

  return(graph)
}

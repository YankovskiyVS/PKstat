#' Parameters of the pharmacokinetics
#'
#' Use the model of 1 compartment
#'
#' @param data  data frame, which include time and concentrations points
#' @param time  time row from data frame
#' @param conc  concentration row from data frame
#' @examples
#' \dontrun{
#' long_running_function()
#' }
#' @return PK parameters%: A, alpha, R^2 of the regression line, AUC, CL1, V1, K10, T1/2, Cmax
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs
#' @importFrom stats lm coef predict
#' @importFrom utils head tail
#' @export
par_1c <- function(data, time, conc) {
  # Calculate the log-transformed concentration
  data$log_conc <- log(data[[conc]])

  # Fit the linear model to the log-transformed concentration data
  model <- lm(log_conc ~ data[[time]], data = data)

  # Extract coefficients
  intercept <- coef(model)[1]
  slope <- coef(model)[2]

  # Transform the intercept to the original scale
  A <- exp(intercept)
  alpha <- -slope

  # Calculate predicted concentrations
  data$predicted_conc <- A * exp(-alpha * data[[time]])

  # Calculate AUC using the trapezoidal rule
  auc <- sum(diff(data[[time]]) * (head(data$predicted_conc, -1) + tail(data$predicted_conc, -1)) / 2)

  # Calculate pharmacokinetic parameters
  CL1 <- alpha / A             # Clearance
  V1 <- 1 / A                  # Volume of distribution
  K10 <- alpha                 # Elimination rate constant
  T_half <- log(2) / alpha     # Half-life
  Cmax <- max(data$predicted_conc)  # Maximum concentration

  # Compile results into a data frame
  parameters <- data.frame(
    parameter = c("A", "alpha", "R^2", "AUC", "CL1", "V1", "K10", "T1/2", "Cmax"),
    value = c(
      A,
      alpha,
      summary(model)$adj.r.squared,
      auc,
      CL1,
      V1,
      K10,
      T_half,
      Cmax
    )
  )

  return(parameters)
}
#' Regression model of the pharmacokinetics
#' Use the model of 1 compartment
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
#' @export
graph_1c <- function(data, time, conc){
  # Calculate the log-transformed concentration
  data$log_conc <- log(data[[conc]])

  # Fit the linear model to the log-transformed concentration data
  model <- lm(log_conc ~ data[[time]], data = data)

  # Extract coefficients
  intercept <- coef(model)[1]
  slope <- coef(model)[2]

  # Transform the intercept to the original scale
  A <- exp(intercept)
  alpha <- -slope

  # Calculate predicted concentrations
  data$predicted_conc <- A * exp(-alpha * data[[time]])

  graph <- ggplot2::ggplot(data = data, aes(x = .data[[time]], y = .data[[conc]])) +
    ggplot2::geom_point(color = "black", size = 2) +
    ggplot2::geom_line(aes(y = .data$predicted_conc), color = "red") +
    ggplot2::labs(title = "Exponential Regression",
                  x = colnames(data)[1],
                  y = colnames(data)[2])

  return(graph)
}

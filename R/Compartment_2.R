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
#' @export
par_2c <- function(data, time, conc) {
  fit <- nls(conc ~ A * exp(-alpha * data[[time]]) + B * exp(-beta * data[[time]]), 
             data = data,
             start = list(A = 9, B = 3, alpha = 0.6, beta = 0.2))
  summary(fit)
  A <- coef(fit)[1]
  B <- coef(fit)[2]
  alpha <- coef(fit)[3]
  beta <- coef(fit)[4]
  
  # Calculate predicted concentrations
  data$predicted_conc <- A * exp(-alpha * data[[time]]) + B * exp(-beta*data[[time]])
  
  #R squared calculation
  df$residuals <- df$conc - df$predicted_conc
  SS_residuals <- sum(df$residuals^2)
  SS_total <- sum((df$conc - mean(df$conc))^2)
  R_squared <- 1 - (SS_residuals / SS_total)
  
  # Calculate AUC using the trapezoidal rule
  auc <- sum(diff(data[[time]]) * (head(data$predicted_conc, -1) + tail(data$predicted_conc, -1)) / 2)
  
  # Calculate pharmacokinetic parameters
  K21 <- (A*alpha+B*beta)/(A+B)
  K10 <- (alpha*beta)/K21
  K12 <- alpha + beta - K21 - K10
  V1 <- 1/(A+B)
  V2 <- (V1/K21)*K12
  Vdss <- V1 + V2
  CL1 <- V1*K10
  CL2 <- V2*K12
  T_halph_alpha <- log(2)/alpha
  T_halph_beta <- log(2)/beta
  
  # Compile results into a data frame
  parameters <- data.frame(
    parameter = c("A", "alpha", "B", "beta", "R^2", "AUC", "K21", "K10", 
                  "K12", "V1", "V2", "Vdss", "CL1", "CL2", "T1/2 alpha", "T1/2 beta"),
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
      Vdss,
      CL1,
      CL2,
      T_halph_alpha,
      T_halph_beta
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
#' @export
graph_2c <- function(data, time, conc){
  fit <- nls(conc ~ A * exp(-alpha * data[[time]]) + B * exp(-beta * data[[time]]), 
             data = data,
             start = list(A = 9, B = 3, alpha = 0.6, beta = 0.2))
  
  # Calculate predicted concentrations
  data$predicted_conc <- A * exp(-alpha * data[[time]]) + B * exp(-beta*data[[time]])
  
  graph <- ggplot2::ggplot(data = data, aes(x = .data[[time]], y = .data[[conc]])) +
    ggplot2::geom_point(color = "black", size = 2) +
    ggplot2::geom_line(aes(y = .data$predicted_conc), color = "red") +
    ggplot2::labs(title = "Exponential Regression",
                  x = colnames(data)[1],
                  y = colnames(data)[2])
  
  return(graph)
}

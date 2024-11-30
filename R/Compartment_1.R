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
par_compartment_1 <- function(data, time, conc){
  df <- data
  df$log_conc <- log(df$conc)
  model <- lm(log_conc ~ time, data = df)
  summary(model)
  intercept <- coef(model)[1]
  slope <- coef(model)[2]

  A <- exp(intercept)
  alpha <- -slope
  df$predicted_conc <- A * exp(-alpha * df$time)

  auc <- sum(diff(df$time) * (head(df$predicted_conc, -1) + tail(df$predicted_conc, -1)) / 2)
  CL1 <-alpha/A
  V1 <- 1/A
  K10 <- alpha
  t <- log(2)/alpha
  Cmax <- predict(model)
  parameters = data.frame(parameter = c("A", "alpha", "R^2", "AUC", "CL1", "V1", "K10", "T1/2", "Cmax"),
                          value = c(A, alpha, summary(model)$adj.r.squared, auc, CL1, V1, K10, t, exp(Cmax[1])))
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
graph_compartment_1 <- function(data, time, conc){
  df <- data
  df$log_conc <- log(df$conc)
  model <- lm(log_conc ~ time, data = df)
  summary(model)
  intercept <- coef(model)[1]
  slope <- coef(model)[2]

  A <- exp(intercept)
  alpha <- -slope
  df$predicted_conc <- A * exp(-alpha * df$time)

  graph <- ggplot2::ggplot(df, aes(x = time, y = conc)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(aes(y = predicted_conc), color = "red") +
    ggplot2::labs(title = "Exponential Regression", x = "Time, min", y = "Concentration, ng/ml")
  return(graph)
}

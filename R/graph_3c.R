#' Regression model of the pharmacokinetics
#'
#' Use the model of 3 compartment
#'
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
graph_3c <- function(data, time, conc) {
  data$conc <- data[[conc]]
  fit <- minpack.lm::nlsLM(conc ~ A*exp(-alpha*data[[time]]) + B*exp(-beta*data[[time]]) + C*exp(-gamma*data[[time]]),
                           data = data,
                           start = list(A=10, B=5, C=1, alpha=0.5, beta=0.1, gamma=0.01),
                           control = list(maxiter = 200))
  A <- coef(fit)[1]
  B <- coef(fit)[2]
  C <- coef(fit)[3]
  alpha <- coef(fit)[4]
  beta <- coef(fit)[5]
  gamma <- coef(fit)[6]

  # Calculate predicted concentrations
  data$predicted_conc <- A * exp(-alpha * data[[time]]) + B * exp(-beta*data[[time]]) +
    C*exp(-gamma*data[[time]])
  #R squared calculation
  data$residuals <- data[[conc]] - data$predicted_conc
  SS_residuals <- sum(data$residuals^2)
  SS_total <- sum((data[[conc]] - mean(data[[conc]]))^2)
  R_squared <- 1 - (SS_residuals / SS_total)
  graph <- ggplot2::ggplot(data = data, aes(x = .data[[time]], y = .data[[conc]])) +
    ggplot2::geom_point(color = "black", size = 2) +
    ggplot2::geom_line(aes(y = .data$predicted_conc), color = "red") +
    ggplot2::labs(title = "Exponential Regression",
                  x = colnames(data)[1],
                  y = colnames(data)[2])

  return(graph)
}

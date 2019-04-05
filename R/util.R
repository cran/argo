#' logit function
#' @param x numeric value for logit transformation
#' @examples
#' logit(0.5)
#' @export
logit <- function(x) log(x) - log(1-x)

#' inverse logit function
#' @param x numeric value for inverse logit transformation
#' @examples
#' logit_inv(0)
#' @export
logit_inv <- function(x) {
  exp(x) / (1 + exp(x))
}

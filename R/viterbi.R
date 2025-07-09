#' Viterbi algorithm for Hidden Markov Models (R implementation)
#'
#' R implementation of the Viterbi algorithm for finding the most likely
#' sequence of hidden states in a Hidden Markov Model.
#'
#' @param emission Matrix of emission probabilities (log scale).
#'   Rows correspond to observations, columns to states.
#' @param transition Matrix of transition probabilities (log scale).
#'   Element [i,j] is probability of transitioning from state i to state j.
#' @param observations Vector of observation indices (1-based)
#'
#' @return Vector of most likely hidden states (0-based indexing)
#'
#' @details
#' This is the R implementation of the Viterbi algorithm. For better performance,
#' use the C++ implementation via the \code{viterbi} function (from RcppExports).
#' The algorithm uses dynamic programming to find the most probable path through
#' the hidden states given the observed data.
#'
#' @examples
#' # Simple 2-state HMM example
#' # States: 0 (low), 1 (high)
#' emission <- matrix(c(-1, -3, -3, -1), nrow = 2, ncol = 2)  # log probs
#' transition <- matrix(c(-0.1, -2.3, -2.3, -0.1), nrow = 2, ncol = 2)  # log probs
#' observations <- c(1, 1, 2, 2, 1)
#'
#' states <- viterbiR(emission, transition, observations)
#'
#' @seealso \code{\link{viterbi}} for the C++ implementation
#'
#' @export
viterbiR <- function(emission, transition, observations) {
  emission <- t(emission)
  initial <- log(rep(1 / length(emission[, 1]), length(emission[, 1])))

  numStates <- nrow(transition)
  numObs <- length(observations)

  T1 <- matrix(data = 0, nrow = numStates, ncol = numObs)
  T2 <- matrix(data = 0, nrow = numStates, ncol = numObs)
  firstObs <- observations[1]

  T1[, 1] <- initial + emission[, observations[1]]

  for (j in 2:length(observations)) {
    for (i in 1:numStates) {
      probs <- T1[, j - 1] + transition[, i] + emission[i, observations[j]]
      T1[i, j] <- max(probs)
      T2[i, j] <- which.max(probs)
    }
  }

  # MLP = most likely path
  MLP <- numeric(numObs)

  MLP[numObs] <- T2[which.max(T1[, numObs]), numObs]

  # backtrace
  for (i in numObs:2) {
    zm <- which.max(T1[, i])
    MLP[i - 1] <- T2[zm, i]
  }

  return(MLP - 1)
}

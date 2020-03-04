#' Compute a generatic transition probability function
#'
#' A closure that creates a probability transition function based on the
#' transition passed in args 1 and 2
#' @param state0 First state in the transtition
#' @param state1 Second state in the transition
#' @details The resulting function takes a trace of binary states.
#' @details If using the column asleep from a behavr object, 1 encodes asleep and 0 awake
#' @return A function that computes P(state0|state1) in a sleep trace
#' @export
generic_transition <- function(state0, state1) {
  transition_prob <- function(asleep_sequence) {
    transitions <- asleep_sequence[1:(length(asleep_sequence)-1)] - asleep_sequence[-1]
    p <- sum(transitions == state0 - state1) / sum(asleep_sequence[1:(length(asleep_sequence)-1)] == state0)
    return(p)
  }
  return(transition_prob)
}

#' Compute P(D|W), the probability that the fly changes from a wake to a doze state
#'
#' Takes as input a sequence of binary states (1/0 or T/F) assumming 1 is asleep or doze state
#' @export
p_doze <- generic_transition(0, 1)

#' Compute P(W|D), the probability that the fly changes from a doze to a wake state
#'
#' Takes as input a sequence of binary states (1/0 or T/F) assumming 1 is asleep or doze state
#' @export
p_wake <- generic_transition(1, 0)

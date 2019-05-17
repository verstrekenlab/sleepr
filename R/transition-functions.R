#' Compute a generatic transition probability fucntion
#'
#' A closure that creates a probability transition function based on the
#' transition passed in args 1 and 2
#' @param state0 First state in the transtition
#' @param state1 Second state in the transition
#' @details The resulting function takes a trace of states where 1 encodes asleep and 0 awake
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

p_doze <- generic_transition(0, 1)
p_wake <- generic_transition(1, 0)

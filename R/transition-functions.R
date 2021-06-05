#' Compute a generic transition probability function
#'
#' A function that yields another function which computes the fraction of transitions
#' starting in state0 that end in state 1
#' Based on math from https://www.pnas.org/content/117/18/10024
#' @param state0 First state in the transition
#' @param state1 Second state in the transition
#' @details The resulting function takes a trace of binary states.
#' @details If using the column asleep from a behavr object, 1 encodes asleep and 0 awake
#' @param asleep_sequence A sequence of T/F where every entry is the behavior in a bin and TRUE represents sleeping behavior
#' @return A function that computes P(state0|state1) in a sleep trace
#' @export
generic_transition <- function(state0, state1) {
  transition_prob <- function(asleep_sequence) {

    transitions <- diff(asleep_sequence)
    p <- sum(transitions == (state1 - state0)) / sum(asleep_sequence[1:(length(asleep_sequence)-1)] == state0)
    return(p)
  }
  return(transition_prob)
}

#' Compute P(D|W), the probability that the fly changes from a wake to a doze state
#'
#' Takes as input a sequence of binary states (1/0 or T/F) assumming 1 is asleep or doze state
#' @export
#' @rdname generic_transition
p_doze <- function(asleep_sequence) {}
p_doze <- generic_transition(0, 1)

#' Compute P(W|D), the probability that the fly changes from a doze to a wake state
#'
#' Takes as input a sequence of binary states (1/0 or T/F) assumming 1 is asleep or doze state
#' @export
#' @rdname generic_transition
p_wake <- function(asleep_sequence) {}
p_wake <- generic_transition(1, 0)

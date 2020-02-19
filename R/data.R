#' A Pair of Sampled flowSets
#'
#' @name example_data
#'
#' @description
#' These flowSets originated from data in
#' \url{https://flowrepository.org/id/FR-FCM-ZZGS}.
#'
#' The data were first compensated and then transformed using a custom biexponential transform
#' using \href{https://github.com/rogerswt/wadeTools}{wadeTools}.
#' They were then gated on singlet live events.  Then, each flowFrame in the
#' flowSet was sampled to 10,000 events.  Finally, only young (age <= 35 yrs)
#' was retained.
#'
#' All this was done to reduce the size of the dataset for demonstration purposes.
#'
#' @format FlowSets with 10 sampled flowFrames each.
#'
#' @source \url{https://flowrepository.org/id/FR-FCM-ZZGS}
NULL



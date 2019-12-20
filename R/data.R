#' A sampled flowSet
#' @title A sampled flowSet
#' @description
#' This flowSet originated from data in
#' \url{https://flowrepository.org/id/FR-FCM-ZZGS}.
#' The data were first compensated and then transformed using a custom biexponential transform.
#' They were then gated on singlet live events.  Then, each flowFrame in the
#' flowSet was sampled to 10,000 events.  Finally, only young (age <= 35 yrs)
#' was retained.
#'     All this was done to reduce the size of the dataset for demonstration purposes.
#'
#' @format A flowSet with 43 flowFrames.
#'
#' @source \url{https://flowrepository.org/id/FR-FCM-ZZGS}
sampled_flowset_young = "inst/extdata/sampled_flowset_young.rda"


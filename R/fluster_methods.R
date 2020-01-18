#
# fluster_methods.R
#
# declares exposed functions.
#
# 2019-12-18 WTR
#
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomics LLC 2019.                 ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without Still Pond Cytomics' express written consent.            ##
################################################################################
################################################################################

##### NOTE ######
# I think that there should be an @import igraph here, but it breaks the build
# Not sure why it's working without it...


#' @importFrom igraph normalize tree parent
#' @import igraph
#' @import flowFP
#' @import cluster
#' @import wadeTools

#' @title fluster
#' @description This function wraps most of what's needed.  It computes a fingerprint
#' model, calculates the multivariate bin centers, clusters them using
#' hierarchical agglomerative clustering (cluster::ages()), and
#' creates graph strucures for visualization.
#' @param fcs The data (either a flowFrame or a flowSet)
#' @param parameters The parameters in fcs used for analysis
#' @param nRecursions The number of recursions in calculating the fingerprint (default = 12)
#' @param nclust The number of clusters you want panoply to make.  If NULL, fluster
#' will make a guess as to the "best" number of clusters
#' @description Fluster (**F**ingerprint-based c**luster**ing)
#'  implements a workflow of doing Cytometric Fingerprint (CF) binning of
#' a flowFrame (or flowSet)
#' using flowFP to a relatively high resolution (default nRecursions of 12 results
#' in 4096 bins).  The bin centroids are then computed, using the median value of all
#' of the events in the bin for each of the included parameters.  Next, these
#' multivariate bin centroids are clustered using agglommerative hierarchical clustering
#' (cluster::agnes()).  The resulting data are represented in a graph structure
#' for visualization and interpretation.
#' @return An object of class 'fluster', with the following elements:
#' \describe{
#'   \item{mod}{The flowFPModel generated}
#'   \item{mfi}{A list representation of the bin centers}
#'   \item{centers}{The t-SNE map.  Dots represent bins}
#'   \item{clustering}{A Vector containing cluster membership of the bins}
#'   \item{graph}{A graph that can be used for visualization}
#' }
#' @examples
#' load(system.file("extdata", "sampled_flowset_young.rda", package = "fluster"))
#' flust_params = c(7:9, 11:22)
#' flust_obj = fluster(fs_young, parameters = flust_params)
#' @export
fluster = function(fcs, parameters = NULL, nRecursions = 12, nclust = NULL) {
  if (is(fcs, "flowFrame")) {
    ff = fcs
  } else if (is(fcs, "flowSet")) {
    ff = suppressWarnings(as(fcs, "flowFrame"))
    flowCore::exprs(ff) = flowCore::exprs(ff)[,which(flowCore::colnames(flowCore::exprs(ff)) != "Original")]
  } else {
    stop("Argument fcs must either be a flowFrame or a flowSet\n")
  }
  # check parameters
  if (is.null(parameters)) {
    stop("Parameters must be either a numeric or character vector\n")
    if (is.numeric(parameters)) {
      parameters = flowCore::colnames(ff)[parameters]
    }
  }

  message("computing fingerprint bins...")
  mod = flowFP::flowFPModel(ff, parameters = parameters, nRecursions = nRecursions)
  fp = flowFP::flowFP(ff, mod)
  message("calculating bin centers...")
  res = calculate_bin_phenotypes(fp = fp, fs = ff, method = "median")
  mfi = as.list(data.frame(t(res$center)))

  # agnes on bin centers
  g = build_graph(mfi = mfi)
  g = add_mfi_vertex_attributes(g, mfi)
  mst = igraph::mst(g)
  message("clustering bins...")
  mat = t(res$center)
  ag = agnes(mat)

  # check nclust
  if (is.null(nclust)) {
    nclust = advise_n_clust(ag, show = FALSE)$n1
    message("advising ", nclust, " clusters")
  }



  clusters = cutree(as.hclust(ag), k = nclust)

  # set up to plot as a graph
  message("forming a graph representation of clusters...")
  cag = agnes_to_community(ag, nclust = nclust)
  gcomm = make_graph_from_community(comm = cag, g = mst)
  gcomm = attach_layout_fr(gcomm)

  fluster = list(mod = mod, centers = mat, graph = gcomm)
  class(fluster) = "fluster"

  fluster
}

#' @title plot_fluster
#' @description Draw a picture using graph-based representation of clusters of the
#' result of fluster.
#' @param fluster The result of running fluster
#' @param vs The max size of nodes in the graph
#' @param ms The minimum size of nodes in the graph
#' @param log_size If true, scale node sizes logarithmically
#' @param vertex_frame Logical.  Should we draw a frame.
#' @param cex.main Scale factor for titles of the individual markers
#' @return N/A.
#' @examples
#' plot_fluster(fluster_obj)
#' @export
plot_fluster = function(fluster, markers = colnames(fluster$mat), vs = 10, ms = 5, log.size = FALSE, vertex.frame = TRUE, cex.main = 2, cex.lab = 2) {
  plot_comm_spread(fluster$graph, markers = colnames(fluster$mat), vs = vs, ms = ms,
                   log.size = log.size, vertex.frame = vertex.frame, cex.main = cex.main)
  draw_color_scale(cex.lab = cex.lab)
}






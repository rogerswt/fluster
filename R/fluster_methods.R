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
#' @import Rtsne
#' @import diptest
#' @import KernSmooth
#' @import igraph
#' @import fields
#' @import flowFP
#' @import cluster
#' @import wadeTools
NULL

#' @title fluster
#' @description This function wraps most of what's needed.  It computes a fingerprint
#' model, calculates the multivariate bin centers, clusters them using
#' hierarchical agglomerative clustering (cluster::ages()), and
#' creates graph strucures for visualization.
#' @param fcs The data (either a flowFrame or a flowSet)
#' @param parameters The parameters in fcs used for analysis
#' @param nRecursions The number of recursions in calculating the fingerprint (default = 12)
#' @param nclust The number of clusters you want panoply to make.  If NULL, fluster
#' @param merge Logical: should we merge initial clusters based on categorical similarity?
#' will make a guess as to the "best" number of clusters
#' @param graph Logical: should we compute the MST for visualization?
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
#'   \item{centers}{Multivariate bin centers}
#'   \item{graph}{A graph that can be used for visualization}
#'   \item{clustering}{A named list containing cluster membership of the bins}
#' }
#' @examples
#' load(system.file("extdata", "sampled_flowset_young.rda", package = "fluster"))
#' flust_params = c(7:9, 11:22)
#' flust_obj = fluster(fs_young, parameters = flust_params)
#' @export
fluster = function(fcs, parameters = NULL, nRecursions = 12, nclust = NULL, merge = TRUE, graph = TRUE, tsne = TRUE) {
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
  message("clustering bins...")
  mat = t(res$center)
  ag = agnes(mat)

  # check nclust
  if (is.null(nclust)) {
    nclust = advise_n_clust(ag, show = FALSE)$n1
    message("advising ", nclust, " clusters...")
  }

  clusters = cutree(as.hclust(ag), k = nclust)

  c_index = list()
  for (i in 1:nclust) {
    idx = which(clusters == i)
    c_index[[i]] = idx
  }
  clst = list(clst = clusters, c_index = c_index, c_centers = NULL)

  # determining modality
  modality = parameter_modality(ff, parameters = parameters)

  fluster_obj = list(mod = mod, centers = mat, agnes_obj = ag, graph = NULL, tsne = NULL, clustering = clst, modality = modality)
  class(fluster_obj) = "fluster"

  # merging clusters
  if (merge) {
    fluster_obj = merge_categorical_clusters(fluster_obj = fluster_obj, parameters = parameters)
    nclust = max(fluster_obj$clustering$clst)
    message("merging clusters, resulting in ", nclust, " clusters...")

  }

  # calculate cluster centers, mostly for visualization
  c_centers = matrix (NA, nrow = 0, ncol = length(parameters))
  colnames(c_centers) = parameters
  for (i in 1:nclust) {
    idx = which(fluster_obj$clustering$clst == i)
    c_centers = rbind(c_centers, distributions_bins(fluster_obj, parameters, bin_indices = idx)$med)
  }
  fluster_obj$clustering$c_centers = c_centers

  if (graph) {
    # set up to plot as a graph
    message("building a MST representation of clusters...")
    fluster_obj = fluster_add_mst(fluster_obj)
  }

  if (tsne) {
    message("builing a tSNE representation of clusters...")
    fluster_obj = fluster_add_tsne(fluster_obj)
  }

  fluster_obj
}

#' @title Add a minimum spanning tree representation to the fluster object
#' @description Add a minimum spanning tree representation to the fluster object
#' #param fluster_obj The result of running fluster()
#' @return A fluster object with the graph slot populated.
#' @export
fluster_add_mst = function(fluster_obj) {
  mfi = as.list(data.frame(fluster_obj$centers))
  g = build_graph(mfi = mfi)
  g = add_mfi_vertex_attributes(g, mfi)
  mst = igraph::mst(g)
  ag = fluster_obj$agnes_obj
  n_clust = length(fluster_obj$clustering$c_index)
  cag = agnes_to_community(ag, nclust = n_clust)   # BUGBUGBUG - not correct with merging
  gcomm = make_graph_from_community(comm = cag, g = mst)
  gcomm = attach_layout_fr(gcomm)
  fluster_obj$graph = gcomm

  fluster_obj
}

#' @title Add a tSNE representation to the fluster object
#' @description Add a tSNE representation to the fluster object
#' #param fluster_obj The result of running fluster()
#' @return A fluster object with the tsne slot populated.
#' @export
fluster_add_tsne = function(fluster_obj) {
  centers = fluster_obj$clustering$c_centers
  n_items = nrow(centers)
  perplexity = min((n_items - 1) / 3, 30)
  set.seed(137)   # so we'll get the same map for the same data
  res = Rtsne(dist(centers), perplexity = perplexity)$Y
  colnames(res) = c("tsne_1", "tsne_2")
  fluster_obj$tsne = res

  fluster_obj
}


#' @title plot_fluster_graph
#' @description Draw a picture of the result of fluster using graph-based representation of clusters.
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
plot_fluster_graph = function(fluster, markers = colnames(fluster$centers), vs = 10, ms = 5, log.size = FALSE, vertex.frame = TRUE, cex.main = 2, cex.lab = 2) {
  plot_comm_spread(fluster$graph, markers = colnames(fluster$mat), vs = vs, ms = ms,
                   log.size = log.size, vertex.frame = vertex.frame, cex.main = cex.main)
  draw_color_scale(cex.lab = cex.lab)
}

#' @title plot_fluster_tsne
#' @description Draw a picture of the result of fluster using tsne representation of clusters.
#' @param fluster The result of running fluster
#' @param markers Markers to include in the spread
#' @param mode Compute colors using either arithmetic or robust average
#' @param cex Scale factor for node size
#' @param proportional Logical.  Scale by the number of events in the cluster
#' @param emph Logical.  Emphasize each blob with a black line.
#' @param cex.lab Scale factor for titles of the individual markers
#' @param highlight_clusters IF not NULL, a collection of cluster indices to highlight.
#' @param legend Logical.  Whether or not to draw a legend.
#' @return N/A.
#' @examples
#' plot_fluster(fluster_obj)
#' @export
plot_fluster_tsne = function(fluster, markers = colnames(fluster$centers), mode = c("arithmetic", "robust"),
                             cex = 20.0, proportional = TRUE, emph = TRUE, cex.lab = 2,
                             highlight_clusters = NULL, legend = TRUE, show_cluster_numbers = NULL) {
  plot_tsne_spread(fluster, markers, mode, cex, proportional, emph, highlight_clusters, show_cluster_numbers)
  if (legend) {
    if (markers[1] != 'categorical') {
      draw_color_scale(cex.lab = cex.lab)
    } else {
      draw_cluster_legend(fluster_obj = fluster, cex.text = cex.lab)
    }
  }
}

#' title map_functional_names
#' @description Based on a user-specified table, assign symbolic names to clusters
#' @param fluster_obj The result of running fluster()
#' @param defs_file A file containing the functional definitions.
#' @return A decorated fluster object
#' @export
map_functional_names = function(fluster_obj, defs_file) {
  fd = retrieve_categorical_definitions(defs_file)
  fluster_obj = assign_functional_names(fluster_obj, fd)

}

#' @title Map a Sample to a Fluster Model
#' @description This function determines, for a single sample, the number of cells in each cluster.
#' @param ff A sample flowFrame
#' @param fluster_obj An object of type "fluster", the result of running fluster()
#' @return  Per-cluster counts and fractions
#' @export
fluster_map_sample = function(ff, fluster_obj) {
  # apply the flowFPModel to the sample
  fp = flowFP(fcs = ff, model = fluster_obj$mod)

  # get the vector of event bin membership
  btag = tags(fp)[[1]]

  # assign event cluster membership
  nclust = max(fluster_obj$clustering$clst)
  nevents = nrow(ff)
  c_count = vector('numeric', length = nclust)
  for (i in 1:nclust) {
    bidx = fluster_obj$clustering$c_index[[i]]
    eidx = which(btag %in% bidx)
    c_count[i] = length(eidx)
  }

  # convert to percentages of total cells
  c_pctg = c_count / nevents

  # return the result
  return(list(counts = c_count, fractions = c_pctg))
}

#' @title Visualize Cluster Phenotypes
#' @description Draw a "phenobar" representation of a cluster phenotype.  Bars have
#' a height equalt to the medial value of the parameter and are
#' color-coded.  Error flags represent first and third quartiles of the bin centers
#' belonging to the cluster.
#' @param fluster_obj An object of type "fluster", the result of running fluster()
#' @param parameters Which parameters to include in the plot (default = all parameters)
#' @param cluster Which cluster to plot.
#' @param bin_indices Instead of bin indices in a cluster, specify them directly.
#' @param mode Use either arithmetic (mean/sd) or robust (median/quarties).
#' @param plot_global_flag Indicate the global distributions.
#' @export
fluster_phenobars = function(fluster_obj,
                             parameters = colnames(fluster_obj$centers),
                             cluster = 1, bin_indices = NULL,
                             mode = c("arithmetic", "robust"),
                             plot_global_flag = FALSE,
                             main = paste("Cluster", cluster)) {
  mode = match.arg(mode)

  # make an empty plot
  plot(0, 0, pch = '', xlim = c(0, bx(262143)), ylim = c(1 - .3, length(parameters) + .3),
       xaxt = 'n', yaxt = 'n',
       xlab = '', ylab = '',
       main = main)
  wadeTools::ax(1, type = 'biexp')
  axis(side = 2, labels = parameters, at = 1:length(parameters), las = 1)

  centers = fluster_obj$centers

  # get the bin indices of the cluster
  if(is.null(bin_indices)) {
    idx = which(fluster_obj$clustering$clst == cluster)
  } else {
    idx = bin_indices
  }

  val = distributions_bins(fluster_obj, parameters, bin_indices = idx)

  if (mode == "robust") {
    # draw the median
    col = pcolor(val$med, min_value = 0, max_value = 5)
    add_bars(vals = val$med, yvals = 1:length(parameters), col = col)

    # draw the flags
    for (i in 1:length(parameters)) {
      draw_flag(y = i, q1 = val$q1[i], q3 = val$q3[i], med = NA, cex = 2, lwd = 2)
    }

    # add global flags
    if(plot_global_flag) {
      n_bins = 2^nRecursions(fluster_obj$mod)
      qglbl = distributions_bins(fluster_obj, parameters = parameters, bin_indices = 1:n_bins)
      for (i in 1:length(parameters)) {
        draw_flag(y = i-.2, q1 = qglbl$q1[i], q3 = qglbl$q3[i], med = qglbl$med[i], cex = 2, lwd = 2, col = 'gray')
      }
    }
  } else {
    # draw the mean
    col = pcolor(val$med, min_value = 0, max_value = 5)
    add_bars(vals = val$mn, yvals = 1:length(parameters), col = col)

    # draw the flags
    for (i in 1:length(parameters)) {
      draw_flag(y = i, q1 = val$lo[i], q3 = val$hi[i], med = NA, cex = 2, lwd = 2)
    }

    # add global flags
    if(plot_global_flag) {
      n_bins = 2^nRecursions(fluster_obj$mod)
      qglbl = distributions_bins(fluster_obj, parameters = parameters, bin_indices = 1:n_bins)
      for (i in 1:length(parameters)) {
        draw_flag(y = i-.2, q1 = qglbl$lo[i], q3 = qglbl$hi[i], med = qglbl$mn[i], cex = 2, lwd = 2, col = 'gray')
      }
    }
  }
}


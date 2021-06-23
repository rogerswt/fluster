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

#' @title Fingerprint-based Clustering
#' @param fcs The data (either a flowFrame or a flowSet)
#' @param parameters The parameters in fcs used for analysis.  Default is all
#' parameters in the input data.
#' @param nRecursions The number of recursions in calculating the fingerprint (default = 12)
#' @param nclust The number of clusters you want panoply to make.  If NULL, fluster
#' will decide based on analysis of the clustering dendrogram.
#' @param merge Logical: should we merge initial clusters based on categorical similarity?
#' will make a guess as to the "best" number of clusters
#' @param graph Logical: should we compute the MST for visualization?
#' @param tsne Logical: should we compute a tSNE embedding for visualization?
#' @param manual_thresholds A named vector of one or more values that will override
#' the internally-calculated positivity thresholds.  If NULL (default) use the internally-calculated
#' thresholds.
#' @param modality If previously computed with examine_positivity_thresholds, you can pass
#' in the resulting modality object.  If NULL, modality is computed internally.
#' @param sd_fac A factor by which to multiply per-marker cluster standard deviation to determine
#' positivity.  Smaller values tend to force a marker to be declared either positive
#' or negative, rather than "undetermined".  Default: 1.0.
#' @description Fluster (**F**ingerprint-based c**luster**ing)
#'  implements a workflow starting with Cytometric Fingerprint (CF) binning of
#' a flowFrame (or flowSet)
#' using flowFP to a relatively high resolution (default nRecursions of 12 results
#' in 4096 bins).  The bin centroids are then computed, using the median value of all
#' of the events in the bin for each of the included parameters.  Next, these
#' multivariate bin centroids are clustered using agglommerative hierarchical clustering
#' \code{cluster::agnes()}.  Bins in a cluster carry their cellular cargo into
#' the cluster.  The resulting data are represented in a graph structure and/or as a tSNE embedding
#' for visualization and interpretation.
#'
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
fluster = function(fcs, parameters = NULL, nRecursions = 12, nclust = NULL, merge = TRUE, graph = TRUE, tsne = TRUE,
                   manual_thresholds = NULL, modality = NULL, sd_fac = 1.0) {

  call.args = list(parameters = parameters, nRecursions = nRecursions,
                   nclust = nclust, merge = merge, graph = graph, tsne = tsne,
                   manual_thresholds = manual_thresholds,
                   modality = modality, sd_fac = sd_fac
  )

  if (is(fcs, "flowFrame")) {
    ff = fcs
  } else if (is(fcs, "flowSet")) {
    message("Aggregating the flowSet...")
    ff = suppressWarnings(as(fcs, "flowFrame"))
    flowCore::exprs(ff) = flowCore::exprs(ff)[,which(flowCore::colnames(flowCore::exprs(ff)) != "Original")]
  } else {
    stop("Argument fcs must either be a flowFrame or a flowSet\n")
  }
  # check parameters
  if (is.null(parameters)) {
    parameters = flowCore::colnames(ff)
  }
  if (is.numeric(parameters)) {
    parameters = flowCore::colnames(ff)[parameters]
  }

  message("computing fingerprint bins...")
  mod = flowFP::flowFPModel(ff, parameters = parameters, nRecursions = nRecursions)
  fp = flowFP::flowFP(ff, mod)

  message("calculating bin centers...")
  res = calculate_bin_phenotypes(fp = fp, fs = ff)
  mat = t(res$center)
  variance = t(res$variance)

  # agnes on bin centers
  message("clustering bins...")
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

  fluster_obj = list(call.args = call.args, fcs = ff, parameters = parameters, mod = mod, fp = fp, centers = mat, bvar = variance, agnes_obj = ag,
                     graph = NULL, tsne = NULL, clustering = clst, modality = NULL)
  class(fluster_obj) = "fluster"

  # determining modality and calculating positivity thresholds
  if (is.null(modality)) {
    message("calculating positivity thresholds using cluster extremes...")
    fluster_obj = positivity_thresholds(fluster_obj, ff, manual_thresholds = manual_thresholds)
  } else {
    message("Using previously calculated thresholds")
  }

  # merging clusters
  if (merge) {
    fluster_obj = merge_categorical_clusters(fluster_obj = fluster_obj, sd_fac = sd_fac)
    nclust = max(fluster_obj$clustering$clst)
    message("merging clusters, resulting in ", nclust, " clusters...")

  }

  # calculate cluster centers, mostly for visualization
  message("Calculating cluster centers...")
  c_centers = matrix(NA, nrow = 0, ncol = length(parameters))
  colnames(c_centers) = parameters
  for (i in 1:nclust) {
    idx = which(fluster_obj$clustering$clst == i)
    c_centers = rbind(c_centers, distributions_bins(fluster_obj, bin_indices = idx)$mn)
  }
  fluster_obj$clustering$c_centers = c_centers

  if (graph) {
    # set up to plot as a graph
    message("building a MST representation of clusters...")
    fluster_obj = fluster_add_mst(fluster_obj)
  }

  if (tsne) {
    message("building a tSNE representation of clusters...")
    fluster_obj = fluster_add_tsne(fluster_obj)
  }

  fluster_obj
}

#' @title Make a Plot of Modality
#' @description Positivity thresholds are calculated by finding the compact extremes
#' of marker expressions in various clusters.  Thresholds are the mid-point between
#' extremes.
#' @param fluster_obj A fluster object.
#' @export
display_modality = function(fluster_obj) {
  modality = fluster_obj$modality
  parameters = names(modality$thresh)
  kde = modality$kde

  # for displaying as polygons make sure first and last points of the kde's are 0
  klen = length(kde[[1]][[1]]$y)
  for (marker in parameters) {
    kde[[marker]][["global"]]$y[1] = kde[[marker]][["global"]]$y[klen] = 0
    kde[[marker]][["lo"]]$y[1] = kde[[marker]][["lo"]]$y[klen] = 0
    kde[[marker]][["hi"]]$y[1] = kde[[marker]][["hi"]]$y[klen] = 0
  }

  # calculate plot layout
  n = length(parameters) + 1
  sq = sqrt(n)
  frac = sq - floor(sq)
  if (frac == 0) {
    ac = dn = floor(sq)
  } else {
    ac = floor(sq) + 1
    dn = ceiling(n / ac)
  }

  opar = par(mfrow = c(dn, ac), mar = c(2, 2, 2, 1))
  for (marker in parameters) {
    plot(kde[[marker]][["global"]], type = 'l', lwd = 1, xaxt = 'n', xlab = '', ylab = '', main = marker)
    ax()
    bcol = "dodgerblue2"
    rcol = "indianred2"
    polygon(kde[[marker]][["global"]], col = make_transparent("black", alpha = .25), border = "black")
    polygon(kde[[marker]][["lo"]], col = make_transparent(bcol), border = bcol)
    polygon(kde[[marker]][["hi"]], col = make_transparent(rcol), border = rcol)
    xline(modality$med_lo[marker], lty = 'dotdash', col = "dodgerblue2")
    xline(modality$med_hi[marker], lty = 'dotdash', col = "indianred2")
    xline(modality$thresh[marker], lty = 'dotdash', col = "black", lwd = 2)
  }
  par(opar)
}

#' @title Merge Categorically Similar Clusters
#' @description Clusters are labeled with a categorical vector in which each
#' marker is either "hi" or "lo" with respect to a threshold.  If a marker is not unambiguously
#' either hi or lo, it's labeled as "un" for "unknown.  To receive hi (lo), the
#' cluster center must be sufficiently above (below) the threshold in units of
#' the standard deviation of that marker.
#' @param fluster_obj A fluster object.
#' @param sd_fac A factor multiplying the standard deviation to determine if that
#' marker is sufficiently above (below) the threshold in order to labeled "hi" ("lo").
#' @return A fluster object after categorical merging.
#' @export
merge_categorical_clusters = function(fluster_obj, sd_fac = 1.0) {
  n_clust = max(fluster_obj$clustering$clst)

  # overwrite affected call parameters in fluster_obj
  fluster_obj$call.args$sd_fac = sd_fac
  fluster_obj$call.args$merge = TRUE

  # get the categorical mapping
  categ = list()
  # for (i in 1:n_clust) {
    categ = categorical_phenotype_all_clusters(fluster_obj = fluster_obj, sd_fac = sd_fac)
  # }

  # roll through and create clusters of clusters
  cmerge = list()
  phenotype = list()
  # cvec is a vector of cluster indices.  When a cluster joins a merge, it's removed from this vector
  cvec = 1:n_clust
  k = 1

  while (length(cvec) > 0) {
    # get the head of the list of remaining clusters
    ith = cvec[1]
    cmerge[[k]] = ith                          # add ith to the next merge
    phenotype[[k]] = categ[[ith]]              # record the phenotype
    cvec = cvec[which(cvec != ith)]            # remove it from cvec
    j = 1
    while (j <= length(cvec)) {
      jth = cvec[j]
      if (compare_categories(categ[[ith]], categ[[jth]])) {
        cmerge[[k]] = append(cmerge[[k]], jth)     # add jth cluster to cmerge
        cvec = cvec[which(cvec != jth)]            # remove jth cluster from cvec
        j = j - 1                                  # don't skip next element
      }
      j = j + 1
    }
    k = k + 1
  }

  # replace clustering slot with the merged result
  if (is.null(fluster_obj[["original_clustering"]])) {
    orig_clustering = fluster_obj$clustering
  }
  c_index = list()
  for (i in 1:length(cmerge)) {
    c_index[[i]] = vector(mode = 'numeric')
    for (j in 1:length(cmerge[[i]])) {
      c_index[[i]] = append(c_index[[i]], orig_clustering$c_index[[cmerge[[i]][j]]])
    }
  }
  nbins = 2 ^ nRecursions(fluster_obj$mod)
  clst = rep(NA, length = nbins)
  for (i in 1:length(c_index)) {
    clst[c_index[[i]]] = i
  }
  clustering = list(clst = clst, c_index = c_index, phenotype = phenotype, func_phenotype = NULL)
  fluster_obj$orig_clustering = orig_clustering
  fluster_obj$clustering = clustering

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
  cag = agnes_to_community(ag, nclust = n_clust)
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
#' @param plot_global_flag Indicate the global distributions.
#' @param show_thresholds Logical. Show per-parameter thresholds.
#' @param show_sd_fac Logical.  Superimpose modified error flags if sd_fac != 1.0.
#' @param main Title of plot.
#' @export
fluster_phenobars = function(fluster_obj,
                             parameters = fluster_obj$parameters,
                             cluster = 1, bin_indices = NULL,
                             plot_global_flag = FALSE,
                             show_thresholds = TRUE,
                             show_sd_fac = TRUE,
                             main = paste("Cluster", cluster)) {

  # make an empty plot
  plot(0, 0, pch = '', xlim = c(0, bx(262143)), ylim = c(1 - .3, length(parameters) + .3),
       xaxt = 'n', yaxt = 'n',
       xlab = '', ylab = '',
       main = main)
  wadeTools::ax(1, type = 'biexp')
  axis(side = 2, labels = parameters, at = 1:length(parameters), las = 1)

  centers = fluster_obj$centers

  # get the bin indices of the cluster
  if (is.null(bin_indices)) {
    idx = which(fluster_obj$clustering$clst == cluster)
  } else {
    idx = bin_indices
  }

  val = distributions_bins(fluster_obj, bin_indices = idx)

  # draw the mean
  col = pcolor(val$mn, min_value = 0, max_value = 5)
  add_bars(vals = val$mn, yvals = 1:length(parameters), col = col)

  # draw the flags
  for (i in 1:length(parameters)) {
    draw_flag(y = i, q1 = val$lo[i], q3 = val$hi[i], med = NA, cex = 2, lwd = 2)
  }

  if (show_sd_fac) {
    sd_fac = fluster_obj$call.args$sd_fac
    if (sd_fac != 1) {
      for (i in 1:length(parameters)) {
        sdev = val$hi[i] - val$lo[i]
        mod_sdev = sd_fac * sdev
        draw_flag(y = i, q1 = val$mn[i] - .5 * mod_sdev, q3 = val$mn[i] + 0.5 * mod_sdev, med = NA, cex = 2, lwd = 2, col = 'gray')
      }
    }
  }

  if (show_thresholds) {
    for (i in 1:length(parameters)) {
      p = parameters[i]
      draw_thresh(y = i, thresh = fluster_obj$modality$thresh[p], len = .7, col = 'blue', lwd = 2)
    }
  }

  # add global flags
  if (plot_global_flag) {
    n_bins = 2^nRecursions(fluster_obj$mod)
    qglbl = distributions_bins(fluster_obj, bin_indices = 1:n_bins)
    for (i in 1:length(parameters)) {
      draw_flag(y = i - .2, q1 = qglbl$lo[i], q3 = qglbl$hi[i], med = qglbl$mn[i], cex = 2, lwd = 2, col = 'gray')
    }
  }
}

#' @title Gate Events in Cluster(s)
#' @description Given a fluster model and a list of its clusters, gate a flowFrame
#' or flowSet so as to include events that below to these clusters.
#' @param fluster_obj An object of type "fluster", the result of running fluster()
#' @param fcs_obj A flowFrame or flowSet to be gated.
#' @param clusters A list of clusters we want to gate.
#' @return The gated flowFrame or flowSet
#' @export
fluster_gate_clusters = function(fluster_obj, fcs_obj, clusters) {
  # get a list of bins that belong to the indicated clusters
  bvec = which(fluster_obj$clustering$clst %in% clusters)

  # gate the fcs_obj
  if (is(fcs_obj) == "flowFrame") {
    ff = fcs_obj

    # apply the flowFP model to the fcs_obj
    ff = fcs_obj
    fp = flowFP(ff, fluster_obj$mod)
    ev = which(tags(fp)[[1]] %in% bvec)
    exprs(ff) = exprs(ff)[ev, ]
    return(ff)

  } else {
    fs = fcs_obj
    fp = flowFP(fs, fluster_obj$mod)
    for (i in 1:length(fs)) {
      ev = which(tags(fp)[[i]] %in% bvec)
      exprs(fs[[i]]) = exprs(fs[[i]])[ev, ]
    }
    return(fs)
  }
}


# New method to find positivity thresholds.  Given a (high resolution) clustering,
# calculate, for each marker, the median and spread of each cluster.  Then, identify
# the most negative and most positive (and also tight) clusters.  The positivity
# threshold will be the mid-point between them.
#
# We prioritize tight clusters over diffuse clusters by computing a figure of merit,
# which is the normalized median value divided by the spread.
#
# NOTE: One of fluster_obj or ff must not be null.  If fluster_obj is null, then
# it is calculated from ff.  If fluster_obj is not null, and its internal fcs
# object also isn't null, it is used as ff, and the ff argument is ignored.






#' @title  Positivity Thresholds
#' @param fluster_obj An existing fluster object, the result of running fluster().
#' @param fcs A flowFrame or flowSet.
#' @param parameters The parameters you wish to consider.
#' @param manual_thresholds A named vector of one or more values that will override
#' the internally-calculated positivity thresholds.
#' @param show Logical.  Should we display the result graphically?
#' @param show_range The range of values to include in calculating the kernel
#' density estimates.
#' @description Given a (high resolution) clustering, calculate, for each marker,
#' the median and spread of each cluster.  Then, identify the most negative and most
#' positive (and also tight) clusters.  The positivity threshold will be between
#' them, closer to the tighter of the two.
#'
#' We prioritize tight clusters over diffuse clusters by computing a figure of merit,
#' which is the normalized median value divided by the spread.
#'
#' NOTE: One of fluster_obj or ff must not be null.  If fluster_obj is null, then
#' it is calculated from ff.  If fluster_obj is not null, and its internal fcs
#' object also isn't null, it is used as ff, and the ff argument is ignored.
#'
#' ALSO NOTE: This function internally de-rails and de-negs the data so that the
#' kernel density estimates are not confused by rail or excessivley negative events.
#' This is an experimental feature.
#'
#' @return A fluster object that includes a modality slot.
#'
#' @export
positivity_thresholds = function(fluster_obj = NULL, fcs = NULL,
                                 parameters = NULL,
                                 manual_thresholds = NULL,
                                 show = FALSE, show_range = NULL) {
  if (is.null(fluster_obj)) {
    if (is.null(fcs)) {
      stop("fluster_obj and fcs can't both be null.")
    }

    if (is(fcs, "flowSet")) {
      fcs = suppressWarnings(as(fcs, "flowFrame"))
      flowCore::exprs(fcs) = flowCore::exprs(fcs)[,which(flowCore::colnames(flowCore::exprs(fcs)) != "Original")]
    }
    if (is.null(parameters)) {
      stop("Need to supply parameters if we're starting with FCS data (and not a fluster object).")
    }
    if (is.numeric(parameters)) {
      parameters = flowCore::colnames(fcs)[parameters]
    }
    fluster_obj = fluster(fcs = fcs, parameters = parameters, merge = FALSE, graph = FALSE, tsne = FALSE)
  }

  # sanity check that names of thresholds matches names in parameters
  if (!is.null(manual_thresholds)) {
    n_manual = length(manual_thresholds)
    if (length(which(names(manual_thresholds) %in% parameters)) < n_manual) {
      message("WARNING: names of manual thresholds do not match parameters.  Ignoring...")
      manual_thresholds = NULL
    }
  }

  ff = fluster_obj$fcs
  parameters = fluster_obj$parameters

  # # experiment with derail and deneg to avoid misleading threshold determination
  # ff = derail(ff, parameters = parameters)
  # ff = deneg(ff, parameters = parameters, min.value = bx(-1000))
  # fp = flowFP(ff, model = fluster_obj$mod)
  # # fp = fluster_obj$fp
  #
  # tg = tags(fp)[[1]]
  # tmp_flus = fluster_obj
  # tmp_flus$fcs = ff
  # tmp_flus$fp = fp
  #
  # res = spreads_and_meds(tmp_flus)

  tg = tags(fluster_obj$fp)[[1]]
  res = spreads_and_meds(fluster_obj)

  kde = list()    # for visualization
  if (is.null(show_range)) {
    idx = which(flowCore::colnames(ff) %in% parameters)
    bot = min(parameters(ff)$minRange[idx])
    top = max(parameters(ff)$maxRange[idx])
    show_range = c(bot, top)
  }
  med_lo = med_hi = thresh = rep(NA, length = length(parameters))
  names(med_lo) = names(med_hi) = names(thresh) = parameters
  for (marker in parameters) {
    idx = select_extremes(med = res$med[, marker], spread = res$spread[, marker])
    kde[[marker]] = list()
    kde[[marker]][["global"]] = normalize.kde(bkde(exprs(ff)[, marker], gridsize = 1001, range.x = show_range))
    ev_lo = which(tg %in% fluster_obj$clustering$c_index[[idx$lo]])
    kde[[marker]][["lo"]] = normalize.kde(bkde(exprs(ff)[ev_lo, marker], gridsize = 1001, range.x = show_range))
    wid_lo = IQR(exprs(ff)[ev_lo, marker])
    ev_hi = which(tg %in% fluster_obj$clustering$c_index[[idx$hi]])
    kde[[marker]][["hi"]] = normalize.kde(bkde(exprs(ff)[ev_hi, marker], gridsize = 1001, range.x = show_range))
    wid_hi = IQR(exprs(ff)[ev_hi, marker])

    med_lo[marker] = median(exprs(ff)[ev_lo, marker])
    med_hi[marker] = median(exprs(ff)[ev_hi, marker])
    pos_fac = wid_lo / (wid_lo + wid_hi)
    thresh[marker] = med_lo[marker] + pos_fac * (med_hi[marker] - med_lo[marker])
  }

  modality = list(kde = kde, med_lo = med_lo, med_hi = med_hi, thresh = thresh)

  # attach modality to the fluster object
  fluster_obj$modality = modality

  if (show) {
    display_modality(fluster_obj)
  }

  # override manually supplied thresholds
  for (i in 1:length(manual_thresholds)) {
    fluster_obj$modality$thresh[names(manual_thresholds)[i]] = manual_thresholds[i]
  }

  return(fluster_obj)
}



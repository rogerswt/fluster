#
# fluster_utils.R
#
# declares utility functions, not intended to be exposed to the user.
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

# input is a flowFP (fp) and the corresponding flowSet (fs)
calculate_bin_phenotypes = function(fp, fs) {
  parameters = parameters(fp)
  n_bins = 2 ^ nRecursions(fp)
  n_parameters = length(parameters)
  center = matrix(NA, nrow = n_parameters, ncol = n_bins)
  rownames(center) = parameters
  variance = matrix(NA, nrow = n_parameters, ncol = n_bins)
  rownames(variance) = parameters

  # for convenience, lump the frames in fs together and recalculate the fingerprint
  ff = as(fs, "flowFrame")
  fp = flowFP(fcs = ff, model = as(fp, "flowFPModel"))
  for (i in 1:n_bins) {
    idx = which(tags(fp)[[1]] == i)
    if (length(idx) == 0) {
      next
    }

    for (j in 1:n_parameters) {
      p = parameters[j]
      vals = exprs(ff)[idx, p]

      # aiming to limit huge error bars due to excessive negative spread
      lowest_possible_mfi = bx(-1000)
      vidx = which(vals > lowest_possible_mfi)
      if (length(vidx) > 0) {
        center[j, i] = median(vals[vidx], na.rm = TRUE)
        variance[j, i] =  var(vals[vidx], na.rm = TRUE)   # hack
      } else {
        center[j, i] = lowest_possible_mfi
        variance[j, i] = 0
      }
    }
  }
  return(list(center = center, variance = variance))
}

# assumes biexp vert scale
draw_color_scale = function(min_col_value = 0, max_col_value = 5, transformation, ...) {
  requireNamespace("wadeTools")
  ll = -0.5
  ul = bx(262143)

  vec = seq(ll, ul, length.out = 500)
  cols = pcolor(pvalue = vec, min_value = min_col_value, max_value = max_col_value)

  opar = par(mar = c(0, 5, 0, 0) + .1)
  plot(0, 0, pch = '', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "Signal",
       xlim = c(0, 5), ylim = c(ll, ul), ...
  )
  for (i in 1:length(vec)) {
    y = vec[i]
    segments(x0 = 0, y0 = y, x1 = 1, y1 = y, col = cols[i], lwd = 3)
  }
  wadeTools::ax(axis = 2, type = transformation, ...)
  par(opar)
}

draw_per_marker_scale = function(dyn_range = 2, ...) {
  ll = -dyn_range
  ul = dyn_range

  vec = seq(ll, ul, length.out = 500)
  cols = red_blue_fun(vec, dyn_range = dyn_range)

  opar = par(mar = c(0, 5, 0, 0) + .1)
  plot(0, 0, pch = '', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "Rel. Signal",
       xlim = c(0, 5), ylim = c(ll, ul), ...
  )
  for (i in 1:length(vec)) {
    y = vec[i]
    segments(x0 = 0, y0 = y, x1 = 1, y1 = y, col = cols[i], lwd = 3)
  }
  at = -dyn_range:dyn_range
  lab = at
  lab[which(at == 0)] = "Thresh"
  axis(side = 2, at = at, labels = lab)
  par(opar)
}

# draw legend for categorical labeling
draw_cluster_legend = function(fluster_obj, cex.text = 1.25) {
  opar = par(mar = c(.1, .1, .1, .1))
  plot(0, 0, pch = '', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xlim = c(0, 1), ylim = c(0, 1), bty = 'n')
  xcol = 0.05
  wcol = .1
  xtext = 0.18
  subsets = fluster_obj$labels$names
  cols = fluster_obj$labels$colors

  n_items = length(subsets)
  y = seq(.05, .95, length.out = n_items)
  ht = .8 * .5 / n_items
  for (i in 1:n_items) {
    rect(xleft = xcol, ybottom = y[i] - ht, xright = xcol + wcol, ytop = y[i] + ht, col = cols[i])
    text(x = xtext, y = y[i], labels = subsets[i], pos = 4, cex = cex.text)
  }
  par(opar)
}

# define a color function to transform a parameter value
pcolor = function(pvalue, min_value = 0, max_value = 4) {
  top_col = 'red'
  bot_col = 'darkgreen'
  mid_col = 'yellow'
  zero_col = 'darkgray'
  zero_col = bot_col

  len = length(pvalue)
  if (length(which(pvalue < min_value)) == len) {
    col_values = rep(zero_col, len)
  } else if (length(which(pvalue > max_value)) == len) {
    col_values = rep(top_col, len)
  } else {
    pvalue[pvalue < min_value] = min_value
    pvalue[pvalue > max_value] = max_value

    col_values = fields::color.scale(
      z = pvalue, col = fields::two.colors(
        n = 100,
        start = bot_col,
        end = top_col,
        middle = mid_col),
      zlim = c(min_value, max_value)
    )
    col_values[which(pvalue <= min_value)] = zero_col
  }

  col_values
}

# standardize cluster centers by subtracting thresh, dividing by sdev
standardize_to_threshold = function(vals, thresh) {
  sdev = sd(vals)
  vals = (vals - thresh) / sdev

  vals
}

red_blue_fun = function(vals, dyn_range = 2) {
  top_col = "darkred"
  pos_mid = "pink"
  bot_col = "darkblue"
  neg_mid = "lightblue"
  mid_col = "white"

  len = length(vals)

  pos_idx = which(vals >= 0)
  neg_idx = which(vals <= 0)

  vals[vals >  dyn_range] =  dyn_range
  vals[vals < -dyn_range] = -dyn_range
  pos_cols = fields::color.scale(z = vals[pos_idx], fields::two.colors(n = 51, start = mid_col, end = top_col, middle = pos_mid), zlim = c(0, dyn_range))
  neg_cols = fields::color.scale(z = vals[neg_idx], fields::two.colors(n = 51, start = bot_col, end = mid_col, middle = neg_mid), zlim = c(-dyn_range, 0))

  col_values = rep(NA, len)
  col_values[pos_idx] = pos_cols
  col_values[neg_idx] = neg_cols

  col_values
}

 # calculate red/blue colors using per-marker order statistics
perMarkerColor = function(vals, thresh) {

  vals = standardize_to_threshold(vals, thresh)

  col_values = red_blue_fun(vals = vals)

  col_values
}

add_bars = function(vals, yvals, col) {
  hw = 0.4
  if (length(col) == 1) {
    col = rep(col, length(vals))
  }
  for (i in 1:length(vals)) {
    rect(xleft = 0, ybottom = yvals[i] - hw, xright = vals[i], ytop = yvals[i] + hw, col = col[i], border = 'black')
  }
}

draw_flag = function(y, q1, q3, med = NA, ...) {
  segments(x0 = q1, y0 = y, x1 = q3, y1 = y, ...)
  if (!is.na(med)) {points(med, y, pch = 20, ...)}
}

draw_thresh = function(y, thresh, len, ...) {
  segments(x0 = thresh, y0 = y - .5 * len, x1 = thresh, y1 = y + .5 * len, ...)
}


################################################################################
################################################################################
# here starts igraph stuff
################################################################################
################################################################################
build_graph = function(mfi) {
  # make mfi into a distance matrix
  mat = matrix(NA, nrow = length(mfi[[1]]), ncol = 0)

  for (i in 1:length(mfi)) {
    mat = cbind(mat, mfi[[i]])
  }
  colnames(mat) = names(mfi)
  dst = as.matrix(stats::dist(mat))
  g = graph_from_adjacency_matrix(adjmatrix = dst, mode = 'undirected', weighted = TRUE)

  # some algorithms interpret weight as distance, and others as strength.  Let's
  # preserve both interpretations so they can be conveniently swapped.  Initial weight
  # will be distance.  This is appropriate for MST.
  # modularity treats weights as strengths
  # betweenness treats weights as distances
  dstnce = E(g)$weight
  strngt = 1 / dstnce
  edge_attr(g, 'distance') <- dstnce
  edge_attr(g, 'strength') <- strngt

  g
}

# create a new graph from a communities object.  g is the MST graph from which the
# communities object 'comm' was induced.
# note:  marker expressions should already have been added to g.
make_graph_from_community = function(comm, g) {
  markers = vertex_attr_names(g)
  markers = markers[-which(markers == "name")]

  n_vertices = max(comm$membership)
  gcomm = make_empty_graph(n = n_vertices, directed = FALSE)
  sizes = vector('numeric')
  indices = list()
  for (i in 1:n_vertices) {
    indices[[i]] = idx = which(comm$membership == i)
    V(gcomm)[i]$size = length(indices[[i]])
    for (j in 1:length(markers)) {
      vertex_attr(gcomm, markers[j])[i] <- mean(vertex_attr(g, markers[j])[idx])
    }
  }

  # add edges.  Edges are the cummulative strengths between communities.
  k = 1
  for (i in 1:(n_vertices - 1)) {
    for (j in (i + 1):n_vertices) {
      res = edges_between_communities(comm, g, indices[[i]], indices[[j]])
      if (res$strength != 0) {
        # add an edge
        gcomm = add_edges(gcomm, c(i, j))
        E(gcomm)$strength[k] = res$strength
        k = k + 1
        # cat(i, j, "\n")
      }
    }
  }
  E(gcomm)$distance = 1/E(gcomm)$strength
  E(gcomm)$weight = E(gcomm)$strength
  mst(gcomm)
}

# c1 is the vector of indices of community 1, c2 of community 2
edges_between_communities = function(comm, g, c1, c2) {
  edges = E(g)[c1 %--% c2]
  strength = sum(E(g)$strength[edges])

  return(list(edges = edges, strength = strength))
}

# if the graph is induced from a distance matrix, then the edge weights are
# equal to distances.  However, the strongest edges should be the ones closest
# in metric distance, so invert the weights in this case.
set_weight_as = function(g, weight = c("distance", "strength")) {
  weight = match.arg(weight)
  E(g)$weight = edge_attr(g, weight)

  g
}

# best so far
attach_layout_fr = function(g, overwrite = TRUE) {
  set.seed(137)
  g = add_layout_(g, with_fr(), overwrite = overwrite)

  g
}

# also nice
attach_layout_lgl = function(g, overwrite = TRUE) {
  set.seed(137)
  g = add_layout_(g, with_lgl(), overwrite = overwrite)

  g
}

attach_layout = function(g, layfun, overwrite = TRUE) {
  set.seed(137)
  g = add_layout_(g, layfun, overwrite = overwrite)

  g
}

add_mfi_vertex_attributes = function(g, mfi) {
  # add marker expression as vertex attributes
  markers = names(mfi)
  for (i in 1:length(markers)) {
    vertex_attr(g, markers[i]) <- mfi[[i]]
  }

  g
}



plot_community_graph = function(fluster_obj, marker, mode = c("global", "per-marker"), vs = 3, ms = 0, log.size = TRUE, vertex.frame = TRUE, cex.main) {
  if (is.null(fluster_obj$graph)) {
    stop("You must first compute the MST graph using fluster_add_graph")
  }
  g = fluster_obj$graph
  mode = match.arg(mode)
  n_clust = length(fluster_obj$clustering$c_index)
  if (marker == "categorical") {
    cols = fluster_obj$clustering$func_color
  } else {
    if (mode == 'global') {
      cols = pcolor(vertex_attr(g, marker))
    } else {
      cols = perMarkerColor(vertex_attr(g, marker), thresh = fluster_obj$modality$thresh[marker])
    }
  }

  g = set_weight_as(g, "distance")

  if (log.size) {
    vsize = vs * log10(V(g)$size)
    vsize[vsize < ms] = ms
  } else {
    vsize = V(g)$size
    mx = max(vsize)
    mn = min(vsize)
    vsize = vs * vsize / mx
    if (ms != 0) {
      vsize[vsize < ms] = ms
    }
  }
  if (vertex.frame) {
    vfc = 'black'
  } else {
    vfc = NA
  }
  plot(g, vertex.size = vsize, vertex.color = cols,
       vertex.frame.color = vfc, vertex.label = NA)
  title(main = marker, cex.main = cex.main)
}


plot_comm_spread = function(fluster_obj, mode = c("global", "per-marker"), markers = NULL, vs = 3, ms = 1, log.size = TRUE, vertex.frame = FALSE, cex.main) {
  g = fluster_obj$graph
  if (is.null(markers)) {
    markers = vertex_attr_names(g)
    markers = markers[-which(markers == "size")]
  }

  # calculate plot layout
  n = length(markers) + 1
  sq = sqrt(n)
  frac = sq - floor(sq)
  if (frac == 0) {
    ac = dn = floor(sq)
  } else {
    ac = floor(sq) + 1
    dn = ceiling(n / ac)
  }

  par(mfrow = c(dn, ac), mar = c(0, 0, 2, 0))
  for (i in 1:length(markers)) {
    marker = markers[i]
    plot_community_graph(fluster_obj, marker, mode = mode, vs = vs, ms = ms, log.size = log.size, vertex.frame = vertex.frame, cex.main)
  }
}

plot_tsne = function(fluster_obj, marker, mode = c("global", "per-marker"),
                     box = TRUE, cex = 50.0, proportional = TRUE, emph = TRUE,
                     highlight_clusters = NULL, show_cluster_number = NULL) {
  if (is.null(fluster_obj$tsne)) {
    stop("You must first compute the tSNE embedding using fluster_add_tsne")
  }
  mode = match.arg(mode)
  n_clust = length(fluster_obj$clustering$c_index)
  if (marker == "categorical") {
    cols = fluster_obj$clustering$func_color
  } else {
    if (mode == 'global') {
      cols = pcolor(fluster_obj$clustering$c_centers[, marker], min_value = 0.0, max_value = 5.0)
    } else {
      vals = fluster_obj$clustering$c_centers[, marker]
      thresh = fluster_obj$modality$thresh[marker]
      cols = perMarkerColor(vals, thresh)
    }
  }

  map = fluster_obj$tsne
  size = rep(cex, n_clust)
  if (proportional) {
    for (i in 1:n_clust) {
      size[i] = length(fluster_obj$clustering$c_index[[i]])
    }
    size = sqrt(size / (2 ^ nRecursions(fluster_obj$mod)))
  }
  cex = cex * size
  cexbg = 1.05 * cex


  if (!is.null(highlight_clusters)) {
    hcol = rep('black', length = n_clust)
    hcol[highlight_clusters] = 'magenta'
    cexbg[highlight_clusters] = 2.0 * cexbg[highlight_clusters]
  } else {
    hcol = rep('black', length = n_clust)
  }

  # plot largest first
  srt = sort(size, decreasing = TRUE, index.return = TRUE)$ix
  map = map[srt, ]
  cols = cols[srt]
  cex = cex[srt]
  cexbg = cexbg[srt]
  hcol = hcol[srt]

  bty = ifelse(box, 'o', 'n')
  if (emph) {
    xlim = c(min(map[, 1]), max(map[, 1]))
    ylim = c(min(map[, 2]), max(map[, 2]))
    plot(0, 0, pch = '', , bty = bty, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', xlim = xlim, ylim = ylim, main = marker)
    for (i in 1:nrow(map)) {
      points(x = map[i, 1], y = map[i, 2], pch = 20, col = hcol[i], cex = cexbg[i])
      points(x = map[i, 1], y = map[i, 2], pch = 20, col = cols[i], cex = cex[i])
      if (!is.null(show_cluster_number)) {
        text(x = map[i, 1], y = map[i, 2], labels = srt[i], adj = 0.5, cex = show_cluster_number)
      }
    }
  } else {
    plot(map, bty = bty, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', pch = 20, col = cols, cex = cex, main = marker)
  }
}

plot_tsne_spread = function(fluster_obj, markers = NULL, mode = c("global", "per-marker"), cex = 50.0, proportional = TRUE, emph = TRUE, highlight_clusters, show_cluster_number) {
  mode = match.arg(mode)

  # calculate plot layout
  n = length(markers) + 1
  sq = sqrt(n)
  frac = sq - floor(sq)
  if (frac == 0) {
    ac = dn = floor(sq)
  } else {
    ac = floor(sq) + 1
    dn = ceiling(n / ac)
  }

  par(mfrow = c(dn, ac), mar = c(1, 1, 2, 1))
  for (i in 1:length(markers)) {
    plot_tsne(fluster_obj = fluster_obj, marker = markers[i], mode = mode, cex = cex, proportional = proportional, emph = emph, highlight_clusters = highlight_clusters, show_cluster_number = show_cluster_number)
  }
}

################################################################################
################################################################################
# conventional clustering analyses
################################################################################
################################################################################

# ag is an agnes object
# look for a change in slope of nclust vs cut height
# breaks is the number of height breaks
advise_n_clust = function(ag, breaks = 500, show = FALSE, ...) {
  obj = as.hclust(ag)
  ht = seq(min(ag$height), max(ag$height), length.out = breaks)
  nclust = vector('numeric')
  for (i in 1:breaks) {
    nclust[i] = max(cutree(obj, h = ht[i]))
  }
  df = data.frame(nclust = nclust, height = ht)
  df = df[breaks:1, ]

  # eliminate redundant heights
  nclust = unique(df$nclust)
  ht = vector('numeric')
  for (i in 1:length(nclust)) {
    ht[i] = min(df$height[df$nclust == nclust[i]])
  }
  df2 = data.frame(nclust = nclust, height = ht)


  # fit the asymptotic behavior, and look for where it starts to diverge

  # find a pretty flat region and fit it
  kde = df2
  colnames(kde) = c("x", "y")
  der = deriv1.kde(kde)
  crit_slope = .001
  idx = which(abs(der$y) < crit_slope)
  ignore_last = nrow(df2) - max(idx)
  npts = length(idx)
  # ignore_last = 3
  # npts = 20

  # find the steep part and fit it
  crit_steep = .01
  tmp = which(abs(der$y) > crit_steep)
  tmp = tmp[2:length(tmp)] - tmp[1:(length(tmp)-1)]
  idx_steep = which(tmp == 1)
  ln1 = predict(lm(height ~ nclust, data = df2, subset = idx_steep), newdata = df2)


  idx2 = (nrow(df2) - npts - ignore_last):(nrow(df2) - ignore_last)
  ln2 = predict(lm(height ~ nclust, data = df2, subset = idx2), newdata = df2)
  dfh = data.frame(x = nclust, y = ln2)

  resid = df2$height - dfh$y

  # where does residual exceed more than x% of max?
  crit1 = 0.10
  crit2 = 0.05
  idx1 = min(which(resid < crit1 * max(resid)))
  nclust1 = df2$nclust[idx1]
  idx2 = min(which(resid < crit2 * max(resid)))
  nclust2 = df2$nclust[idx2]
  if (show) {
    plot(df2, ...)
    lines(dfh$x, dfh$y, col = 'blue')

    segments(x0 = df2$nclust[idx1], y0 = 0, x1 = df2$nclust[idx1], y1 = df2$height[idx1])
    segments(x0 = df2$nclust[idx2], y0 = 0, x1 = df2$nclust[idx2], y1 = df2$height[idx2])

    points(df2$nclust[idx1], df2$height[idx1], pch = 20, col = 'red', cex = 2)
    points(df2$nclust[idx2], df2$height[idx2], pch = 20, col = 'red', cex = 2)
    points(df2$nclust[idx1], df2$height[idx1], cex = 2)
    points(df2$nclust[idx2], df2$height[idx2], cex = 2)

    x = 10
    y = 0.75 * max(df2$height)
    text(x = x, y = y, labels = sprintf("Recommend between %d and %d clusters", nclust1, nclust2), pos = 4, cex = 1.5)
  }
  return(list(n1 = nclust1, n2 = nclust2))
}

# calculate the within-cluster sum of squares
# https://discuss.analyticsvidhya.com/t/what-is-within-cluster-sum-of-squares-by-cluster-in-k-means/2706
calculate_wss = function(clustering, max_clust = 100) {
  mat = clustering$data
  ndim = ncol(mat)              # number of dimensions

  wss = vector('numeric')
  for (k in 1:max_clust) {      # for all clusterings
    idx = cutree(as.hclust(clustering), k = k)
    d4 = 0
    for (i in 1:k) {       # for each cluster
      pts = which(idx == i)
      npts = length(pts)
      d3 = 0
      for (j in 1:ndim) {         # for each dimension of data
        xbar = mean(mat[pts, j])
        d2 = sum((mat[pts, j] - xbar) ^ 2)
        d3 = d3 + d2
      }
      d4 = d4 + d3
    }
    wss[k] = d4
  }

  wss
}

# estimate the first derivative using the symmetric difference quotient
lslope = function (kde, normalize=TRUE) {
  x = kde[,1]
  y = kde[,2]
  npts = length(y)
  yp = vector('numeric')
  for (i in 2:(npts - 1)) {
    yp[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
  }
  yp[1] = yp[2]
  yp[npts] = yp[npts - 1]
  yp[yp == -Inf] = 0.0

  if (normalize) {
    yp = yp / max (abs(yp))
  }
  res = list(x=x, y=yp)
  res
}

#construct a community object based on hierarchical (e.g. agnes) clustering
agnes_to_community = function(ag, nclust) {
  comm = list()
  vcount = nrow(ag$data)
  comm$vcount = vcount
  comm$names = as.character(1:vcount)

  membership = cutree(as.hclust(ag), k = nclust)
  comm$membership = membership
  class(comm) = "communties"

  comm
}

distributions_bins = function(fluster_obj, bin_indices) {
  parameters = fluster_obj$parameters
  fcs = fluster_obj$fcs   # flowFrame
  idx = bin_indices

  # collect all events in the bins in bin_indices
  data = exprs(fcs)
  fp = fluster_obj$fp
  idx = which(tags(fp)[[1]] %in% bin_indices)   # select all events in bins

  mn_vec = lo_vec = hi_vec = vector(mode = 'numeric')
  for (i in 1:length(parameters)) {
    mn_vec[i] = mean(data[idx ,parameters[i]])
    sdev = sd(data[idx,parameters[i]])
    lo_vec[i] = mn_vec[i] - 0.5 * sdev
    hi_vec[i] = mn_vec[i] + 0.5 * sdev
  }
  names(mn_vec) = names(lo_vec) = names(hi_vec) = parameters

  return(list(mn = mn_vec, lo = lo_vec, hi = hi_vec))
}


# this function calculates bin centers and standard deviations directly on
# events belonging to each cluster.
distributions_bins_all_clusters = function(fcs, fluster_obj) {
  parameters = fluster_obj$parameters

  if (is(fcs) == "flowSet") {
    fcs = as(fcs, "flowFrame")
  }
  fp = fluster_obj$fp
  tg = tags(fp)[[1]]
  nclust = length(fluster_obj$clustering$c_index)

  # set up matrices for results
  center = sdev = matrix(NA, nrow = nclust, ncol = length(parameters))
  colnames(center) = colnames(sdev) = parameters
  for (i in 1:nclust) {
    ev = which(tg %in% fluster_obj$clustering$c_index[[i]])
    for (p in parameters) {
      # calculate cluster center
      center[i, p] = mean(exprs(fcs)[ev, p])
      # calculate cluster sdev
      sdev[i, p] = sd(exprs(fcs)[ev, p])
    }
  }
  return(list(center = center, sdev = sdev))
}

categorical_phenotype_all_clusters = function(fluster_obj, sd_fac) {
  parameters = fluster_obj$parameters
  fcs = fluster_obj$fcs
  if (is(fcs) == "flowSet") {
    fcs = as(fcs, "flowFrame")
  }
    dbins = distributions_bins_all_clusters(fcs, fluster_obj)

    # label each parameter as either lo or hi
    nclust = length(fluster_obj$clustering$c_index)
    categ = list()
    for (i in 1:nclust) {
      p_category = rep("un", length = length(parameters))
      names(p_category) = parameters
      for (p in parameters) {
        thr = fluster_obj$modality$thresh[p]
        if (dbins$center[i, p] - 0.5 * sd_fac * dbins$sdev[i, p] > thr) {p_category[p] = 'hi'}
        if (dbins$center[i, p] + 0.5 * sd_fac * dbins$sdev[i, p] < thr) {p_category[p] = 'lo'}
      }
      p_category = factor(p_category, levels = c("lo", "un", "hi"))
      categ[[i]] = p_category
    }

  categ
}

# label a cluster categorically
# If the median for the cluster is above the parameter threshold, label it hi, otherwise lo
categorical_phenotype = function(fluster_obj, cluster = 1, sd_fac = 1.0) {
  parameters = fluster_obj$parameters
  idx = fluster_obj$clustering$c_index[[cluster]]
  res = distributions_bins(fluster_obj, bin_indices = idx)

  # label each parameter as either lo or hi
  p_category = rep("un", length = length(parameters))
  names(p_category) = parameters
  for (i in 1:length(parameters)) {
    # recover the standard deviation from the result of distributions_bins()
    sdev = 0.5 * (res$hi[i] - res$lo[i])
    if (res$mn[i] - sd_fac * sdev > fluster_obj$modality$thresh[i]) {p_category[i] = "hi"}
    if (res$mn[i] + sd_fac * sdev < fluster_obj$modality$thresh[i]) {p_category[i] = "lo"}
  }
  p_category = factor(p_category, levels = c("lo", "un", "hi"))

  p_category
}

compare_categories = function(c1, c2) {
  n_param = length(c1)
  if (length(which(c1 == c2)) == n_param) {
    res = TRUE
  } else {
    res = FALSE
  }
  res
}



# make a named color transparent (e.g. "dodgerblue2)
# alpha is in the interval [0, 1]
make_transparent = function(col, alpha = .5) {
  rgb = col2rgb(col)
  r = format(as.hexmode(rgb[1]), width = 2)
  g = format(as.hexmode(rgb[2]), width = 2)
  b = format(as.hexmode(rgb[3]), width = 2)
  a = format(as.hexmode(as.integer(alpha * 256)), width = 2)

  ret = paste0("#", r, g, b, a)

  ret
}

# calculate per-marker spreads and medians
spreads_and_meds = function(fluster_obj) {
  n_clust = max(fluster_obj$clustering$clst)
  ff = fluster_obj$fcs
  parameters = colnames(fluster_obj$centers)

  med = spread = dip = matrix(NA, nrow = n_clust, ncol = length(parameters))
  colnames(med) = colnames(spread) = colnames(dip) = parameters

  fp = fluster_obj$fp
  tg = tags(fp)[[1]]

  for (i in 1:n_clust) {
    n_min = 100    # BUGBUGBUG: think about this
    ev = which(tg %in% fluster_obj$clustering$c_index[[i]])
    for (j in 1:length(parameters)) {
      if (length(ev) >= n_min) {
        med[i, j] = median(exprs(ff)[ev, parameters[j]])
        spread[i, j] = IQR(exprs(ff)[ev, parameters[j]])
      } else {
        med[i, j] = 0
        spread[i, j] = Inf
      }
    }
  }
  return(list(med = med, spread = spread))
}

# select the dimmest and brightest tight clusters for a marker
select_extremes = function(med, spread) {
  tmp = med
  tmp = (tmp - min(tmp)) / (max(tmp) - min(tmp))
  tmp = 2 * (tmp - 0.5)
  tmp = tmp / spread

  idx = sort(tmp, index = TRUE)$ix
  len = length(idx)
  lo = idx[1]
  hi = idx[len]
  return(list(lo = lo, hi = hi))
}

# Use a spreadsheet to define functional phenotypes.
# Rows are functional phenotypes (for example, CD4_EM)
# columns are markers (e.g. CD4, CD8, CD3, ...)
# Entries are either 'lo', 'hi', or 'dc' (indicating don't care)
# First column contains the subset names
# colnames MUST match the marker labeling in the data set
# optional column labeled "Colors" will contain color coding of subsets
retrieve_categorical_definitions = function(file) {
  tab = read.csv(file, row.names = NULL, as.is = TRUE)

  # now turn the columns into properly ordered factors
  icol = which(tolower(colnames(tab)) == "color")
  idx = (1:ncol(tab))[-icol]
  idx = idx[-1]     # first col is fnames
  fnames = tab[ ,1]
  # for (i in idx) {
  #   tab[, i] = factor(tab[, i], levels = c("lo", "hi", "dc"))
  if (length(icol) == 1) {
    colors = tab[, icol]
  } else {
    colors = rep(NA, nrow(tab))
  }
  # }

  return(list(fnames = fnames, defs = tab[, idx], colors = colors))
}

assign_functional_names = function(fluster_obj, functional_definitions) {
  fd = as.matrix(functional_definitions$defs)
  fc = functional_definitions$colors
  fn = functional_definitions$fnames
  # add to the fluster object
  fluster_obj$labels$definitions = fd
  fluster_obj$labels$colors = fc
  fluster_obj$labels$names = fn

  n_clust = length(fluster_obj$clustering$c_index)
  func_phenotype = rep('unassigned', length = n_clust)
  func_color = rep("black", length = n_clust)
  pheno = fluster_obj$clustering$phenotype
  for (i in 1:n_clust) {
    for (j in 1:nrow(fd)) {
      idx = which(fd[j, ] != "dc")
      cmarkers = colnames(fd)[idx]
      def = fd[j, idx]
      # extract values for cluster
      ctype = as.character(pheno[[i]][cmarkers])
      if (length(which(ctype == def)) == length(idx)) {
        func_phenotype[i] = fn[j]
        func_color[i] = fc[j]
      }
    }
  }

  # tack onto fluster_obj
  fluster_obj$clustering$func_phenotype = func_phenotype
  fluster_obj$clustering$func_color = func_color

  fluster_obj
}

# calculate skewness of a distribution
dist.moment = function(x, r = 3) {
  # calculate 3rd moment
  n = length(x)
  mn = mean(x)
  sm = sum((x - mn) ^ r)
  mom = sm / n

  standardized.moment = mom / (sd(x) ^ r)

  standardized.moment
}



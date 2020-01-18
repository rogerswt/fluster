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
#  method = "median" returns median +- quartiles.
#  method = "mean" returns mean +- standard deviation
calculate_bin_phenotypes = function(fp, fs, method=c("median", "mean")) {
  parameters = parameters(fp)
  n_bins = 2 ^ nRecursions(fp)
  n_parameters = length(parameters)
  center = matrix(NA, nrow = n_parameters, ncol = n_bins)
  rownames(center) = parameters
  range = matrix(NA, nrow = n_parameters, ncol = n_bins)
  rownames(range) = parameters

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
      if (method == "median") {
        center[j, i] = median(vals, na.rm = TRUE)
        range[j, i] = (quantile(x = vals, probs = 0.75, na.rm = TRUE) -
                         quantile(x = vals, probs = 0.25, na.rm = TRUE)) / 2
      } else {

      }
    }
  }
  return(list(center = center, range = range))
}

# assumes biexp vert scale
draw_color_scale = function(min_col_value = 0, max_col_value = 5, ...) {
  ll = -0.5
  ul = bx(262143)

  vec = seq(ll, ul, length.out = 500)
  cols = pcolor(pvalue = vec, min_value = min_col_value, max_value = max_col_value)

  opar = par(mar = c(0, 5, 0, 0) + .1)
  plot(0, 0, pch = '', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "Fluorescence Intensity",
       xlim = c(0, 5), ylim = c(ll, ul), ...
  )
  for (i in 1:length(vec)) {
    y = vec[i]
    segments(x0 = 0, y0 = y, x1 = 1, y1 = y, col = cols[i], lwd = 3)
  }
  ax(axis = 2, type = 'biexp', ...)
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
      if(res$strength != 0) {
        # add an edge
        gcomm = add_edges(gcomm, c(i, j))
        E(gcomm)$strength[k] = res$strength
        k = k + 1
        cat(i, j, "\n")
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



plot_community_graph = function(g, marker, vs = 3, ms = 0, log.size = TRUE, vertex.frame = TRUE, cex.main) {
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
  plot(g, vertex.size = vsize, vertex.color = pcolor(vertex_attr(g, marker)),
       vertex.frame.color = vfc, vertex.label = NA)
  title(main = marker, cex.main = cex.main)
}


plot_comm_spread = function(g, markers = NULL, vs = 3, ms = 1, log.size = TRUE, vertex.frame = FALSE, cex.main) {
  if (is.null(markers)) {
    markers = vertex_attr_names(g)
    markers = markers[-which(markers == "size")]
  }

  # calculate plot layout
  n = length(markers) + 1
  sq = sqrt(n)
  frac = sq - floor(sq)
  if(frac == 0) {
    ac = dn = floor(sq)
  } else {
    ac = floor(sq) + 1
    dn = ceiling(n / ac)
  }

  par(mfrow = c(dn, ac), mar = c(0, 0, 2, 0))
  for (i in 1:length(markers)) {
    marker = markers[i]
    plot_community_graph(g, marker, vs = vs, ms = ms, log.size = log.size, vertex.frame = vertex.frame, cex.main)
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
  ignore_last = 3
  npts = 20

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




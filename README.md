# fluster
R package, Fingerprint-based clustering for flow (or mass) cytometry

## Installing the package
*Please note that the vignette takes a few minutes
to build because the example dataset is a bit on the large side (42 mb).*
```
devtools::install_github("rogerswt/fluster", build_vignettes = TRUE)
```

*Also note that some platforms are having trouble building the vignette.  If this is
your situation, you can omit the vignette.*  
```
devtools::install_github("rogerswt/fluster")
```
## Background

The prevalence of high-dimensional flow (or mass) cytometry data is increasing at a rapid pace.
The difficulties and limitations of manual analysis of such data are obvious to anyone who has
tried to do that.  The alternative is to find or create automated analysis approaches that
facilitate this process.

The *fluster* package is closely related to another package, *panoplyCF*.  Both packages 
rely on the speed and efficiency of reducing the number of items to cluster using 
Cytometric Fingerprinting (via package *flowFP*).  The difference between the two is that,
in the case of panoplyCF, clustering is done in a 2-dimensional space after 
manifold learning dimensionality reduction using the t-SNE algorithm.  Fluster on the
other hand performs hierarchical agglomerative clustering directly on the bin centroids.
For visualization fluster provides a minimal spanning tree (MST) graph representation of
the multivariate clusters to aid in interpretation.

## Fluster
Fluster first computes high-resolution 
Cytometric Fingerprinting (CF) bins using FlowFP.  It then computes the centroids of the bins
in high-dimensional space by taking the medians for each independent variable for all
events contained in each bin. The final step is to perform conventional hierarchical agglomerative clustering of the
bins, using their high-dimensional centroids.  The goal is to create clusters
that are homogeneous in high-dimensional space as determined by the tightness of 
the distribution of the independent variables within each cluster.

Data can be easily mapped to the clusters via the tagging that flowFP does on instances,
carried through the bin indices in the fluster object.

## Citation
Fluster is as yet unpublished.  Please acknowledge me if you use fluster for publication,
and drop me a note as well!

[Wade Rogers](mailto:wade.rogers@spcytomics.com)

[Still Pond Cytomics](https://spcytomics.com)

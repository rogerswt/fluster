#
# Transform values in a flowFrame object for "biexponential" plotting
#

############################################################
#	The following functions are utilities
#	and should NOT be exposed
############################################################

kernel_function <- function (x, a=0.002) {
    return (log(0.5*(x + sqrt((x)^2 +1/a^2)))/log(10))
}

biexp.transform  <- function (x, a=0.002, full_scale=262143, jitter=TRUE) {

  if (!is.numeric(x)) {
    stop ("x must be numeric")
  }
    # add a small amount of random noise to smooth things out at low vals
	if (jitter) {
		noise <- rnorm (length(x), sd=1)
		x <- x + noise
	}

	tmp <- kernel_function (x, a)

	# offset downwards, then compensate by multiplying asymptotic value
	offset <- kernel_function(0, a)
	intercept = log10(full_scale)
	fac <- intercept / (kernel_function (full_scale, a) - offset)

	tmp <-  fac * (tmp - offset)
	return (tmp)
}


inv_kernel_function <- function (x, a=0.002) {
	out <- .25 * 10^(-x) * (-1/a^2 + 4 * 10^(2*x))
	out
}

inv.biexp.transform <- function (x, a=.002, full_scale=262143) {

	# reverse the maneuvers
	offset <- kernel_function(0, a)
	intercept = log10(full_scale)
	fac <- intercept / (kernel_function (full_scale, a) - offset)

	tmp <- inv_kernel_function(x/fac + offset, a)
	tmp
}

############################################################
#	This function should be exposed
############################################################
biexp <- function (f, a=0.002, params) {

	# if using symbolic names for params, convert to numeric
	if (is.character(params)) {
		params <- which(colnames(f) %in% params)
	}

   # make a copy of the input object
   fout <- f

   # compute transformed object
   for (parm in params) {
      exprs(fout)[,parm] <- biexp.transform (exprs(f)[,parm], a)
   }
   
   # set the minRange and maxRange values appropriately
   parameters(fout)$minRange[params] <- biexp.transform(parameters(fout)$minRange[params])
   parameters(fout)$maxRange[params] <- biexp.transform(parameters(fout)$maxRange[params])

   return (fout)
}


############################################################
#	This function should be exposed
############################################################
inv.biexp <- function (f, a=0.002, params) {
	# if using symbolic names for params, convert to numeric
	if (is.character(params)) {
		params <- which(colnames(f) %in% params)
	}

   # make a copy of the input object
   fout <- f

   # compute transformed object
   for (parm in params) {
      exprs(fout)[,parm] <- inv.biexp.transform (exprs(f)[,parm], a)
   }
   
   # set the minRange and maxRange values appropriately
   parameters(fout)$minRange[params] <- inv.biexp.transform(parameters(fout)$minRange[params])
   parameters(fout)$maxRange[params] <- inv.biexp.transform(parameters(fout)$maxRange[params])

   return (fout)
}


############################################################
#	This function should be exposed - flowCore compatibility
############################################################
biexpTransform <- function (transformId="mybiexp", a=0.002, full_scale=262143, jitter=TRUE) {
	t <- new ("transform", .Data = function (x) biexp.transform(x, a, full_scale, jitter))
	t@transformationId = transformId
	t
}

# shorthand form
bx = function(x) {
  idx = which(x == -Inf)
  res = biexp.transform(x, jitter = FALSE)
  if(length(idx) > 0) {
    res[idx] = -Inf
  }
  res
}

ibx = function(x) {
  res = inv.biexp.transform(x)
  
  res
}


######################################################
# Generates axes used in pplot
#
#
#######################################################
ax <- function (axis=1, instrument, type, ticksize = 2,...) {
  
  # instrument = match.arg(instrument)
  # type       = match.arg(type)
  
  if (instrument == "diva") {
    channels=262143
    decades=8
  } else {
    channels=10000
    decades=4
  }
  all.ticks<-NULL
  
  if (!(type %in% c("biexp", "log", "linear"))) {
    cat ("invalid axis type\n")
    return()
  }
  if (type == "log") {
    start_decade=0
    major<-(10^(start_decade:decades))
    for(i in 1:(length(major)-1)){
      all.ticks<-c(all.ticks,seq(major[i],major[i+1],l=10))
    }
    all.ticks <- log(all.ticks, base=10); #Log base 10
    major     <- log(major, base=10);
  }
  if (type == "biexp") {
    start_decade=2
    major<-(10^(start_decade:decades))
    for(i in 1:(length(major)-1)){
      all.ticks<-c(all.ticks,seq(major[i],major[i+1],l=10))
    }
    all.ticks <- biexp.transform(all.ticks);
    major     <- biexp.transform(major);
  }
  if (type == "linear") {
    axis (side=axis)
    return()
  }
  axis(side=axis, all.ticks, label=FALSE, tcl=-0.25)
  axis(side=axis, at=major, lwd.ticks=ticksize, label=make.labels(start_decade, decades),...)
}

make.labels <- function (first_decade, last_decade) {
  label <- vector('character')
  n_labels <- last_decade - first_decade + 1
  for (i in 1:n_labels) {
    label[i] <- paste ("10^", first_decade + i - 1, sep="")
  }
  ex = parse(text=label)
  
  return(ex)
}

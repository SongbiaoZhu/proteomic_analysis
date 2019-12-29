checkCPUs <- function(cpus) {
  ## Determine appropriate number of CPUs for parallelization
  maxCores <- detectCores()
  if (is.numeric(cpus) && (cpus > maxCores)) {
    warning("Selected number of cores (", cpus, ") exceeds those available on your device (", maxCores, "). Using ", maxCores, " cores.")
    cpus <- maxCores
  } else if (identical(cpus, "max")) {
    cpus <- maxCores
  } else if (!is.numeric(cpus)){
    stop(paste("Invalid argument",cpus, "for 'cpus'."))
  }
  return(cpus)
}
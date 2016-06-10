#'  Cut a "slice" of a netCDF variable (usually 2D)
#'
#'  @param nc netCDF handle, from \code{\link[ncdf4]{nc_open}}
#'  @param varname variable name as a character string
#'  @param ... indexes to extract for each dimension; either named arguments, in which case the unspecified dimensions are fully extracted; or unnamed arguments, in which case dimensions are assumed to be in the same order as in the netCDF variable and NA elements are fully extracted
#'
#' @export
#'
#' @importFrom ncdf4 ncvar_get
nc_slice <- function(nc, varname, ...) {

  # check arguments
  if ( ! "ncdf4" %in% class(nc) ) {
    stop("nc needs to be a ncdf4 handle, created with nc_open() from package ncdf4")
  }
  varname <- match.arg(varname, names(nc$var))

  # get dimensions of this variable
  dims <- nc$var[[varname]]$dim
  dimNames <- sapply(dims, function(x) {x$name})
  names(dims) <- dimNames   # makes it easier to extract them

  # detect whether data is extracted based on dimension *names* or positions and construct the extraction list appropriately
  # the extraction list is a list of indexes to extract, one element for each dimension of the chosen variable
  extract <- list(...)
  extractNames <- names(extract)

  if ( ! is.null(extractNames) ) {
    # Use the names

    # check that all names are actually dimension names
    if ( ! all(extractNames %in% dimNames) ) {
      notOK <- extractNames[! extractNames %in% dimNames ]
      if (length(notOK) > 1) {
        notOK <- paste(notOK, collapse=", ")
        mess <- " are not dimensions of "
      } else {
        mess <- " is not a dimension of "
      }
      stop(notOK, mess, varname)
    }

    # detect matching names
    matching <- match(dimNames, extractNames)

    # create the extraction list with NA when the dimension is not specified
    extract <- extract[matching]
    extract[is.na(matching)] <- NA  # by default these are NULL

  } else {
    # Consider that the list contains all the dimensions in the correct order (or is empty entirely)
    # check that these conditions are satisfied
    
    if ( length(extract) == 0 ) {

      # if the list is empty, just put NAs everywhere: we extract everything
      extract <- as.list(rep(NA, times=length(dims)))
                  
    } else if ( length(extract) != length(dims) ) {
      
      # the extraction arguments provided are not appropriate
      stop("You specified ", length(extract)," dimensions but variable", varname, " has ", length(dims), " dimensions:", paste(dimNames, collapse=", "), ".")

    }
  }
  
  # make sure the extraction list has appropriate names
  names(extract) <- dimNames
  # TODO check wether that is useful

  # replace NA dimensions with all indexes for the corresponding dimensions
  for ( i in 1:length(dims) ) {
    if ( all(is.na(extract[[i]])) ) {
      extract[[i]] <- 1:(dims[[i]]$len)
    }
  }
  
  # compute start and count vectors for ncvar_get
  start <- sapply(extract,min)
  stop  <- sapply(extract,max)
  count <- stop - start+1
  
  # extract slice
  slice <- ncvar_get(nc, varname, start=start, count=count)
  
  # detect which dimensions have an actual dimension in the extracted array i.e. those with a count > 1
  extractedDimNames <- names(count)[which(count > 1)]

  # get the values of the coordinates extracted in these dimensions
  coords <- lapply(extractedDimNames, function(x) {
    dims[[x]]$vals[extract[[x]]]
  })
  names(coords) <- extractedDimNames
  
  # give the slice dimension names, which makes it easier to melt/plot
  dimnames(slice) <- coords
  
  # give it a specific class to make clever stuff with it afterwards
  class(slice) <- c("nc_slice", class(slice))
  
  # add the variable name as attribute
  attr(slice, "varname") <- varname
  
  return(slice)
}

# as.xyz.nc_slice <- function(x) {
#   coords <- dimnames(x)
# }

#'  Plot a "slice" extracted from a NetCDF file
#'
#'  @param x a NetCDF slice object, extracted via \code{\link{nc_slice}}
#'  @param ... passed to \code{image}
#'
#' @export
#' @seealso \code{\link{nc_slice}}
plot.nc_slice <- function(x, ...) {
  
  if ( length(dim(x)) != 2 ) {
    stop("Can only plot a slice with 2 dimensions")
  }
  
  # make a list for image/persp/...
  out <- dimnames(x)
  out <- lapply(out, as.numeric)
  names(out) <- c("x", "y")
  out$z <- x

  # plot an image
  image(out, ...)

  return(invisible(out))
}

#'  Convert a "slice" extracted from a NetCDF file into a data.frame
#'
#'  @param x a NetCDF slice object, extracted via \code{\link{nc_slice}}
#'  @param ... passed to generic method
#'
#' @export
#' @seealso \code{\link{nc_slice}}
#' @importFrom reshape2 melt
as.data.frame.nc_slice <- function(x, ...) {
  out <- melt(x, value.name=attr(x, "varname"))
  out <- as.data.frame(out, ...)
  return(out)
}


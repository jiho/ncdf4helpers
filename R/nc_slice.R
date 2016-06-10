#' Read a slice from a netCDF file
#'
#' Reads an array of data (a slice) from an existing netCDF file. Typically this array is two-dimensional and meant to be plotted.
#'
#' @param nc netCDF handle, from \code{\link[ncdf4]{nc_open}}
#' @param varid the variable to read the data from. Needs to be an actual variable, not a dimension variable. Can be a string with the name of the variable or an object of class \code{ncvar4}. If left unspecified, the function will determine if there is only one variable in the file and, if so, read from that.
#' @param ... ranges of the values to extract, for each dimension. Each range is a vector of numbers, in the unit of the corresponding dimension, and everything between the minimum and maximum of these numbers will be extracted. Ranges can be specified either as named arguments, whose name matches one of the dimensions, or are assumed to be in the same order as the dimensions of the netCDF variable. Unspecified or NA elements are fully extracted.
#'
#' @return A list containing vectors for each dimension (of length > 1) and one array containing the extracted data. For two-dimensional slices, this list is suitable for plotting functions such as \code{\link[graphics]{image}}, \code{\link[graphics]{contour}} and \code{\link[graphics]{persp}} after changing the names to "x", "y", and "z". Convenience functions are provided to make this conversion, convert the slice into a data.frame, or plot the slice directly.
#'
#' @seealso \code{\link{as.xyz.ncslice4}} and \code{\link{as.data.frame.ncslice4}} to convert a slice to various usual formats; \code{\link{plot.ncslice4}} to plot a slice.
#' 
#' @export
#'
#' @examples
#' nc <- nc_open(system.file("extdata", "sresa1b_ncar_ccsm3.nc", package="ncdf4helpers"))
#' # extract everything
#' x <- ncvar_slice(nc, "tas")
#' # inspect the result
#' str(x)
#' # extract some specific ranges
#' x <- ncvar_slice(nc, "tas", lon=c(100, 200))
#' x <- ncvar_slice(nc, "tas", lat=10:80)
#' x <- ncvar_slice(nc, "tas", lon=c(0,300), lat=c(10:20))
#' x <- ncvar_slice(nc, "tas", lon=c(0,50), lat=NA)
#' nc_close(nc)
#'
#' # plot the result
#' image(as.xyz(x))
#' # or simply
#' image(x)
#'
#' # convert the slice into a data.frame for later manipulation
#' head(as.data.frame(x))
ncvar_slice <- function(nc, varid, ...) {
  # check arguments
  # NB: this checks for the class of nc
  #     this uses an internal, not exported, function in ncdf4 and may therefore be fragile
  var <- ncdf4:::vobjtovarid4(nc, varid)
  var_name <- names(nc$var)[var$list_index]
  
  # prepare the extraction list
  # = a list with as many elements as the variables has dimensions, in the order of the variables dimensions and containing min and max of the values to extract for each dimension

  # initialise the extraction list
  extract <- list(...)
  
  # get dimensions of this variable
  dims <- nc$var[[var$list_index]]$dim
  dim_names <- sapply(dims, function(x) {x$name})
  names(dims) <- dim_names   # makes it easier to extract the values later

  # check the number of extracted dimensions
  if (length(extract) > length(dims)) {
    stop("You specified ", length(extract)," ranges to extract but `", var_name, "` has only ", length(dims), " dimensions: ", paste(dim_names, collapse=", "), ".")
  }

  # fill missing names assuming the variables are in order
  if (is.null(names(extract))) {
    # all missing names
    names(extract) <- dim_names[seq(along=extract)]
    # NB: using seq(along=...) here allows to deal with empty ...
  } else {
    # some missing names
    missing_names <- names(extract) == ""
    names(extract)[missing_names] <- dim_names[missing_names]
  }
  
  # check that provided names are OK
  invalid_names <- setdiff(names(extract), dim_names)
  if (length(invalid_names) != 0) {
    stop(paste(invalid_names, collapse=", "), " not among valid dimension names. Should be one of: ", paste(dim_names, collapse=", "))
  }
  
  # reorder dimensions and add missing dimensions
  extract <- extract[dim_names]
  names(extract) <- dim_names
  
  # compute the extracted range for each dimension
  # this also replaces NA or NULL dimensions with a full range ]-Inf,+Inf[
  extract <- lapply(extract, function(x) {
    suppressWarnings(sort(range(x, na.rm=T)))
  })
  
  # for each dimension, get the coordinates of the extracted points and compute start and count for ncvar_get
  for (dim in dim_names) {
    # get the range and values of the dimensions
    r <- extract[[dim]]
    vals <- nc$dim[[dim]]$vals
    # compute the indexes inside the range
    idx <- which(vals >= r[1] & vals <= r[2])
    # extract the coordinates, start, stop and compute count
    vals <- vals[idx]
    start <- min(idx)
    stop  <- max(idx)
    count <- stop - start + 1
    # store result as a list
    extract[[dim]] <- list(vals=vals, start=start, count=count)
  }

  # prepare the arguments for ncvar_get
  start <- sapply(extract, `[[`, "start")
  count <- sapply(extract, `[[`, "count")

  # extract the data
  dat <- ncdf4::ncvar_get(nc, varid, start=start, count=count)
  
  # store coordinates and data in a list (such as the ones used for image, persp, etc.)
  coords <- lapply(extract, `[[`, "vals")
  slice <- coords[count > 1] # eliminate coordinates with one extracted element
  slice[[var_name]] <- dat
  
  # give it a specific class to make clever stuff with it afterwards
  class(slice) <- c("ncslice4", class(slice))
  
  return(slice)
}


#' Coerce to a x,y,z list
#' 
#' Generic function to convert an object into an x,y,z list suitable for \code{\link[graphics]{image}}, \code{\link[graphics]{contour}}, \code{\link[graphics]{persp}}, etc.
#'
#' @param x an R object
#'
#' @export
as.xyz <- function(x) UseMethod("as.xyz", x)

#' Convert a NetCDF slice into xyz list
#'
#' Convert a slice read from a NetCDF file with \code{\link{ncvar_slice}} into a x,y,z list suitable for \code{\link[graphics]{image}}, \code{\link[graphics]{contour}}, \code{\link[graphics]{persp}}, etc.
#'
#' @param x a NetCDF slice object, extracted via \code{\link{ncvar_slice}}
#'
#' @seealso \code{\link{ncvar_slice}} to read the slice and \code{\link{as.data.frame.ncslice4}} to convert it into another format.
#'
#' @export
#'
#' @examples
#' nc <- nc_open(system.file("extdata", "sresa1b_ncar_ccsm3.nc", package="ncdf4helpers"))
#' x <- ncvar_slice(nc, "tas")
#' nc_close(nc)
#' # inspect and transform the result
#' str(x)
#' x <- as.xyz(x)
#' str(x)
as.xyz.ncslice4 <- function(x) {
  if (length(x) > 3) {
    stop("Can only handle a 2-dimensional slice. This one has ", length(x)-1)
  }
  # change the names for persp and the like
  names(x) <- c("x", "y", "z")
  # remove the class to avoid infinite recursion
  x <- unclass(x) 
  return(x)
}


#' Convert a NetCDF slice into a data.frame
#'
#' Convert a slice read from a NetCDF file with \code{\link{ncvar_slice}} into a data.frame for further manipulation of the data.
#'
#' @param x a NetCDF slice object, extracted via \code{\link{ncvar_slice}}.
#' @param row.names \code{NULL} or a character vector giving the row names for the data frame.  Missing values are not allowed.
#' @param optional logical. If \code{TRUE}, setting row names and converting column names (to syntactic names: see \code{\link[base]{make.names}}) is optional.
#' @param ... ignored.
#'
#' @seealso \code{\link{ncvar_slice}} to read the slice and \code{\link{as.xyz.ncslice4}} to convert it into another format.
#'
#' @export
#' @import reshape2
#'
#' @examples
#' nc <- nc_open(system.file("extdata", "sresa1b_ncar_ccsm3.nc", package="ncdf4helpers"))
#' x <- ncvar_slice(nc, "tas")
#' nc_close(nc)
#' head(as.data.frame(x))
as.data.frame.ncslice4 <- function(x, row.names=NULL, optional=FALSE, ...) {
  n <- length(x)
  # convert the data
  d <- reshape2::melt(x[[n]])
  # replace indexes by actual values for each dimension
  for (i in 1:(n-1)) {
    d[,i] <- x[[i]][d[,i]]
  }
  names(d) <- names(x)
  # handle as.data.frame generic arguments
  d <- as.data.frame(d, row.names=row.names, optional=optional)
  return(d)
}


#' Plot a slice extracted from a NetCDF file
#'
#' @param x a NetCDF slice object, extracted via \code{\link{ncvar_slice}}
#' @param ... passed to the plotting function
#'
#' @seealso \code{\link{ncvar_slice}}
#'
#' @examples
#' nc <- nc_open(system.file("extdata", "sresa1b_ncar_ccsm3.nc", package="ncdf4helpers"))
#' x <- ncvar_slice(nc, "tas", lon=10:30, lat=30:50)
#' nc_close(nc)
#' image(x)
#' contour(x, add=TRUE)
#' persp(x, theta=210, phi=30)
#'
#' @export
#' @import graphics
#' @rdname plot.ncslice4
image.ncslice4 <- function(x, ...) { graphics::image(as.xyz(x), ...) }
#' @export
#' @rdname plot.ncslice4
plot.ncslice4 <- image.ncslice4
#' @export
#' @rdname plot.ncslice4
contour.ncslice4 <- function(x, ...) { graphics::contour(as.xyz(x), ...) }
#' @export
#' @rdname plot.ncslice4
persp.ncslice4 <- function(x, ...) { graphics::persp(as.xyz(x), ...) }

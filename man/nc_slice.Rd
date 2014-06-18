\name{nc_slice}
\alias{nc_slice}
\title{Cut a "slice" of a netCDF variable (usually 2D)}
\usage{
nc_slice(nc, varname, ...)
}
\arguments{
  \item{nc}{netCDF handle, from
  \code{\link[ncdf4]{nc_open}}}

  \item{varname}{variable name as a character string}

  \item{...}{indexes to extract for each dimension; either
  named arguments, in which case the unspecified dimensions
  are fully extracted; or unnamed arguments, in which case
  dimensions are assumed to be in the same order as in the
  netCDF variable and NA elements are fully extracted}
}
\description{
Cut a "slice" of a netCDF variable (usually 2D)
}


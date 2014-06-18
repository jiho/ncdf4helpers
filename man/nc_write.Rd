\name{nc_write}
\alias{nc_write}
\title{Write a data.frame containing data on a grid as a netCDF file}
\usage{
nc_write(d, file, dimensions, variables = NULL, units = NULL,
  attributes = NULL, global.attributes = NULL, missval = -99999)
}
\arguments{
  \item{d}{data.frame containing coordinates (along
  dimensions) and variables}

  \item{file}{output file path}

  \item{dimensions}{names or indexes of columns in d
  representing the data dimensions}

  \item{variables}{names or indexes of columns in d
  represented the measured variables. If NULL, all columns
  which are not dimensions are considered as variables}

  \item{units}{named vector or list containing character
  strings which specify the units of dimensions and
  variables. The names must be those of columns in d.
  Defaults to no unit}

  \item{attributes}{list of lists containing attributes
  (second-level lists) for each dimension and variable
  (first level list). The list must be named with names
  which are dimensions or variables (i.e. columns of d).
  Each second-level list is also a named list which stores
  attributes as key/value pairs, where the key is the name
  of the element of the list and the value is the content
  of that element}

  \item{global.attributes}{list of key/value pairs similar
  to the variables/dimensions attributes above but stored
  as global attributes of the file}

  \item{missval}{value used to code missing values in the
  netCDF file}
}
\description{
Write a data.frame containing data on a grid as a netCDF
file
}
\examples{
# create data, on a grid, stored in a data.frame
x <- 1:3
y <- 1:2
d <- expand.grid(x=x, y=y)
d$V1 <- runif(nrow(d))
d$V2 <- rnorm(nrow(d), mean=90)

file_name <- tempfile()
nc_write(d, file=file_name,
  # dimensions are specified, variables are not
  # => variables are all other columns (V1 and V2)
  dimensions = c("x", "y"),
  # units are a named list
  units = c(x="m", y="m", V1="degrees Celsius", V2="volts"),
  # variables attributes are specified as a list of lists
  attributes = list(
    x = list(long_name="x coordinate", standard_name="x"),
    y = list(long_name="y coordinate", standard_name="y"),
    V1 = list(long_name="Temperature at surface", standard_name="temp"),
    V2 = list(long_name="Tension between poles", standard_name="tension")
  ),
  # global attributes are not associated with any particular variable
  global.attributes = list(
    author = "Pr Foobar",
    convention = "CF-1.2",
    date_created = "2014-06-19"
  )
)
}


#' Write a data.frame containing data on a grid as a netCDF file
#'
#' @param d data.frame containing coordinates (along dimensions) and variables
#'
#' @param file output file path
#'
#' @param dimensions names or indexes of columns in d representing the data dimensions
#'
#' @param variables names or indexes of columns in d represented the measured variables. If NULL, all columns which are not dimensions are considered as variables
#' 
#' @param units named vector or list containing character strings which specify the units of dimensions and variables. The names must be those of columns in d. Defaults to no unit
#'
#' @param attributes list of lists containing attributes (second-level lists) for each dimension and variable (first level list). The list must be named with names which are dimensions or variables (i.e. columns of d). Each second-level list is also a named list which stores attributes as key/value pairs, where the key is the name of the element of the list and the value is the content of that element
#'
#' @param global.attributes list of key/value pairs similar to the variables/dimensions attributes above but stored as global attributes of the file
#'
#' @param missval value used to code missing values in the netCDF file
#'
#' @export
#'
#' @importFrom ncdf4 ncdim_def ncvar_def nc_create ncvar_put ncatt_put nc_close
#' @importFrom plyr llply l_ply
#' @importFrom reshape2 acast
#'
#' @examples
#' # create data, on a grid, stored in a data.frame
#' x <- 1:3
#' y <- 1:2
#' d <- expand.grid(x=x, y=y)
#' d$V1 <- runif(nrow(d))
#' d$V2 <- rnorm(nrow(d), mean=90)
#' 
#' file_name <- tempfile()
#' nc_write(d, file=file_name,
#'   # dimensions are specified, variables are not
#'   # => variables are all other columns (V1 and V2)
#'   dimensions = c("x", "y"),
#'   # units are a named list
#'   units = c(x="m", y="m", V1="degrees Celsius", V2="volts"),
#'   # variables attributes are specified as a list of lists
#'   attributes = list(
#'     x = list(long_name="x coordinate", standard_name="x"),  
#'     y = list(long_name="y coordinate", standard_name="y"),  
#'     V1 = list(long_name="Temperature at surface", standard_name="temp"),  
#'     V2 = list(long_name="Tension between poles", standard_name="tension") 
#'   ),
#'   # global attributes are not associated with any particular variable
#'   global.attributes = list(
#'     author = "Pr Foobar",
#'     convention = "CF-1.2",
#'     date_created = "2014-06-19"
#'   )
#' )
#'
nc_write <- function(d, file, dimensions, variables=NULL, units=NULL, attributes=NULL, global.attributes=NULL, missval=-99999) {

   # check that dimensions are in the data
   if ( is.character(dimensions) & ! all(dimensions %in% names(d)) ) {
      stop("All dimensions must be columns of d")
   }
   if (is.numeric(dimensions)) {
      dimensions <- names(d)[dimensions]
   }

   # check that all variables are in the data
   if ( is.null(variables) ) {
      variables <- setdiff(names(d), dimensions)
   } else {
      if ( is.character(variables) & ! all(variables %in% names(d)) ) {
         stop("All variables must be columns of d")
      }
      if (is.numeric(variables)) {
         variables <- names(d)[variables]
      }
   }

   # prepare the units list
   # defaults
   allVars <- c(dimensions, variables)
   u <- as.list(rep("", times=length(allVars)))
   names(u) <- allVars
   # include those specified as arguments
   units <- as.list(as.character(units))
   if ( ! all(names(units) %in% names(d)) ) {
      stop("All units must relate to columns of d")
   }
   u[names(units)] <- units
   # NB: work even when units is NULL
   units <- u
   
   # prepare attributes list
   a <- as.list(rep(NA, times=length(allVars)))
   names(a) <- allVars
   # include those specified as arguments
   if ( ! all(names(attributes) %in% names(d)) ) {
      stop("All attributes must relate to columns of d")
   }
   a[names(attributes)] <- attributes
   attributes <- a

   # define dimensions
   dims <- llply(dimensions, function(dim, units) {
       x <- d[,dim]
       vals <- sort(unique(x))
       ncdim_def(name=dim, units=units[[dim]], vals=vals)
   }, units=units)

   # define variables
   vars <- llply(variables, function(var, units, missval) {
      type <- class(d[,var])
      precision <- switch(type,
         numeric = "double",
         factor = "char",
         integer = "integer",
         "double"
      )
      ncvar_def(name=var, units=units[[var]], dim=dims, missval=missval, prec=precision)
   }, units=units, missval=missval)

   # create the file based on the variable definition
   # NB: automatically creates the dimensions and the dimensions variables
   nc <- nc_create(filename=file, vars=vars)

   # reformat variables as matrices and store them in the file
   l_ply(variables, function(var, missval) {
      x <- acast(d, formula=as.list(dimensions), value.var=var, fill=missval)
      ncvar_put(nc, varid=var, vals=x)
   }, missval=missval)

   # add attributes
   for (varName in names(attributes)) {
      if ( is.list(attributes[[varName]]) ) {
         for (attName in names(attributes[[varName]])) {
            ncatt_put(nc, varName, attname=attName, attval=attributes[[varName]][[attName]])
         }
      }
   }

   # add global attributes
   if ( ! is.null(global.attributes)) {
      for (attName in names(global.attributes)) {
         ncatt_put(nc, 0, attname=attName, attval=global.attributes[[attName]])
      }
   }

   # close the file
   nc_close(nc)

   return(invisible(file))
}


#' Write a data.frame containing data on a geographical grid as a netCDF file
#'
#' @inheritParams nc_write
#'
#' @param crs Coordinate Reference System specification. Currently only WGS84 is supported.
#'
#' @param ... passed to nc_write
#'
#' @export
#'
#' @importFrom ncdf4 nc_open nc_close ncatt_put ncvar_def ncvar_add
#'
#' @examples
#' # load altitude data over the Mediterranean Sea
#' data(med)
#' file_name <- tempfile()
#' nc_write_map(med, file=file_name,
#'   # dimensions are lon and lat by default
#'   # variables are all other columns (altitude)
#'   # units of lon and lat are set automatically
#'   units = c(altitude="m"),
#'   # as are attributes
#'   attributes = list(
#'     altitude = list(long_name="Altitude relative to mean sea level", standard_name="altitude") 
#'   ),
#'   # some global attributes are set automatically. Those will be added
#'   global.attributes = list(
#'     author = "Pr Foobar",
#'     date_created = "2014-06-19"
#'   )
#' )
#'
nc_write_map <- function(d, file, dimensions=c("lon", "lat"), units=NULL, attributes=NULL, crs="WGS84", ...) {
   
   # check arguments
   crs <- match.arg(tolower(crs), choices="wgs84")
   
   # add lon and lat characteristics to the units and attributes arguments
   units <- c(units, list(lon="degrees_east", lat="degrees_north"))
   
   attributes <- c(attributes, list(
      lon=list(long_name="Longitude coordinate", standard_name="longitude", axis="X"),
      lat=list(long_name="Latitude coordinate", standard_name="latitude", axis="Y")
   ))
   
   # create the file
   file <- nc_write(d=d, file=file, dimensions=dimensions, units=units, attributes=attributes, ...)
   
   # add the definition of a CRS
   nc <- nc_open(file, write=TRUE)
   
   # specify the mapping for all variables
   vars <- names(nc$var)   
   for (var in vars) {
      ncatt_put(nc, var, attname="grid_mapping", attval="crs")
   }
   
   # define the Coordinate Reference System variable (empty with only attributes)
   crsVar <- ncvar_def(name="crs", units="", dim=list(), prec="integer")
   ncvar_add(nc, crsVar)
   nc_close(nc)   # actually add the variable to the file
                  # NB: nc_sync is not enough apparently

   # actually define the CRS through its attributes
   nc <- nc_open(file, write=TRUE)
   if (crs == "wgs84") {
      ncatt_put(nc, crsVar, attname="grid_mapping_name", attval = "latitude_longitude")
      ncatt_put(nc, crsVar, attname="longitude_of_prime_meridian", attval=0.0)
      ncatt_put(nc, crsVar, attname="semi_major_axis", attval=6378137.0)
      ncatt_put(nc, crsVar, attname="inverse_flattening", attval=298.257223563)
   } else {
      stop("CRS unknown")
   }
   
   # specify the convention as a global attribute 
   ncatt_put(nc, 0, attname="Conventions", attval="CF-1.5")

   nc_close(nc)
   
   return(invisible(file))
}

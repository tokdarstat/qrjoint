\name{redmaple}
\alias{redmaple}
\docType{data}
\title{Basal Areas of Red Maple Trees}
\description{
Site-level basal areas of red maple trees (\emph{Acer rubrum}) at 608 unmanaged and forested sites in Connecticut, Massachusetts, and Rhode Island. Data are aggregated from the Forest Inventory Analysis (FIA) of the US Forest Service. Geographical regions of the sites are added using Enviromental Protection Agency (EPA) shapefiles.
}
\usage{data(redmaple)}
\format{
  A data frame with 608 observations on the following variables for each site:
  \describe{
    \item{\code{plotID}}{Unique identifier}
    \item{\code{elev}}{Elevation, measured in feet}
    \item{\code{slope}}{Slope, measured in degrees}
    \item{\code{aspect}}{Aspect, measured in degrees proceeding from North clockwise around a compass. For sites with zero or near-zero slopes, aspect is recorded as 0. North is recorded as 360.}
    \item{\code{Region}}{EPA Level-III geographical region}
    \item{\code{region}}{EPA Level-III geographical region, shortened}
    \item{\code{State}}{State}
    \item{\code{baRedMaple}}{Total basal area of red maple trees, measured in square feet per acre}
  }
}

\details{
This three-state subset from the FIA is intended to illustrate the capabilities of the qrjoint package in flexibly modeling excess-boundary zeros, using its censoring option. This subset of variables should not be construed as a comprehensive list of factors influencing red maple basal area growth. All sites in the sample are of equivalent area.
}

\source{Forest Inventory and Analysis Database, St. Paul, MN: U.S. Department of Agriculture, Forest Service, Northern Research Station. \url{https://apps.fs.usda.gov/fia/datamart/datamart.html}}

\examples{
data(redmaple)
}

\keyword{datasets}

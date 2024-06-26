% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PPutilities.R
\name{loo_mean}
\alias{loo_mean}
\title{efficient leave one out mean and standard deviation}
\usage{
loo_mean(var, year, ...)
}
\arguments{
\item{var}{column name of the variance to take the mean of.}

\item{year}{just pass year as in the example above in order to access the year-column in the data table.}

\item{...}{passed on to mean (or sd)}
}
\description{
These only work if, after grouping, you have only one observation per year. When this is the case, loo-mean and sd can be computed
without looping over years.
This is typically the case if you consider ensemble means, but not if one column is ensemble member.
In this case, if you e.g. want to derive the climatology of an nwp system, you want to compute the mean accross all members leaving out all members
of the year under consideration. For doing this correctly you should use the slower loyo instead!
loo_mean and loo_sd should be used inside a data table, e.g. dt\link{,clim := loo_mean(obs,year)}
}

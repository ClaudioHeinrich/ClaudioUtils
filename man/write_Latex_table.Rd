% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_utilities.R
\name{write_Latex_table}
\alias{write_Latex_table}
\title{Write a LaTeX table from a matrix}
\usage{
write_Latex_table(mat, fn = "", dec = ".", na = "-", replace_underscores = " ")
}
\arguments{
\item{mat}{The matrix to be written as a LaTeX table.}

\item{fn}{Optional. The filename to save the LaTeX table. Defaults to print in console.}

\item{dec}{Optional. The decimal separator to use. Default is ".".}

\item{na}{Optional. The character to represent missing values. Default is "-".}
}
\value{
None
}
\description{
This function writes a LaTeX table from a given matrix.
}
\details{
This function writes the provided matrix as a LaTeX table by using the write.table function.
The resulting LaTeX table is saved in a file specified by the \code{fn} parameter. If no filename is provided,
the LaTeX table is returned.

Note: If there is a file named "temp.txt" in the working directory, this function will stop execution and throw an error.
}

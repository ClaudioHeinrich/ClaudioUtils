
#' Function to convert an array to a data table
#'
#' Lists an array of any dimension as data table, with one column per dimension and one column with the array values.
#' A particular feature is that dimensions can be passed as names of vector-type objects in the parent environment, which is frequently useful, see examples.
#' That only works by choosing object names that are not also the names of imported functions. E.g, you shouldn't call an object 'date'.
#'
#' @param ar array.
#' @param dims Either a list of vector-type objects giving the dimensions, or a vector of names of such objects in the parent environment.
#' @param val_cn Name of the new column
#' @examples
#' dates = as.Date(c('2024-01-01','2020-01-02'))
#' x = 1:3
#' y = 9:10
#'
#' array_data = runif(2*3*2) |> array(dim = c(2,3,2))
#'
#' #call without dimension-variables:
#' dt1 = cube_to_dt(array_data)
#' print(dt1)
#'
#' #call with dimensions given as list:
#' dimlist = list(date = dates,x = x,y = y)
#' dt2 = cube_to_dt(array_data, dims = dimlist)
#' print(dt2)
#'
#' # call with names of vector-objects
#' dt3 = cube_to_dt(array_data, dims = c('dates','x','y'))
#' print(dt3)
#'
#' @export

cube_to_dt = function(ar,dims = NULL,val_cn = 'value'){

  ndims = length(dim(ar))
  # get dimension list if empty:
  if(is.null(dims)){
    dims_new = list()
    for(i in 1:ndims)# what if ar is a 1d-array
    {
      temp_list = list(1:dim(ar)[i])
      names(temp_list) = paste0('Var',i) # use 'Var' in order to be consistent with expand.grid behavior for unnamed lists.
      dims_new = c(dims_new,temp_list)
    }
  } else if(is.character(dims)){
    # get dimension list if provided as names of objects in parent environment
    dims_new = list()
    for(i in seq_along(dims)){
      temp_obj = get(dims[i],envir = parent.env(as.environment(-1)))
      if(!("vector" %in% is(temp_obj) | lubridate::is.Date(temp_obj))){ # is.vector does somehow not recognize arrays that are also vectors, as returned from ncvar_get
        stop(paste0('Found object ',dims[i],' in parent environment, but it is not a vector.'))
      }else if(length(temp_obj) != dim(ar)[i]){
        stop(paste0('Found length-',length(temp_obj),' vector ',dims[i],' in parent environment, but dimension ',i,' of the provided array has length ',dim(ar)[i],'.'))
      }
      temp_list = list( temp_obj)
      names(temp_list) = dims[i]
      dims_new = c(dims_new,temp_list)
    }
  } else{
    dims_new = dims
  }

  dt_temp = as.data.table(expand.grid(dims_new))
  dt_temp[,new_col := as.vector(ar)]
  setnames(dt_temp,'new_col',val_cn)

  return(dt_temp[])
}


#' Write a LaTeX table from a matrix
#'
#' This function writes a LaTeX table from a given matrix.
#'
#' @param mat The matrix to be written as a LaTeX table.
#' @param fn Optional. The filename to save the LaTeX table. Defaults to print in console.
#' @param dec Optional. The decimal separator to use. Default is ".".
#' @param na Optional. The character to represent missing values. Default is "-".
#'
#' @return None
#'
#' @details This function writes the provided matrix as a LaTeX table by using the write.table function.
#' The resulting LaTeX table is saved in a file specified by the \code{fn} parameter. If no filename is provided,
#' the LaTeX table is returned.
#'
#' Note: If there is a file named "temp.txt" in the working directory, this function will stop execution and throw an error.
#'
#'
#'
#' @export
write_Latex_table = function(mat,fn = "",dec = ".",na = "-",replace_underscores = ' ')
{
  rn = !is.null(rownames(mat))

  if(file.exists("temp.txt")) stop('There is a file "temp.txt" in the work directory which would be overwritten.')

  # write out
  write.table(mat, "temp.txt", quote = FALSE, append = FALSE, sep = "  &  ", dec = dec,na = na,eol = " \\\\\n",
              row.names = TRUE, col.names = TRUE)

  temp = readLines("temp.txt")
  file.remove("temp.txt")
  # remove \\ from last line:
  temp[length(temp)] = gsub('\\\\','',temp[length(temp)])
  temp = paste0(temp,collapse = '\n')
  if(rn) tttemp = paste0('\\begin{tabular}{',paste0(rep('c ',ncol(mat) + 1),collapse = ' '),'}\n&',temp,'\n\\end{tabular}')
  if(!rn) tttemp = paste0('\\begin{tabular}{',paste0(rep('c ',ncol(mat)),collapse = ' '),'}\n',temp,'\n\\end{tabular}')
  if(!identical(replace_underscores,FALSE)) tttemp = gsub('_',replace_underscores,tttemp)
  if(fn == "") fn = stdout()
  writeLines(tttemp,fn)
}

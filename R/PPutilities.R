#' fast full-grid bilinear interpolation for different values between the same grid
#'
#' @description Interpolates multiple values in DT onto the grid specified by the coordinates
#' in DT_new. Standard bilinear interpolation functions are looking up the coordinates
#' for every value. This is inefficient if you want to interpolate multiple variables
#' between THE SAME grids, since we need to look up the reference coordinates only once,
#' which is done by this function. A typical scenario is when you have an ensemble forecast
#' with K members and D lead times, and map it onto a new grid. The origin grid is assumed to be
#' regular in the sense that we have lons lon_1,...,lon_n, and lats lat_1,...,lat_m, and the
#' full grid is given by {(lon_i,lat_j)}_{i in {1,...,n}, j in {1,...,m}}
#'
#'  @param DT data table containing old coordinates lon and lat, the value variables to be interpolated
#'  and expand variables, for each level of which the interpolation should be conducted.
#'  In the example above, the expand variables would be 'ensemble member' and 'date'.
#'  @param expand_variable_cols Character vector containing the column names of the expand variables.
#'  @param function_value_cols Character vector containing the column names of the
#'  variables that should be interpolated. For example c('temperature','precip')
#'  @param DT_new Data table containing the new coordinates
#'
#'  @return a data table
#'  @author Claudio
#'  @export


bil_interpol_grid_dt = function(DT,expand_variable_cols,function_value_cols,DT_new)
{

  lon_lat_new = unique(DT_new[,.(lon,lat)])
  ###################

  setkeyv(DT,c(expand_variable_cols,'lon','lat'))

  lon_lat_old = unique(DT[,.(lon,lat)])
  lon_lat_old[,row_index := 1:.N]

  nl_old = lon_lat_old[,.N]
  DT[,row_index := rep(1:nl_old,.N/nl_old)]

  if(! identical(unique(DT[,.(lon,lat,row_index)]),lon_lat_old))
  {
    stop('something is weird with the data, maybe different locations for different levels of the expand variables?')
  }


  # efficient way of finding the square in the old grid containing the new grid

  lons_old = sort(unique(lon_lat_old[,lon]))
  lons_new = unique(lon_lat_new[,lon])

  lats_old = sort(unique(lon_lat_old[,lat]))
  lats_new = unique(lon_lat_new[,lat])


  # check for regular grid:
  if(length(lons_old) * length(lats_old) != lon_lat_old[,.N])
  {
    stop('Only regular grids are supported.')
  }

  # find for each grid point in the new data set the longitude and latitude just below and above it

  lons_comparison = rep(lons_old,length(lons_new)) - rep(lons_new,each = length(lons_old))
  lons_comparison = matrix(lons_comparison,nrow = length(lons_old))

  lon_below_ind = apply(lons_comparison,MARGIN = 2, FUN = function(x){max(which(x <= 0))})
  lon_above_ind = apply(lons_comparison,MARGIN = 2, FUN = function(x){min(which(x > 0))})

  # pad locations that are outside the old grid:
  lon_below_ind[is.infinite(lon_below_ind)] = 1 # minimum longitude index
  lon_above_ind[is.infinite(lon_above_ind)] = length(lons_old) # maximum longitude index

  ####

  lon_below = lons_old[as.vector(lon_below_ind)]
  lon_above = lons_old[as.vector(lon_above_ind)]

  #### same for latitudes: ####

  lats_comparison = rep(lats_old,length(lats_new)) - rep(lats_new,each = length(lats_old))
  lats_comparison = matrix(lats_comparison,nrow = length(lats_old))

  lat_below_ind = apply(lats_comparison,MARGIN = 2, FUN = function(x){max(which(x <= 0))})
  lat_above_ind = apply(lats_comparison,MARGIN = 2, FUN = function(x){min(which(x > 0))})

  # pad locations out of old grid:
  lat_below_ind[is.infinite(lat_below_ind)] = 1 # minimum latitude index
  lat_above_ind[is.infinite(lat_above_ind)] = length(lats_old) # maximum latitude index

  lat_below = lats_old[as.vector(lat_below_ind)]
  lat_above = lats_old[as.vector(lat_above_ind)]

  ########################

  dt_temp = lon_lat_new

  # find Q11-indices (for each new coordinate, the lower left coordinate of the old-coordinate-rectangle containing it)

  Q11_dt = data.table(as.data.table(expand.grid(lon = lon_below,lat = lat_below)),as.data.table(expand.grid(lon_new = lons_new,lat_new = lats_new)))
  setkey(Q11_dt,lon,lat) # required for the following subsetting:

  setkey(DT,lon,lat)
  Q11_dt = merge(DT,Q11_dt, allow.cartesian = TRUE)

  setnames(Q11_dt,function_value_cols,paste0('Q11',function_value_cols))
  setnames(Q11_dt,c('lon','lat'),c('Q11lon','Q11lat'))

  # find Q12-indices

  Q12_dt = data.table(as.data.table(expand.grid(lon = lon_below,lat = lat_above)),as.data.table(expand.grid(lon_new = lons_new,lat_new = lats_new)))
  setkey(Q12_dt,lon,lat) # required for the following subsetting:

  Q12_dt = merge(DT,Q12_dt, allow.cartesian = TRUE)

  setnames(Q12_dt,function_value_cols,paste0('Q12',function_value_cols))
  setnames(Q12_dt,c('lon','lat'),c('Q12lon','Q12lat'))

  # find Q21-indices

  Q21_dt = data.table(as.data.table(expand.grid(lon = lon_above,lat = lat_below)),as.data.table(expand.grid(lon_new = lons_new,lat_new = lats_new)))
  setkey(Q21_dt,lon,lat) # required for the following subsetting:

  Q21_dt = merge(DT,Q21_dt, allow.cartesian = TRUE)

  setnames(Q21_dt,function_value_cols,paste0('Q21',function_value_cols))
  setnames(Q21_dt,c('lon','lat'),c('Q21lon','Q21lat'))

  # find Q22-indices

  Q22_dt = data.table(as.data.table(expand.grid(lon = lon_above,lat = lat_above)),as.data.table(expand.grid(lon_new = lons_new,lat_new = lats_new)))
  setkey(Q22_dt,lon,lat) # required for the following subsetting:

  Q22_dt = merge(DT,Q22_dt, allow.cartesian = TRUE)

  setnames(Q22_dt,function_value_cols,paste0('Q22',function_value_cols))
  setnames(Q22_dt,c('lon','lat'),c('Q22lon','Q22lat'))

  # get data table containing everything necessary:

  DT_new = data.table(Q11_dt,
                      Q12_dt[,.SD,.SDcols = c(paste0('Q12',function_value_cols),'Q12lon','Q12lat')],
                      Q21_dt[,.SD,.SDcols = c(paste0('Q21',function_value_cols),'Q21lon','Q21lat')],
                      Q22_dt[,.SD,.SDcols = c(paste0('Q22',function_value_cols),'Q22lon','Q22lat')]
  )

  # get the bilinear interpolation value, see wikipedia
  for(fv in function_value_cols)
  {
    factor1 = DT_new[,1/((Q22lon - Q11lon) * (Q22lat - Q11lat))]
    f11 = DT_new[,(Q22lon - lon_new)*(Q22lat - lat_new) * get(paste0('Q11',fv))]
    f21 = DT_new[,(lon_new - Q11lon)*(Q22lat - lat_new) * get(paste0('Q21',fv))]
    f12 = DT_new[,(Q22lon - lon_new)*(lat_new - Q11lat) * get(paste0('Q12',fv))]
    f22 = DT_new[,(lon_new - Q11lon)*(lat_new - Q11lat) * get(paste0('Q22',fv))]
    DT_new[,(fv) := factor1 * (f11 + f21 + f12 + f22)]

  }

  DT_new = DT_new[,.SD,.SDcols = c(expand_variable_cols,'lon_new','lat_new',function_value_cols)]

  setnames(DT_new,c('lon_new','lat_new'),c('lon','lat'))
  setkeyv(DT_new,c(expand_variable_cols,c('lon','lat')))

  return(DT_new)
}



bootstrap_MSE = function(prec_dt,pred_colnames,obs_colname = 'obs',bycols = c('lon','lat','year'),RR = 1000)
{
  # get errors for all predictive models

  if('member' %in% colnames(prec_dt))
  {
    temp_dt = prec_dt[,lapply(.SD,mean,na.rm = T),
                      .SDcols = c(pred_colnames,obs_colname),
                      by = bycols]
  }

  MSEs = temp_dt[,lapply(.SD,function(x){(x-get(obs_colname))^2}),
                 .SDcols = pred_colnames,
                 by = bycols]


  ### bootstrap everything ###

  statistic = function(x,inds){mean(x[inds],na.rm = T)}

  bootstrap_all_region_all_years = as.data.table(expand.grid(model = pred_colnames,R = 1:RR))

  for(pred_mod in pred_colnames)
  {
    print(pred_mod)
    data = MSEs[,get(pred_mod)]
    bootstrap_all_region_all_years[model == pred_mod,boots := boot::boot(data,statistic = statistic,R = RR)$t]
  }

  return(bootstrap_all_region_all_years)
}


#' Claudios path to SFE directory
#' @export
claudio_sfe_dir = function(){
  return("/nr/user/claudio/bigdisk/SFE/")
}


#' @export
global_confer_lat_subset = function(){
  lat_subset = seq(-25,35,0.5)
  return(lat_subset)
}


#' @export
global_confer_lon_subset = function(){
  lon_subset = seq(5,65,0.5)
  lon_subset[lon_subset < 0] = lon_subset[lon_subset < 0] + 360
  return(lon_subset)
}


#' specification of smallest grid fully containing the greater-horn-of-Africa-region, as considered in CONFER
#' @export
global_gha_lat_subset = function(){
  lat_subset = seq(-12,23.5,0.5)
  return(lat_subset)
}

#' specification of smallest grid fully containing the greater-horn-of-Africa-region, as considered in CONFER
#' @export
global_gha_lon_subset = function(){
  lon_subset = seq(21.5,51.5,0.5)
  return(lon_subset)

}




#' @export
get_cov_mat = function(dt,cn, bycols = c('lon','lat'), levelcols = c('year','member'),collapse = F)
{
  temp_dt = dt[,.SD,.SDcols = c(bycols,levelcols,cn)]
  if(collapse) temp_dt = unique(temp_dt)

  #temp_dt = temp_dt[!is.na(get(cn))]
  temp_dt[,paste0(cn,'_mu') := mean(get(cn)), by = bycols]
  temp_dt[,(cn) := get(cn) - get(paste0(cn,'_mu'))]
  temp_dt[,paste0(cn,'_mu') := NULL]
  setkeyv(temp_dt,bycols)

  p = unique(temp_dt[,.SD,.SDcols = bycols])[,.N]
  n = unique(temp_dt[,.SD,.SDcols = levelcols])[,.N]
  if(p*n != temp_dt[,.N]) stop('dimensions are not matching! If you want to collapse along a column, set collapse to TRUE (e.g. if you consider observations for a data table that contains a member-column)')
  mm = matrix(temp_dt[,get(cn)],byrow = T,nrow = p)

  cvm = 1/(n-1) * mm %*% t(mm)

  return(cvm)
}


#' @export

get_spatial_EOFs = function(dt,cn,bycols = c('lon','lat'), levelcols = intersect(c('year','member'),names(dt)),collapse = F,ret_irlba = F,...)
{
  cov_mat = get_cov_mat(dt,cn,bycols,levelcols,collapse)

  tr = sum(diag(cov_mat))

  svd = irlba(cov_mat,...)


  percentages_of_variance = 100*sum(svd$d)/tr
  print(paste0('The ',length(svd$d),' returned EOFs explain ',paste0(round(percentages_of_variance,1),sep = ', '),'% of the variance.'))

  ret_dt = unique(dt[,.SD,.SDcols = bycols])
  setkeyv(ret_dt,bycols)

  EOF_dt = as.data.table(svd$u)
  setnames(EOF_dt,paste0(cn,'_eof_',1:length(svd$d)))

  if(!ret_irlba) return(data.table(ret_dt,EOF_dt))
  if(ret_irlba) return(svd)
}



#' Get the Greater Horn of Africa Region considered in CONFER
#'
#'
#' @param resolution Either 'half_deg' or 'full_deg', depending on whether you want all half-degree-gridpoints in the region or only full-degree.
#' @param file_name where is the ICPAC_region information stored?
#'
#' @examples
#' \dontrun{ #for subsetting a data table dt to the ICPAC locations do this:
#' setkey(dt,lon,lat)
#' dt_sub = dt[GHA_locs()]}
#'
#' @export


GHA_locs = function(resolution = 'half_deg',file_name = '/nr/project/stat/CONFER/Data/ICPAC_region.RData')
{
  load(file_name)
  if(resolution == 'half_deg')
  {
    return(half_deg_locs)
  }
  if(resolution == 'full_deg'){
    return(full_deg_locs)
  }
}



#' convenient wrapper for loading predictions and observations into a data table
#'
#' WARNING: This version is outdated and should not be used, use load_cds_monthly_data instead!
#'
#' #'
#' @param var_name short name for variable to be loaded
#' @param area should the variable be only returned for a certain area, e.g. 'ICPAC' returns it only for the ICPAC region (still, the entire globe is loaded, so it doesn't speed the function up).
#' @param system_name For 'era' the ERA-reanalysis data is loaded, otherwise the prediction from the NWP with corresponding name is loaded.
#' @param target_mons which months do you want?
#' @param init_mons If system_name !- 'era', which initialization months for the predictions do you want to consider?
#' @param load_dir path to the files.
#' @param resolution 'half_degree' or 'full_degree'

load_dt_monthly = function(var_name = 'prec',area = NULL,system_name = 'ERA',target_mons = NULL,init_mons = NULL,
                           load_dir = ifelse(system_name %in% c('ERA','era'),
                                             yes = '~/bigdisk/SFE/ERA_monthly_compiled/',
                                             no = '~/bigdisk/CONFER/Systems_monthly_compiled/'),
                           resolution = 'half_degree')
{
  # translate short variable names into the variable name used in the file name:
  if(var_name %in% c('prec','precipitation','Prec','Precipitation')) var_name = 'total_precipitation'
  if(var_name %in% c('sst','SST')) var_name = 'sea_surface_temperature'

  # define area:
  if(is.null(area))
  {
    if(resolution == 'full_degree')
    {
      subset_locs = as.data.table(expand.grid(lon = -179:180,lat = -90:90))
    }
    if(resolution == 'half_degree')
    {
      subset_locs = as.data.table(expand.grid(lon = seq(-179.5,180,0.5),lat = seq(-90,90,0.5)))
    }
    setkey(subset_locs,lon,lat)

  }else if(area == 'ICPAC')
  {
    load('/nr/project/stat/CONFER/Data/ICPAC_region.RData')

    if(resolution == 'full_degree')
    {
      subset_locs = full_deg_locs
    }
    if(resolution == 'half_degree')
    {
      subset_locs = half_deg_locs
    }
  }else if(area == 'Nino3.4')
  {
    if(resolution == 'full_degree')
    {
      subset_locs = as.data.table(expand.grid(lat = -5:5, lon = (-170:-120)))
    }
    if(resolution == 'half_degree')
    {
      subset_locs = as.data.table(expand.grid(lat = seq(-5,5,0.5), lon = seq(-170,-120,0.5)))
    }
    setkey(subset_locs,lon,lat)
  }else if(area == 'IOD')
  {
    if(resolution == 'full_degree')
    {
      subset_locs1 = as.data.table(expand.grid(lat = -10:10, lon = (50:70)))
      subset_locs2 = as.data.table(expand.grid(lat = -10:0, lon = (90:110)))
      subset_locs = rbindlist(list(subset_locs1,subset_locs2))
    }
    if(resolution == 'half_degree')
    {
      subset_locs1 = as.data.table(expand.grid(lat = seq(-10,10,0.5), lon = seq(50,70,0.5)))
      subset_locs2 = as.data.table(expand.grid(lat = seq(-10,0,0.5), lon = seq(90,110,0.5)))
      subset_locs = rbindlist(list(subset_locs1,subset_locs2))
    }
    setkey(subset_locs,lon,lat)
  }

  if(is.null(target_mons)) target_mons = 1:12

  if(system_name %in% c('ERA','era'))
  {
    obs_dt = data.table()
    for(mm in target_mons)
    {
      print(mm)
      temp = fread(paste0(load_dir,var_name,'_era_compiled_',mm,'.csv'))
      temp[lon > 180,lon:=lon - 360]
      setkey(temp,lon,lat)
      obs_dt = rbindlist(list(obs_dt,temp[subset_locs,][,month := mm]))
      rm(temp)
      gc()
    }

    obs_dt[lon > 180,lon:=lon - 360]

    obs_dt[,obs := climatology + anomaly]

    return(obs_dt)
  }

  if(!(system_name %in% c('ERA','era')))
  {
    pred_dt = data.table()
    for(in_mm in init_mons)
    {
      for(tar_mm in target_mons)
      {
        print(c(in_mm,tar_mm))
        fn = paste0(load_dir,var_name,'_',system_name,'_',in_mm,'_',tar_mm,'.csv')
        if(file.exists(fn))
        {
          temp = fread(fn)
          temp[lon > 180,lon:=lon - 360]
          setkey(temp,lon,lat)
          pred_dt = rbindlist(list(pred_dt,temp[subset_locs,][,month := tar_mm]))
          gc()
        }

      }
    }
    return(pred_dt)
  }

}


#### leave-one-out tools: ####

#' efficient leave one out mean and standard deviation
#'
#' These only work if, after grouping, you have only one observation per year. When this is the case, loo-mean and sd can be computed
#' without looping over years.
#' This is typically the case if you consider ensemble means, but not if one column is ensemble member.
#' In this case, if you e.g. want to derive the climatology of an nwp system, you want to compute the mean accross all members leaving out all members
#' of the year under consideration. For doing this correctly you should use the slower loyo instead!
#' loo_mean and loo_sd should be used inside a data table, e.g. dt[,clim := loo_mean(obs,year)]
#' @param var column name of the variance to take the mean of.
#' @param year just pass year as in the example above in order to access the year-column in the data table.
#' @param ... passed on to mean (or sd)
#'
#' @export

loo_mean = function(var,year,...)
{
  ny = length(unique(year))
  return(ny/(ny-1) * (mean(var,...) - var/ny))
}

#' same as loo_mean, only for standard deviation
#' @export

loo_sd = function(var,year,...)
{
  ny = length(unique(year))
  loo_sq_sum = sum(var^2,...) - var^2
  loo_sum = sum(var,...) - var

  remove_negatives = (loo_sq_sum - loo_sum^2/(ny - 1)) # this should never be negative
  remove_negatives[remove_negatives <0] = 0
  ret_val = sqrt(1/(ny-2)*(loo_sq_sum - loo_sum^2/(ny - 1)))
  return(ret_val)
}

#' leave-one-year-out application of any function
#'
#' new columns are generated in the data table containing the leave-one-out-values. Requires a function with only one argument. If you require multiple arguments use loyoma.
#' @param dt Data table with columns 'year', SDcols and bycols
#' @param FUN function to apply
#' @param SDcols character vector containing the column names of the columns FUN sould be applied to.
#' @param bycols character vector containing the column names of the grouping variables, pass NULL if you don't want to group.
#' @param mc.cores how many cores to use in parallelization.
#' @param newnamesstring The results of the leave-one-out computation are appended to dt as columns with the names \code{SDcols}_\code{newnamesstring}.
#' @param args Character vector. Allows to pass additional columns as arguments to the function. The argument name within FUN needs to match the column name.
#' @param overwrite Logical. If you want to use loyo to overwrite existing columns in the data table, this needs to be set to TRUE.
#' @param ... Arguments passed on to FUN.
#'
#' @return If bycols is not NULL, the data table dt with appended columns, else a data table with the SDcols and a singlerow for each year
#'
#' @import data.table
#' @importFrom parallel mclapply
#' @export

loyo = function(dt, FUN,SDcols,bycols = c('lon','lat'),mc.cores = 4,newnamesstring = 'new',args = NULL,overwrite = F,...)
{
  if(identical(bycols,'')) bycols = NULL

  # duplicate columns in dt are an issue for merging, and should not be there anyway:
  if(!identical(colnames(dt),unique(colnames(dt))))
  {
       dt = dt[, .SD, .SDcols = unique(names(dt))]
       warning('dt had duplicate column names. The duplicates were removed (by selecting the first column).')
  }

  owcols = intersect(paste0(SDcols,'_',newnamesstring),colnames(dt))
  if(length(owcols) >0)
  {
    if(overwrite)
    {
      dt[,(owcols):= NULL]
    } else {
      stop("This would overwrite existing columns in the data table. If that's ok, set overwrite to TRUE.")
    }
  }

  apply_by_year = function(yy)
  {
    dt_temp = dt[year != yy]
    if(!is.null(args))
    {
      argstring = paste(args,args,sep = ' = ', collapse = ', ')
      lapplyparse = paste0('lapply(X = .SD,FUN = FUN,',argstring,',...)')
    } else {
      lapplyparse = 'lapply(X = .SD,FUN = FUN,...)'
    }
    dt_temp = dt_temp[,eval(parse(text = lapplyparse)),.SDcols = SDcols,by = bycols]
    setnames(dt_temp,SDcols,paste0(SDcols,'_',newnamesstring))

    if(is.null(bycols))
    {
      dt_new = data.table(dt[year == yy],dt_temp)
    } else{
      dt_new = merge(dt[year == yy],dt_temp,by = bycols)
    }

    return(dt_new)
  }

  ret_dt = rbindlist(parallel::mclapply(X = unique(dt[,year]),
                                        FUN = apply_by_year,
                                        mc.cores = mc.cores))
  return(ret_dt)
}





#' leave-one-year-out standardization
#'
#' Takes a data table and standardizes one or several columns in a leave-one-year-out fashion.
#' For standardization it is important to use leave-one-year-out if one of your columns is 'member', otherwise your mean estimate is in-sample.
#'
#' @param dt Data table with columns 'year', SDcols and bycols
#' @param SDcols character vector containing the column names of the columns that should be standardized.
#' @param bycols character vector containing the column names of the grouping variables, pass NULL if you don't want to group.
#' If you have several ensemble members, do not group by member!
#' @param mc.cores how many cores to use in parallelization.
#' @param newnamesstring Either NULL or character. If NULL, a data table of the same form as dt is returned with the SDcols
#' overwritten by standardized versions (default). If character, new columns are attached and their names are
#' \code{paste0(SDcols,'_',newnamesstring)}.
#' @param na.rm Logical. Should NAs be removed when the means and variances are computed?
#'
#' @return a copy of dt with standardized columns.
#' @export
#'
#' @author Claudio

loyo_standardize = function(dt,
                            SDcols,
                            bycols = c('lon','lat'),
                            mc.cores = 4,
                            newnamesstring = NULL,
                            na.rm = T)
{
  dt_temp = copy(dt)

  # get loyo mean and sd

  dt_temp = loyo(dt_temp,
                 FUN = mean,
                 SDcols = SDcols,
                 bycols = bycols,
                 mc.cores = mc.cores,
                 newnamesstring = 'mean',
                 overwrite = T,
                 na.rm = na.rm)

  dt_temp = loyo(dt_temp,
                 FUN = sd,
                 SDcols = SDcols,
                 bycols = bycols,
                 mc.cores = mc.cores,
                 newnamesstring = 'sd',
                 overwrite = T,
                 na.rm = na.rm)

  for(cn in SDcols)
  {
    standardized_vals = dt_temp[,(get(cn) - get(paste0(cn,'_mean')))/get(paste0(cn,'_sd'))]
    if(is.null(newnamesstring)) dt_temp[,(cn) := standardized_vals]
    if(!is.null(newnamesstring)) dt_temp[,(paste0(cn,'_',newnamesstring)) := standardized_vals]
  }

  trashcols = c(paste0(SDcols,'_mean'),paste0(SDcols,'_sd'))
  dt_temp[,(trashcols) := NULL]

  return(dt_temp)
}



#' function for masking out areas with low precip
#' @param dt data table containing precip data
#' @param threshold areas with average precip below threshold are masked out
#' @param value_col name of the column where precip is stored.
#' @param bycols masking is done grouped by these columns (should at least contain 'lon','lat')
#' @return data table dt, with a logical column called 'mask'.
#'
#' @author Claudio
#' @export
#'

mask_precip = function(dt,threshold = 0.2,value_col = 'prec',bycols = intersect(c('month','lon','lat'),names(dt)))
{
  dt[,mask := mean(get(value_col))<=threshold,by = bycols]
  return(dt)
}





MSEdiffplot = function(dt,col1,col2,obscolname = 'obs',rr = NULL)
{
  dt_temp = dt[,lapply(.SD,mean,na.rm = T),.SDcols = c(col1,col2,obscolname),by = c('year','lon','lat')]

  MSEs = dt_temp[,lapply(.SD,function(x){mean((x - get(obscolname))^2)}),.SDcols = c(col1,col2),by = c('lon','lat')]

  MSEs[,MSE_diff := get(col1) - get(col2)]
  if(is.null(rr))
  {
    rr = max(abs(MSEs[,range(MSE_diff)]))
    rr = c(-rr,rr)
  }
  plot_diagnostic(MSEs,'MSE_diff',mn = paste0('MSE difference ',col1,' - ',col2),rr = rr)
}




#' function for running a permutation test
#' @param a scores of first model
#' @param b scores of second model
#' @param N numebr of permutations used
#'
#' @return list with three elements: d_bar is the mean score difference, D is a vector with the N mean score differences,
#'         and p_val is the p-value of the permutation test for the hypothesis that a is worse or equally good as b.
#'
#'
#'
#' @author Alex
#'
#' @export


permutation_test_difference = function(a,
                                       b,
                                       N = 5e3){
  n = length(a)
  d = a - b
  d_bar = mean(d)
  D = NULL
  for(i in 1:N){
    swap = rbinom(n,1,0.5)
    w_swap = which(swap == 1)
    d_i = d
    d_i[w_swap] = -d_i[w_swap]
    D[i] = mean(d_i)
  }

  p_val = sum(d_bar > sort(D))/N + sum(d_bar == sort(D))/(2*N)

  return(list(d_bar = d_bar, D = D,p_val = p_val))
}


#' function for reducing a data table containing ensemble predictions to the ensemble means:
#' @param dt data table, needs to contain a column called 'member'
#' @param pred_vars character vector containing the names of the columns containing the ensemble predictions.

take_ens_mean = function(dt,pred_vars)
{
  if(!('member' %in% colnames(dt))) stop('dt needs to contain a column called "member".')
  other_colnames = setdiff(names(dt),c(pred_vars,'member'))
  ret_dt = dt[,lapply(.SD,mean,na.rm = T),.SDcols = pred_vars,by = other_colnames]
  return(ret_dt)
}

#' function for time-lagging monthly data
#' @param dt a data table with years and months.
#' @param lag how many months do you want to lag?
#'
#' @return The lagged data table.
#' @author Claudio
#'
#' @export


time_lag_monthly_dt = function(dt,lag = 0)
{
  dt[,month := month + lag]

  while(dt[month < 1,.N] > 0)
  {
    dt[month < 1, year := year - 1]
    dt[month < 1, month := month + 12]
  }

  while(dt[month > 12,.N] > 0)
  {
    dt[month > 12, year := year + 1]
    dt[month > 12, month := month - 12]
  }

  return(dt)

}

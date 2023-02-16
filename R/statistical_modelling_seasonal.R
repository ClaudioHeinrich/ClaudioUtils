

############## Contains functions useful for (spatial) regression modelling using data tables. #############




#' @description Runs a (generalized) linear model on a data table. Allows grouping without looping.
#'
#' Simply runs a glm on a data table containing all covariates and the dependent variable. Allows to group by columns in the data table.
#' For example you can regress precip on pred_precip and estimate new regression coefficients for each location by running 
#' \code{glm_dt(dt,prec ~ pred_prec, bycols = c('lon','lat'))}. This function wraps \code{glm} but can also be used to fit linear models.
#'
#'  @param dt a data table with columns year, month, lon, lat 
#'  @param formula formula for glm, all covariates need to be column names in dt.
#'  @param ...  passed on to glm. In particular family is passed here when the glm is nonlinear.
#'  @param bycols NULL or a character vector containing the column names of dt to group by: A glm is fit for each level of the grouping variables.
#'  @param coefficients_only logical. If TRUE, a data table containing the coefficients is returned. Otherwise a data table with 
#'  coefficients and predictions is returned.
#'  @param alongcols Only required if !coefficients_only. Which columns contain the dimension variables separating the observation used in the same regression.
#'  In the example above this might be 'year', for example.
#'  
#'  
#'  @return See coefficients_only.
#'    
#'  @author Claudio
#'  @export


glm_dt = function(dt,formula,...,
                  bycols = NULL,
                  coefficients_only = T,
                  alongcols = 'year')
{
  if(!coefficients_only)# crude check if alongcols and bycols span all coordinates in dt
  {
    check = dt[,.N] == unique(dt[,.SD,.SDcols = c(bycols,alongcols)])[,.N]
    if(!check & !coefficients_only) stop('If you want predictions, alongcols and bycols together need to span the data table.')
  }
  
  
  allvars = all.vars(formula)
  response = allvars[1]
  
  if(coefficients_only)
  {
    grab_from_glm = function(gg)
    {
      new_list = as.list(coefficients(gg))
      coeff_names = names(new_list)[!(names(new_list) == '(Intercept)')] 
      new_coeff_names = unlist(lapply(coeff_names,function(x){return(paste0(x,'_coef'))}))
      names(new_list)[!(names(new_list) %in% c('(Intercept)'))] = new_coeff_names
      return(as.data.table(new_list))
    }
  } else {
    grab_from_glm = function(gg)
    {
      new_list = c(as.list(coefficients(gg)),
                   list(fitted = as.numeric(fitted(gg))))
      coeff_names = names(new_list)[!(names(new_list)%in% c('(Intercept)','fitted'))] 
      new_coeff_names = unlist(lapply(coeff_names,function(x){return(paste0(x,'_coef'))}))
      names(new_list)[!(names(new_list)%in% c('(Intercept)','fitted'))] = new_coeff_names
      return(as.data.table(new_list))
    }
  }
  
  setkeyv(dt,bycols)
  
  dt_glm = dt[,grab_from_glm(glm(data = .SD, formula = formula)),
              by = bycols]
    

  if(!coefficients_only)
  {
    dt_glm = data.table(dt_glm, dt[,.SD,.SDcols = alongcols])
    dt_glm[,(allvars) := dt[,.SD,.SDcols = allvars]]
  }
  
  return(dt_glm)
}


#' @description Fit a linear model to a data table in a leave-one-year-out fashion. Allows grouping without looping, see glm_dt
#' 
#' This function currently assumes that for all grouping levels we have the same number of observations. 
#' This is sometimes violated, e.g. when you have several systems and some of them have missing years. In this case the function will
#' stop.
#' @param dt a data table. Needs to contain all covariates as columns and a column 'year'.
#' @param formula formula for the linear model.
#' @param mc_cores Number of cores for parallelization.
#' @param bycols Number of columns to group by.
#' @param ...  passed on to glm_dt. In particular, \code{bycols} can be passed for grouped lms. 
#' Set bycols = c('lon','lat') to estimate different regression coefficients for each gridpoint.
#' @param alongcols the coordinates specifying your repeated observations for the fit. Standard is 'year'. However, if you have a column named 
#' 'member' and you want to fit a regression for an ensemble system using the members as observations (thus allowing for multiple observations per year),
#' you need to pass alongcols = c('year','member') here. This is necessary to ensure that the predictions are sorted correctly into the data.table.
#'  
#' @return a data table containing coefficients and predictions.
#'  
#' @author Claudio

loyo_dt_lm = function(dt,formula,mc_cores = 1,bycols,...,alongcols = 'year')
{ # bycols: a new glm is fitted for each level of these variables
  
  if(length(dt[,.N,by = bycols][,unique(N)]) > 1)
  {
    stop("you have different numbers of observations for different grouping levels. That's currently not supported.")
  }
  
  years = unique(dt[,year])
  
  response = all.vars(formula)[1]
  
  fit_by_year = function(yy)
  {
    
    print(yy)
    
    dt_yy = glm_dt(dt[year != yy],formula,bycols = bycols,...,alongcols = alongcols,
                    coefficients_only = T
                    )
    
    # when you have multiple entries in alongcols, the returned coefficient vector needs to be repeated accordingly, which can be done like this:
    if(!identical(alongcols,'year'))
    {
      dt_yy = merge(dt_yy,dt[year == yy],by = intersect(colnames(dt_yy),colnames(dt[year == yy])),
                    allow.cartesian = T)
    }
    
    # use this model to predict for the left out year
    covariates =  all.vars(formula)[2:length(all.vars(formula))]
    coef_colnames =  paste0(all.vars(formula)[2:length(all.vars(formula))],'_coef')
    
    sdcols = intersect(c('(Intercept)',coef_colnames),colnames(dt_yy))
    coef_dt = dt_yy[,.SD,.SDcols = sdcols]
    
    setkeyv(dt,key(dt_yy))
    covs = dt[year == yy,.SD,.SDcols = covariates]
    if('(Intercept)' %in% colnames(dt_yy)) 
    {
      covs[,intercept := 1]
      setcolorder(covs,'intercept')
    }
    
    predictions = rowSums(as.matrix(covs)*as.matrix(coef_dt))
    
    dt_yy[,prediction := predictions]
    
    dt_yy[,year := yy]
    
    if(!identical(alongcols,'year'))
    {
      dt_yy[,(coef_colnames) := NULL]
    }
    
    setcolorder(dt_yy,'year')
    return(dt_yy)
  }
    
  if(mc_cores == 1)
  {
    ret_dt = data.table()
    
    for(yy in years)
    {
      dt_temp = fit_by_year(yy)
      ret_dt = rbindlist(list(ret_dt,dt_temp))
    }
  }
  
  if(mc_cores > 1)
  {
    
    ret_dt = rbindlist(parallel::mclapply(X = years,
                                          FUN = fit_by_year,
                                          mc.cores = mc_cores))
    
  }
  
  return(ret_dt)
}






#####
# This is a bit work in progress but fits a glm by considering every gridpoint within a fixed distance in the fitting process.

glm_dt_spatial = function(dt,formula,distance = 800,...,byvars = NULL,mc_cores = 8,verbose = F)
{
  
  prep_dt = copy(unique(dt[,.(lon,lat)]))
  prep_dt = as.data.table(expand.grid(lon = unique(prep_dt[,lon]),lat = unique(prep_dt[,lat]),
                                      lon0 = unique(prep_dt[,lon]),lat0 = unique(prep_dt[,lat])))
  
  prep_dt[,dist := geosphere::distHaversine(p1 = matrix(c(lon,lat),ncol = 2),
                                            p2 = matrix(c(lon0,lat0),ncol = 2))/1000]
  
  setkey(prep_dt,lon,lat)
  
  prep_dt[,indicator := dist <= distance]
  
  
  allvars = all.vars(formula)
  response = allvars[1]
  
  
  setkey(dt,lon,lat)
  
  coords = unique(dt[!is.na(get(response)),.(lon,lat)])
  
  if(mc_cores == 1)
  {
    
    dt_glm = data.table()
    for(ind in 1:coords[,.N])
    {
      if(verbose & ind %% 50 == 0) print(paste0(ind,'/',coords[,.N]))
      
      llon = coords[ind,lon]
      llat = coords[ind,lat]
      
      rel_locs = prep_dt[lon0 == llon & lat0 == llat & indicator][,.(lon,lat)]  
      
      dt_sub = merge(dt,rel_locs)
      
      temp_dt = glm_dt(dt_sub,formula,...,byvars = byvars)[lon == llon & lat == llat]
      dt_glm = rbindlist(list(dt_glm,temp_dt))
    }
  }
  
  if(mc_cores > 1)
  {
    fit_glm_parallel = function(ind)
    {
      if(verbose & ind %% 50 == 0) print(paste0(ind,'/',coords[,.N]))
      
      llon = coords[ind,lon]
      llat = coords[ind,lat]
      
      rel_locs = prep_dt[lon0 == llon & lat0 == llat & indicator][,.(lon,lat)]  
      
      dt_sub = merge(dt,rel_locs)
      
      temp_dt = glm_dt(dt_sub,formula,...,byvars = byvars)[lon == llon & lat == llat]
      return(temp_dt)
      
    }
    
    dt_glm = parallel::mclapply(1:coords[,.N],fit_glm_parallel,mc.cores = mc_cores)
    dt_glm = rbindlist(dt_glm)
  }
  
  return(dt_glm)
}



#' @description A function for lagging monthly variables, useful for e.g. regressing Julys temperature on Junes soil moisture.
#'
#' Names of columns (other than lon,lat,year,month) of dt that you want to keep, but shouldn't be lagged, should be included in \code{vars}. Just put a 0 at the corresponding position in lags.
#'  @param dt a data table with columns year, month, lon, lat and \code{vars}
#'  @param vars vector with names of columns of dt that should be lagged
#'  @param  lags integer vector (same length as vars) containing the lags for the variables, can contain zeros and negative numbers
#'  
#'  @return a data table containing columns year, month, lon, lat and \code{vars}. The values in \code{vars} are lagged by \code{lags} months.
#'  Only the time coordinates are kept, for which the lagged time coordinates are still in dt for all lags contained in \code{lags}
#'  @author Claudio
#'  @export

lag_dt = function(dt,vars,lags)
{
  setkey(dt,year,month)
  
  time_coor = unique(dt[,.(year,month)])    
  time_coor[,ind:= 1:.N]
  
  if(length(vars) != length(lags)) stop('the number of columns and provided lags need to match')
  
  for(ii in 1:length(vars))
  {
    temp_dt = lagged_time(time_coor,lags[ii])[,paste0('new_ind_',ii) := 1:.N]
    time_coor =  merge(time_coor,temp_dt,all.x = TRUE)
    
    assign(paste0('time_coor_',ii),temp_dt[,.(year,month)])
    setkey(get(paste0('time_coor_',ii)),year,month)
  }
  
  #find row indices that are present for all lags
  red_indices = Reduce(intersect, time_coor[,3:ncol(time_coor)])
  time_coor_new = time_coor[red_indices,.(year,month)]
  
  # initialize return_dt
  
  return_dt = dt[time_coor_new,.(year,month,lon,lat)]
  
  for(ii in 1:length(vars))
  {
    tc_temp = get(paste0('time_coor_',ii))[red_indices,]
    
    dt_temp = dt[tc_temp,.(year,month,lon,lat,get(vars[ii]))]
    
    return_dt[,(vars[ii]) := dt_temp[[5]]]
  }
  
  return(return_dt)
}

#####################


#' @description Auxiliary function helping to lag variables by a fixed amount of months, without getting confused at break of years.
#'
#'  @param time_coor a data table containing years and months.
#'  @param ll An integer for the lag, can be 0 or negative
#'  
#'  @return a data table contining the lagged time coordinates
#'  @author Claudio

lagged_time = function(time_coor,ll)
{ # given a dt with month and year, it returns a data table with year and month which contains the original months lagged by  ll months
  month = time_coor[,month]
  year = time_coor[,year]
  new_month = month - ll
  new_year = year
  new_year[new_month < 1] = new_year[new_month < 1] - 1
  new_month[new_month < 1] = new_month[new_month < 1] + 12
  
  return(data.table(year = new_year,month = new_month))
}


#' @description Leave-one-out cross validation for glms using the Brier score.
#' 
#' Currently assumes logistic regression.
#' @param dt a data table with columns year, month, lon, lat 
#' @param formula formula for glm, all covariates need to be column names in dt.
#' @param loo column name along which one is left out. Currently this might only work properly if this is 'year'.
#' @param byout Columns by which the output is grouped. Over all other columns the Brier score is averaged. Useful for permutation tests.
#' @param mc_cores Number of cores for parallelization.
#' @param ...  passed on to glm_dt. In particular, \code{family} should be correctly specified and \code{byvars} can be passed for grouped glms.
#'  
#' @return a data table containing Brier scores.
#'  
#' @author Claudio

loocv_glm_BS = function(dt,formula,loo = 'year',byout = NULL,mc_cores = 1,spatial = TRUE,...)
{# byvars: a new glm is fitted for each level of these variables
  #loo: currently supports only a single column name, mostly will be 'year'
  # byout: should the error be grouped by variables?
  levels = unique(dt[,.SD,.SDcols = loo])
  
  
  if(spatial)
  {
    fitting_function = function(dt,formula,...)
    {
      return(glm_dt_spatial(dt,formula,...,mc_cores = 1)) #set mc_cores to 1 to avoid multiple branching
    } 
  }else {
    fitting_function = function(dt,formula,...)
      {
      return(glm_dt(dt,formula,...))
      }
    
  }
  
  response = all.vars(formula)[1]
  if(mc_cores == 1)
  {
    ret_dt = data.table()
    
    for(ind in 1:levels[,.N])
    {
      ll = as.numeric(levels[ind,])
      
      print(ll)
      
      dt_glm = fitting_function(dt[get(loo) != ll],formula,...)
      
      # use this model to predict for the left out year
      covariates =  all.vars(formula)[2:length(all.vars(formula))]
      coef_colnames =  paste0(all.vars(formula)[2:length(all.vars(formula))],'_coef')
      coef_colnames = coef_colnames[which(coef_colnames %in% colnames(dt_glm))]
      
      coef_dt = dt_glm[year == min(year),.SD,.SDcols = coef_colnames]
      Intercept = dt_glm[year == min(year),.SD,.SDcols = '(Intercept)']
      
      covs = dt[year == ll,.SD,.SDcols = covariates]
      
      predictions = Intercept + rowSums(as.matrix(covs)*as.matrix(coef_dt))
      # inverse link function:
      predictions = 1/(1+exp(-predictions))
      
      err_dt = copy(dt[year == ll])[,fitted := predictions]
      err_dt = err_dt[,.SD,.SDcols = c('year','month','lon','lat',response,'fitted')]
      err_dt = err_dt[, mean((get(response) - fitted)^2),by = byout]
      
      setnames(err_dt,ncol(err_dt),'BS')
      
      err_dt[,paste0('lo_',loo) := ll]
      
      ret_dt = rbindlist(list(ret_dt,err_dt))
    }
  }
  
  if(mc_cores > 1)
  {
    dummy_f_par = function(ind)
    {
      ll = as.numeric(levels[ind,])
      
      print(ll)
      
      #dt_glm = fitting_function(dt[get(loo) != ll],formula,...)
      dt_glm = fitting_function(dt[get(loo) != ll],formula,...)
      
      # use this model to predict for the left out year
      covariates =  all.vars(formula)[2:length(all.vars(formula))]
      coef_colnames =  paste0(all.vars(formula)[2:length(all.vars(formula))],'_coef')
      coef_colnames = coef_colnames[which(coef_colnames %in% colnames(dt_glm))]
      
      coef_dt = dt_glm[year == min(year),.SD,.SDcols = coef_colnames]
      Intercept = dt_glm[year == min(year),.SD,.SDcols = '(Intercept)']
      
      covs = dt[year == ll,.SD,.SDcols = covariates]
      
      predictions = Intercept + rowSums(as.matrix(covs)*as.matrix(coef_dt))
      # inverse link function:
      predictions = 1/(1+exp(-predictions))
      
      err_dt = copy(dt[year == ll])[,fitted := predictions]
      err_dt = err_dt[,.SD,.SDcols = c('year','month','lon','lat',response,'fitted')]
      err_dt = err_dt[, mean((get(response) - fitted)^2),by = byout]
      
      setnames(err_dt,ncol(err_dt),'BS')
      
      err_dt[,paste0('lo_',loo) := ll]
      
      return(err_dt)
    }
    
    ret_dt = rbindlist(parallel::mclapply(X = 1:levels[,.N],
                                          FUN = dummy_f_par,
                                          mc.cores = mc_cores))
    
  }
  
  return(ret_dt)
}


#' @description Leave-one-year-out application of the regression specified in \code{spatial_eof_regression_dt}
#'
#'  @param dt Input data table, needs to contain year,lon,lat and all variables from formula
#'  @param formula formula for the linear regression
#'  @param nv How many PCs (EOFs) do you want to use?
#'  @param mc_cores How many cores to use in parallelization.
#'  
#'  @return A data table with columns lon,lat,year, prediction.
#'    
#'  @author Claudio
#'  @export


loyo_eof_regression = function(dt,formula,nv = 1,mc_cores = 1)
{
  fit_by_year = function(yy)
  {
    dt_temp = dt[year != yy]
    eof_reg = spatial_eof_regression_dt(dt_temp,formula,nv = nv)
    prediction = predict_eof_regression(eof_reg,newdata = dt[year == yy])[,.(lon,lat,year,value)]
    setnames(prediction,'value','prediction')
    
    return(prediction)
  }
  
  if(mc_cores == 1)
  {
    ret_dt = data.table()
    
    for(yy in unique(dt[,year]))
    {
      print(yy)
      ret_dt = rbindlist(list(ret_dt,fit_by_year(yy)))
    }
  } else {
    ret_dt = rbindlist(parallel::mclapply(unique(dt[,year]),FUN = fit_by_year,mc.cores = mc_cores))
  }
  
  return(ret_dt)
}



#' @description Predict values with the output from \code{spatial_eof_regression_dt}
#'
#'  @param eof_output Output of the \code{spatial_eof_regression_dt}-function.
#'  @param newdata provide this if new covariate-data should be used in prediction.#'  
#'  
#'  @return A data table. The prediction column is named 'value'
#'    
#'  @author Claudio
#'  @export


predict_eof_regression = function(eof_output,newdata = NULL)
{
  nv = ncol(eof_output$eofs) - 2 # the eof-dt contains the eofs and columns lon and lat.
  
  if(nv == 1)
  {
    covariates = names(eof_output$model$coefficients)
  }
  if(nv > 1)
  {
    covariates = rownames(eof_output$model$coefficients)
  }
  
  intercept_present = '(Intercept)' %in% covariates
  
  if(!is.null(newdata))
  {
    newdata = unique(newdata[,.SD,.SDcols = intersect(c('year',covariates),colnames(newdata))])
    
    if(intercept_present )
    {
      newdata[,'(Intercept)':= 1]
    }
    
    subnewdata = newdata[,.SD,.SDcols = covariates]
    
    
    fl_pred = as.matrix(subnewdata) %*% as.matrix(eof_output$model$coefficients)
    
  } else {
    fl_pred = as.matrix(eof_output$model$fitted.values)
    newdata = as.data.table(eof_output$model$model[[2]])
  }
  
  
  eofs = eof_output$eofs[,3:ncol(eof_output$eofs)]
  
  pred = as.matrix(eofs) %*% t(fl_pred)
  pred = as.data.table(pred)
  pred = data.table(pred,eof_output$eofs[,.(lon,lat)])
  pred = melt(pred,id.vars = c('lon','lat'))
  
  # attach predictor 
  index_vec = pred[,match(variable,unique(pred[,variable]))]
  pred = data.table(pred, newdata[index_vec,])
  pred[,variable := NULL]
  
  return(pred)
}


#' @description Leave-one-out cross validation for a climatological model using the Brier score.
#' 
#' Works only for 0-1 events.
#' 
#' @param dt a data table with columns year, month, lon, lat, and \code{response}. 
#' @param response column name for which climatology should be computed. Should be logical.
#' @param loo column name along which one is left out. Currently this should be 'year'.
#' @param groupby Columns to group by. 
#' @param byout Columns by which the output is grouped. Over all other columns the Brier score is averaged. Useful for permutation tests.
#' @param mc_cores Number of cores for parallelization.
#'  
#' @return a data table containing the Brier scores.
#'  
#' @author Claudio


loocv_clim_BS = function(dt, response, loo = 'year',groupby = c('lon','lat','month'),byout = NULL,mc_cores = 1)
{ # fits climatological model and performs loocv with the Brier score
  #loo: currently supports only a single column name, mostly will be 'year'
  # byout: should the error be grouped by variables?
  levels = unique(dt[,.SD,.SDcols = loo])
  setkey(dt,year,month,lon,lat)
  
  if(mc_cores == 1)
  {
    ret_dt = data.table()
    
    for(ind in 1:levels[,.N])
    {
      ll = as.numeric(levels[ind,])
      
      print(ll)
      
      dt_clim = dt[get(loo) != ll,][,.(year = year,month = month,lon = lon,lat = lat,clim = mean(get(response))),by = groupby]
      setkey(dt_clim,year,month,lon,lat)
      
      err_dt = dt_clim[year == min(year)]
      err_dt[,year:=ll]
      err_dt[,(response) := dt[get(loo) == ll,get(response)]]
      err_dt2 = err_dt[,mean((get(response) - clim)^2),by = byout]
      setnames(err_dt2,ncol(err_dt2),'BS')
      err_dt2[,lo_year := ll]
      ret_dt = rbindlist(list(ret_dt,err_dt))
    }
  }
  
  if(mc_cores > 1)
  {
    dummy_f_par = function(ind)
    {
      ll = as.numeric(levels[ind,])
      
      print(ll)
      
      dt_clim = dt[get(loo) != ll,][,.(year = year,month = month,lon = lon,lat = lat,clim = mean(get(response))),by = groupby]
      setkey(dt_clim,year,month,lon,lat)
      
      err_dt = dt_clim[year == min(year)]
      err_dt[,year:=ll]
      err_dt[,(response) := dt[get(loo) == ll,get(response)]]
      err_dt2 = err_dt[,mean((get(response) - clim)^2),by = byout]
      setnames(err_dt2,ncol(err_dt2),'BS')
      err_dt2[,lo_year := ll]
      
      
      return(err_dt2)
    }
    ret_dt = rbindlist(parallel::mclapply(X = 1:levels[,.N],
                                          FUN = dummy_f_par,
                                          mc.cores = mc_cores))
    
  }
  
  return(ret_dt)
}



#' @description Run an EOF regression on a data table. 
#'
#' This function runs linear regression on the factor loadings of spatial EOFs. 
#' Typical example (and where it was developed) is when you want to predict precipitation (spatially), and your covariate is univariate, like N3.4 index.
#' 
#'
#'  @param dt a data table with columns year, lon, lat, and all variables contained in formula
#'  @param formula formula for lm.
#'  @param nv how many principal components (EOFs) should be considered in the regression?
#'  
#'  
#'  @return a list of length 2: First entry is the model of the (multivariate) linear regression (after projecting into EOF space).
#'  The second page is a data table containing the EOFs.
#'    
#'  @author Claudio
#'  @export


spatial_eof_regression_dt = function(dt,formula,nv = 1)
{
  
  allvars = all.vars(formula)
  response = allvars[1]
  
  eofs = get_spatial_EOFs(dt,cn = response,bycols = c('lon','lat'),nv = nv)
  
  # get factor loadings:
  dt_fls = merge(dt,eofs,by = c('lon','lat'))
  dt_fls = dt_fls[,lapply(.SD,FUN = function(x){sum(x * get(response))}),.SDcols = paste0(response,'_eof_',1:nv),by ='year']
  dt_fls_red = as.matrix(dt_fls[,.SD,.SDcols = paste0(response,'_eof_',1:nv)])
  
  
  regression_dt = unique(dt[,.SD,.SDcols = c('year',allvars[-1])])
  setkey(regression_dt)
  regression_dt_red = as.matrix(regression_dt[,.SD,.SDcols = allvars[-1]])
  
  
  model = lm(dt_fls_red ~ regression_dt_red)
  
  #if either response or covariates are NOT a matrix, the coefficient names are weird, we now correct that:
  if(is.matrix(model$coefficients))
  {
    cns = colnames(model$coefficients)
    rns = rownames(model$coefficients)
    
    cns[cns == 'dt_fls_red'] = paste0(response,'_eof_1')
    rns[rns == 'regression_dt_red'] = allvars[-1]
    
    cutname = function(x){
      splits = unlist(strsplit(x,split = 'regression_dt_red'))
      if(length(splits) > 1)splits = splits[2]
      return(splits)
    }
    rns = unlist(lapply(rns,FUN = cutname))
    
    colnames(model$coefficients) = cns
    rownames(model$coefficients) = rns
  } else {
    ns = names(model$coefficients)
    ns[ns == 'regression_dt_red'] = allvars[-1]
    cutname = function(x){
      splits = unlist(strsplit(x,split = 'regression_dt_red'))
      if(length(splits) > 1)splits = splits[2]
      return(splits)
    }
    ns = unlist(lapply(ns,FUN = cutname))
    names(model$coefficients) = ns
  }
  
  return(list(model = model,eofs = eofs))
}


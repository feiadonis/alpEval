library(ggplot2)
library(icsUtil)


getICarray <- function(x, alphaNames, fwdNames, univName, method = "pearson"){
  periods = name2value(fwdNames)[['v']]
  tmp = array(0, dim = c(length(fwdNames), length(alphaNames)), dimnames = list(periods, alphaNames))
  for ( a in 1:length(alphaNames)){
    for ( ff in 1:length(alphaNames)){
      xx = x[,alphaNames[a]]
      yy = x[,fwdNames[a]]
      xx = xx[which(x[,univName]==1)]
      yy = yy[which(x[,univName]==1)]
      if (period[ff] > 0){ tmp[ff,a] = cor(xx,yy,method = method, use = "pairwise.complete.obs")}
      else { tmp[ff,a] = cor(xx,yy,method = method, use = "pairwise.complete.obs")}
    }
  }
  return(tmp)
}


calcICDateRange <- function( mdFwd, mdAlpha, mdUniv, method = "pearson"){
  md = panel.combine(list(mdFwd, mdAlpha, mdUniv))
  tnames = dimnames(md)[['T']]
  alphaNames = dimnames(mdAlpha)[['V']]
  fwdNames = dimnames(mdFwd)[['V']]
  univName = dimnames(mdUniv)[['V']]
  tuniv = dimnames(mdUniv)[['T']]
  for ( t in tnames){
    md[,,t,univName] = md[,,tuniv,univName]
  }
  md1 = aaply(md, c(2,3), getICarray, alphaNames, fwdNames, univName, method = method, .drop=FALSE)
  md1 = aperm(md1, c(3,1,2,4))
  names(dimnames(md1)) = c("K","D","T","V")
  return(md1)
}

calcICMean <- function(mdICDateRange, addZero = TRUE){
  md = adrop(aaply(mdICDateRange, c(1,4), mean, na.rm = TRUE, .drop = FALSE),3)
  # add zero to ic mean
  if( addZero ){
    c1 = md[which(as.integer(dimnames(md)[['K']]) < 0),,drop=FALSE]
    c2 = md[which(as.integer(dimnames(md)[['K']]) > 0),,drop=FALSE]
    mdICMean = rbind(c1,0,c2)
    names(dimnames(mdICMean)) = c("K","V")
    dimnames(mdICMean)[['K']] = c(dimnames(c1)[['K']],c("0"), dimnames(c2)[["K"]])
  }
  else{
    mdICMean = md
    names(dimnames(mdICMean)) = c("K","V")
  }
  return(mdICMean)
}

calcICYears <- function(mdICDateRange){
  dates = dimnames(mdICDateRange)[["D"]]
  minYear = format(as.Date(min(dates),"%Y%m%d"),"%Y")
  maxYear = format(as.Date(max(dates),"%Y%m%d"),"%Y")
  mlist = llply(c("wholePeriod",minYear:maxYear), function(m){
    if(m == "wholePeriod"){tmpIC = mdICDateRange}
    else{
      stDate = as.integer(paste(m,"0101",sep=""))
      edDate = as.integer(paste(m,"1231",sep=""))
      tmpIC = mdICDateRange[,which(as.integer(dates)>=stDate & as.integer(dates) <= edDate),,,drop=FALSE]
    }
    mdmean = calcICMean(tmpIC)
    tmpdata = array(0,dim = c(dim(mdmean)[1],1,dim(mdmean)[2]),dimnames = list(dimnames(mdmean)[["K"]],m,dimnames(mdmean)[['V']]))
    names(dimnames(tmpdata)) = c("K","D","V")
    tmpdata[,1,] = mdmean
    tmpdata
  })
  mdICYears = panel.combine(mlist)
}

# calc IC T space
calcICMean.T <- function(mdICDateRange, focus.period){
  mdIC = adrop(mdICDateRange[as.character(focus.period),,,,drop=FALSE],1)
  md = adrop(aaply(mdIC,c(2,3),mean,na.rm=TRUE,.drop=FALSE),3)
}

# calc monthly effect
calcMonthlyEffect <- function(mdICDateRange, focus.period){
  mdICDateRange = adrop(aaply(mdICDateRange,c(1,2,4),mean,na.rm=TRUE,.drop=FALSE),4)
  mdIC = adrop(mdICDateRange[as.character(focus.period),,,drop=FALSE],1)
  dates = dimnames(mdIC)[['D']]
  dtmths = (as.integer(dates)%/%100)%%100
  umonths = unique(dtmths)
  umonths = umonths[order(umonths)]
  mlist = llply(umonths, function(m){
    a = aaply(mdIC, 2, function(x){mean(x[dtmths==m],na.rm=TRUE)},.drop=FALSE)
    dimnames(a)[[2]]=m
    names(dimnames(a)) = c("V","D")
    a = aperm(a,c(2,1))
  })
  meffect = panel.combine(mlist)
}

# calc ic IR
calcIR <- function(mdICDateRange){
  period = as.integer(dimnames(mdICDateRange)[['K']])
  mdICDateRange = mdICDateRange[which(period >= 0),,,,drop=FALSE]
  mdMean = calcICMean(mdICDateRange)
  tmp = adrop(aaply(mdICDateRange,c(1,4),sd,na.rm=TRUE,.drop = FALSE),3)
  mdStd = rbind(1,tmp)
  names(dimnames(mdStd)) = c("K","V")
  dimnames(mdStd)[["K"]] = c(c("0"), dimnames(tmp)[["K"]])
  mdIR = mdMean / mdStd
  periods = as.integer(dimnames(mdIR)[["K"]])
  for( p in 2:length(periods)){
    mdIR[p,] = mdIR[p,]*sqrt(242)/sqrt(periods[p])
  }
  return(mdIR)
}

calcIRYears <- function(mdICDateRange){
  dates = dimnames(mdICDateRange)[["D"]]
  minYear = format(as.Date(min(dates),"%Y%m%d"),"%Y")
  maxYear = format(as.Date(max(dates),"%Y%m%d"),"%Y")
  mlist = llply(c("wholePeriod",minYear:maxYear),function(m){
    if(m == "wholePeriod"){tmpIC = mdICDateRange}
    else{
      stDate = as.integer(paste(m,"0101",sep=""))
      edDate = as.integer(paste(m,"1231",sep=""))
      tmpIC = mdICDateRange[,which(as.integer(dates)>=stDate & as.integer(dates)<= edDate),,,drop=FALSE]
    }
    mdIR = calcIR(tmpIC)
    tmpdata = array(0,dim = c(dim(mdIR)[1],1,dim(mdIR)[2]),dimnames = list(dimnames(mdIR)[["K"]],m,dimnames(mdIR)[["V"]]))
    names(dimnames(tmpdata)) = c("K","D","V")
    tmpdata[,1,] = mdIR
    tmpdata
  })
  mdIRYears = panel.combine(mlist)
}

clacRetWithAlphaWeight <- function(x,alphaNames,fwdNames,univName){
  x1 = x[which(x[,univName]==1),c(alphaNames,fwdNames)]
  tmpdata = array(0,dim = c(length(fwdNames),length(alphaNames)),dimnames = list(name2value(fwdNames)[['v']],alphaNames))
  for ( i in 1:length(fwdNames)){
    for( j in 1:length(alphaNames)){
      tmpdata[i,j] = sum(x1[,alphaNames[j]]*x1[,fwdNames[i]],na.rm=TRUE)
    }
  }
  tmpdata
}

#calc alpha
calcAlpha <- function(mdFwd, mdAlpha, mdUniv, method = "pearson", quickTest = FALSE){



}


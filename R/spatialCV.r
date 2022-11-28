
#' @title Split data for k-fold spatially stratified cross validation
#'
#' @description Splitting tool for cross-validation
#'
#' @details
#' See Examples.
#'
#' @param sp.pts a SpatialPoints or SpatialPointsDataFrame object
#' @param nfolds number of desired output folds. Default value of -1 makes a reasonable guess based on sample size.
#' @param nsubclusters intermediate number of clusters randomly split into nfolds. Default value of -1 makes a reasonable guess based on sample size. If you specify this manually, it should be an integer multiple of nfolds.
#' @param diagnostic.plot plot the points colored by fold
# @keywords
#' @export
#'
# @examples
#'
#'
#' @return Returns a data dataframe by adding a column (`fold') to the existing dataframe that indicates which fold each sample is assigned to.
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
#' @seealso plotSDMfolds
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

# TODO
# consider making this a list, or adding the option to ouptut as a list
# env for the plot currently grabbed rfom environment
spatialStratify=function(sp.pts,
												 stats,
												 nfolds=-1,
												 nsubclusters=-1,
												 diagnostic.plot=FALSE,
												 verbose=5){

	out=try({

	species=stats$species[1]
	if(verbose>2) print(paste0(species,': ','     stratifying occurrences '))
	focal.pres=subset(sp.pts,sp==species)
	other.sp.pres=subset(sp.pts,!sp==species)
	#-- don't bother with cross validation if there are too few points
	pres.coords=sp::coordinates(focal.pres)
	other.sp.pres.coords=sp::coordinates(other.sp.pres)

	# charlie
	# folds.tmp=kmeans(pres.coords,2)
	# folds.tmp1=flexclust::as.kcca(folds.tmp, data=pres.coords)
	# focal.pres$folds=folds.tmp$cluster
	
	n=nrow(pres.coords)
	if (n<=5) {
		folds=rep(1, nrow(pres.coords))
		focal.pres$folds=folds
		stop('Get some more data if you want to do spatial CV. All samples assigned to fold 1.')
	}

	if(nfolds==-1) {
		if (n>5 & n<=15) {
			folds.tmp=kmeans(pres.coords,5)
			folds.tmp1=flexclust::as.kcca(folds.tmp, data=pres.coords)
			focal.pres$folds=folds.tmp$cluster
			if(length(other.sp.pres)>0){
				other.sp.pres$folds=predict(folds.tmp1,other.sp.pres.coords)
			}
		}
		if (n>15 & n<=30) {
			folds.tmp=kmeans(pres.coords,10)
			folds.tmp1=flexclust::as.kcca(folds.tmp, data=pres.coords)
			folds=folds2=folds.tmp$clust
			combine.folds=matrix(sample(1:nrow(folds.tmp$centers), 10,replace=F),ncol=2)
			for(ii in 1:nrow(combine.folds)){ folds[folds2%in%c(combine.folds[ii,])]=ii }
			focal.pres$folds=folds
			if(length(other.sp.pres)>0){
				others.folds=others.folds2=predict(folds.tmp1,other.sp.pres.coords)
				for(ii in 1:nrow(combine.folds)){ others.folds[others.folds2%in%c(combine.folds[ii,])]=ii }
				other.sp.pres$folds=others.folds
			}
		}
		if (n>30 & n<=45) {
# 			folds.tmp=kmeans(pres.coords,15)
# 			folds=folds2=folds.tmp$clust
# 			combine.folds=matrix(sample(1:nrow(folds.tmp$centers), 15,replace=F),ncol=3)
# 			for(ii in 1:nrow(combine.folds)){ folds[folds2%in%c(combine.folds[ii,])]=ii
			folds.tmp=kmeans(pres.coords,15)
			folds.tmp1=flexclust::as.kcca(folds.tmp, data=pres.coords)
			folds=folds2=folds.tmp$clust
			combine.folds=matrix(sample(1:nrow(folds.tmp$centers), 15,replace=F),ncol=3)
			for(ii in 1:nrow(combine.folds)){ folds[folds2%in%c(combine.folds[ii,])]=ii }
			focal.pres$folds=folds
			if(length(other.sp.pres)>0){
				others.folds=others.folds2=predict(folds.tmp1,other.sp.pres.coords)
				for(ii in 1:nrow(combine.folds)){ others.folds[others.folds2%in%c(combine.folds[ii,])]=ii }
				other.sp.pres$folds=others.folds
			}
		}
		if (n>45 & n<=60) {
			folds.tmp=kmeans(pres.coords,20)
			folds.tmp1=flexclust::as.kcca(folds.tmp, data=pres.coords)
			folds=folds2=folds.tmp$clust
			combine.folds=matrix(sample(1:nrow(folds.tmp$centers), 20,replace=F),ncol=4)
			for(ii in 1:nrow(combine.folds)){ folds[folds2%in%c(combine.folds[ii,])]=ii }
			focal.pres$folds=folds
			if(length(other.sp.pres)>0){
				others.folds=others.folds2=predict(folds.tmp1,other.sp.pres.coords)
				for(ii in 1:nrow(combine.folds)){ others.folds[others.folds2%in%c(combine.folds[ii,])]=ii }
				other.sp.pres$folds=others.folds
			}
		}

		if (n>60) {
			# folds.tmp=kmeans(pres.coords,25)
# 			folds=folds2=folds.tmp$clust
# 			combine.folds=matrix(sample(1:nrow(folds.tmp$centers), 25,replace=F),ncol=5)
# 			for(ii in 1:nrow(combine.folds)){ folds[folds2%in%c(combine.folds[ii,])]=ii
			folds.tmp=kmeans(pres.coords,25)
			folds.tmp1=flexclust::as.kcca(folds.tmp, data=pres.coords)
			folds=folds2=folds.tmp$clust
			combine.folds=matrix(sample(1:nrow(folds.tmp$centers), 25,replace=F),ncol=5)
			for(ii in 1:nrow(combine.folds)){ folds[folds2%in%c(combine.folds[ii,])]=ii }
			focal.pres$folds=folds
			if(length(other.sp.pres)>0){
				others.folds=others.folds2=predict(folds.tmp1,other.sp.pres.coords)
				for(ii in 1:nrow(combine.folds)){ others.folds[others.folds2%in%c(combine.folds[ii,])]=ii }
				other.sp.pres$folds=others.folds
			}
		}
	} else { # end default splitting rules
		## 7/12/17 not sure what this was for, and haven;t fixed for other species
		# if(nsubclusters==-1) stop('you must specify nsubclusters')
		# folds.tmp=kmeans(pres.coords,nsubclusters)
		# folds=folds2=folds.tmp$clust
		# combine.folds=matrix(sample(1:nrow(folds.tmp$centers), nsubclusters,replace=F),ncol=nfolds)
		# 	for(ii in 1:nrow(combine.folds)){ folds[folds2%in%c(combine.folds[ii,])]=ii }
	  if(nsubclusters==-1) stop('you must specify nsubclusters')
	  folds.tmp=kmeans(pres.coords,nsubclusters)
	  folds=folds2=folds.tmp$clust
	  combine.folds=matrix(sample(1:nrow(folds.tmp$centers), nsubclusters,replace=F),nrow=nsubclusters,ncol=nfolds)
	  for(ii in 1:nrow(combine.folds)){ folds[folds2%in%c(combine.folds[ii,])]=ii }
	  focal.pres$folds=folds #es added to specify #folds manually
	}

	sp.pts=rbind(focal.pres,other.sp.pres)
	#sp.pts$folds=folds !!!!!!

	if(diagnostic.plot){
		plotSDMFolds(sp.pts,env)
	}
	return(sp.pts)
	})
	return(out)
}


# Testing
#points=read.csv('/Users/ctg/Dropbox/Projects/BIEN/bien-range/Test_Settings/for_local_setup/na_trees/Abies_balsamea.csv')
# coordinates(points)=c(2,1)
#
# points=read.csv('/Users/ctg/Dropbox/Projects/Dimensions/cfrsdm/data/Presences/Protea_Atlas/Proteaceae/Aulax_cancellata.csv')
# coordinates(points)=c('lon','lat')
#
# (points.cv=spatialStratify(points,-1,-1))
# (points.cv=spatialStratify(points,5,25))
# if(length(points.cv)==1) status['do.cv']=FALSE
# if(length(points.cv)>1) status['do.cv']=TRUE


#======================================================================
#======================================================================
#======================================================================

#' @title Summarize SDM output from k-fold cross validation for CFR
#'
#' @description Summarize SDM output from k-fold cross validation for CFR
#'
#' @details
#' See Examples.
#'
#' @param dirs a list of workflow directories
#' @param stats a dataframe of model inputs and outputs
#' @param writeSDPred logical; save the raster of SD over folds?
#' @param transfer integer; referencing which env scenario to transfer to. These scenarios should be listed in \code{dirs$otherEnvPred}.
# @keywords
#' @export
#'
# @examples
#'
#'
#' @return NULL
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

summarizeKFold=function(dirs,
												stats,
												writeSDPred,
												transfer=NULL,
												algorithm,
												ensembleWeights=NULL,
												verbose=5){
	#  for testing
	#  writeSDPred=toSave$SDPred; transfer=NULL; algorithm=toDo$misc$algorithm

	out=try({
		species=basename(dirs$sp.pred.path)
		if(verbose>2) print(paste0(species,': ','     summarizing cv models - ',ifelse(is.null(transfer),'',transfer)))
		modelPath=dirs$sp.pred.path # both in and out path
		modelPathSD=dirs$sp.pred.path.sd
		# rules for names from other use cases. it'd be nice if i could just specify the directory
		if(!is.null(transfer)) {
			topDir=allBaseDirs$algorithmDirs[algorithm][[1]]
			modelPath=paste0(topDir$OtherEnvPred[[transfer]]$predictionDir, '/',species)
			#paste0(dirs$otherEnvPred[[transfer]]$maxentPredDir,'/',species)
			if(!file.exists(modelPath)) dir.create(modelPath)
			modelPathSD=paste0(topDir$OtherEnvPred[[transfer]]$sdDir, '/',species)
			#paste0(dirs$otherEnvPred[[transfer]]$maxentSDDir,'/',species)
			if(!file.exists(modelPathSD)) dir.create(modelPathSD)
		}

		all.files=list.files(modelPath,pattern="fold",full.names=T)

		if(length(all.files)==0) {
			cat('no cv to summarize')
			return(list(stats=stats))
		}
			#-- to get models in the same order as stats
		psf=lapply(stats$modelNames,function(x) {all.files[grep(x,all.files)]})
		not.empty=unlist(lapply(psf,function(x) !length(x)==0))

		if(length(psf)>0){
				#-- read in models
			ps=lapply(psf,function(x) if(!length(x)==0) unrtrans(stack(x)))
				#-- get mean and sd of spatial predictions
			# hack solution
			tmp.stats=stats[1:length(stats$modelNames),]
			meanModelPath=paste0(modelPath,'/',species,'__mean__', stats$modelNames,'.tif')
			sdModelPath=paste0(modelPathSD,'/',species,'__sd__', stats$modelNames,'.tif')

			for(l in which(not.empty)){
				#keep=grep(stats$modelNames[l],names(ps))
				sp.mean=rtrans(raster::calc(ps[[l]],mean))
				writeRaster(sp.mean,filename = meanModelPath[l], options = c("COMPRESS=LZW", "PREDICTOR=2"), datatype = "INT2S", overwrite = TRUE)
				bossMaps::rmRaster(sp.mean)
				rm(sp.mean)
				if(writeSDPred){
					sp.sd=rtrans(raster::calc(ps[[l]],sd))
					writeRaster(sp.sd,na.rm=T,filename = sdModelPath[l], options = c("COMPRESS=LZW", "PREDICTOR=2"), datatype = "INT2S", overwrite = TRUE)
				}
				#bossMaps::rmRaster(sp.sd)
				rm(sp.sd)
				gc()
			} # end loop ove modelNames
			#lapply(ps,bossMaps::rmRaster)
			rm(ps)

			stats$meanModelPath=meanModelPath
			stats$sdModelPath=sdModelPath
			stats$meanModelPath[!not.empty]=stats$sdModelPath[!not.empty]=NA
		} else {
			stats$meanModelPath=NA
			stats$sdModelPath=NA
		}# end if no cv to model

		# toss files if no model fit
		stats=stats[!is.na(stats$meanModelPath),]

		return(stats)
	})
	return(out)
}



#======================================================================
#======================================================================
#======================================================================

#' @title Summarize SDM output from k-fold cross validation for CFR
#'
#' @description Summarize SDM output from k-fold cross validation for CFR
#'
#' @details
#' See Examples.
#'
#' @param dirs a list of workflow directories
#' @param stats a dataframe of model inputs and outputs
#' @param writeSDPred logical; save the raster of SD over folds?
#' @param transfer integer; referencing which env scenario to transfer to. These scenarios should be listed in \code{dirs$otherEnvPred}.
# @keywords
#' @export
#'
# @examples
#'
#'
#' @return NULL
#' @author Erica Stuber <efstuber@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

summarizeKFoldFast=function(dirs,
                        stats,
                        writeSDPred,
                        transfer=NULL,
                        algorithm,
                        ensembleWeights=NULL,
                        verbose=5){

  out=try({
    species=basename(dirs$sp.pred.path)
    if(verbose>2) print(paste0(species,': ','     summarizing cv models '))
    modelPath=dirs$sp.pred.path # both in and out path
    modelPathSD=dirs$sp.pred.path.sd
    # rules for names from other use cases. it'd be nice if i could just specify the directory
    if(!is.null(transfer)) {
      topDir=allBaseDirs$algorithmDirs[algorithm][[1]]
      modelPath=paste0(topDir$OtherEnvPred[[transfer]]$predictionDir, '/',species)
      #paste0(dirs$otherEnvPred[[transfer]]$maxentPredDir,'/',species)
      if(!file.exists(modelPath)) dir.create(modelPath)
      modelPathSD=paste0(topDir$OtherEnvPred[[transfer]]$sdDir, '/',species)
      #paste0(dirs$otherEnvPred[[transfer]]$maxentSDDir,'/',species)
      if(!file.exists(modelPathSD)) dir.create(modelPathSD)
    }

    all.files=list.files(modelPath,pattern="fold",full.names=T)

    if(length(all.files)==0) {
      cat('no cv to summarize')
      return(list(stats=stats))
    }
    #-- to get models in the same order as stats
    psf=lapply(stats$modelNames,function(x) {
      #string=strsplit(basename(x),'full_')[[1]][2]
      all.files[grep(x,all.files)]
    })
    not.empty=unlist(lapply(psf,function(x) !length(x)==0))
    #psf=psf[not.empty]

    if(length(psf)>0){
      #-- read in models
      ps=lapply(psf,function(x) if(!length(x)==0) unrtrans(stack(x)))

      #-- get mean and sd of spatial predictions
      # hack solution
      tmp.stats=stats[1:length(stats$modelNames),]
      meanModelPath=paste0(modelPath,'/',species,'_mean_', stats$modelNames,'.tif')
      sdModelPath=paste0(modelPathSD,'/',species,'_sd_', stats$modelNames,'.tif')

      for(l in which(not.empty)){
        #keep=grep(stats$modelNames[l],names(ps))
        sp.mean=rtrans(Mean_SD_raster(ps[[l]])[[1]]) #list is mean, SD
        #sp.mean=rtrans(raster::calc(ps[[l]],mean))
        writeRaster(sp.mean,filename = meanModelPath[l], options = c("COMPRESS=LZW", "PREDICTOR=2"), datatype = "INT2S", overwrite = TRUE);gc()
        if(writeSDPred){
          sp.sd=rtrans(Mean_SD_raster(ps[[l]],writeSDPred=T)[[2]]) #list is mean, SD
          #sp.sd=rtrans(raster::calc(ps[[l]],sd))
          writeRaster(sp.sd,na.rm=T,filename = sdModelPath[l], options = c("COMPRESS=LZW", "PREDICTOR=2"), datatype = "INT2S", overwrite = TRUE);gc()
        }
      } # end loop ove modelNames
      stats$meanModelPath=meanModelPath
      stats$sdModelPath=sdModelPath
      stats$meanModelPath[!not.empty]=stats$sdModelPath[!not.empty]=NA
    } else {
      stats$meanModelPath=NA
      stats$sdModelPath=NA
    }# end if no cv to model

    # toss files if no model fit
    stats=stats[!is.na(stats$meanModelPath),]

    return(stats)
  })
  return(out)
}




#======================================================================
#======================================================================
#======================================================================

#' @title Speedy Mean/SD calculation for raster stack
#'
#' @description Speedy Mean/SD calculation for raster stack
#'
#' @details
#' See Examples.
#'
#' @param writeSDPred logical; save the raster of SD over folds?
# @keywords
#' @export
#'
# @examples
#'
#'
#' @return NULL
#' @author Lifted from Samuel Bosch see: <www.samuelbosch.com/2018/02/speeding-up-mean-and-standard-deviation.html>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

Mean_SD_raster <- function(x, writeSDPred=F) {
  s0 <- nlayers(x)
  s1 <- raster(x, layer=1)
  s2 <- s1^2
  for(ri in 2:s0) {
    r <- raster(x, layer=ri)
    s1 <- s1 + r
    s2 <- s2 + r^2
  }
  if(writeSDPred){
    list(mean=s1/s0, sd=sqrt((s0 * s2 - s1 * s1)/(s0 * (s0 - 1))))
  } else {
    list(mean=s1/s0)
  }

}

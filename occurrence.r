#' @title Get, clean and thin occurrence data
#'
#' @description In preparation for an SDMWorkflow
#'
#' @details
#' See Examples.
#'
#' @param speciesCSV file path. column names in the csv must in include 'lon' and 'lat'
#' @param env raster
#' @param doThin logical; thin the points?
#' @param thinCutoff numeric; minimum distance in km between points
# @keywords
#' @export
#'
# @examples
#'
#'
#' @return Returns a list with presence and background points
#' @author Cory Merow <cory.merow@@gmail.com>
#' @note Thinning is just a wrapper for \code{spThin}
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

# should this keep the extracted layers so i don't have to get them again later?
cleanOcc=function(speciesCSV,
									presSP=NULL,
									env,
									pointsProj=projection(env),
									doThin=FALSE,
									thinCutoff=20,
									doClean=TRUE,
									writeResult=FALSE,
									overwrite=FALSE,
									stats=NULL,
									dirs,
									maxSimilarSpecies=5,
									verbose=5){

	#  for testing
	#  presSP=NULL; doClean=FALSE; writeResult=FALSE; maxSimilarSpecies=toDo$maxSimilarSpecies; overwrite=FALSE	; thinCutoff=20	;pointsProj=projection(env) ;

	out=try({

		species=tools::file_path_sans_ext(basename(speciesCSV))
		if(verbose>2) print(paste0(species,': ','     cleaning occurrences '))
		if(is.null(presSP)){
			if(file.exists(speciesCSV)){
				pres=utils::read.csv(speciesCSV)
				keep=c(grep('on',names(pres)),	grep('at',names(pres)))
				if(length(keep)==0) keep=c(grep('x',names(pres)),grep('y',names(pres)))
				if(length(keep)==0){ keep=c(2,1); warning('assuming presence files have columns lat and lon at 1 and 2')}
				coordinates(pres)=keep
				projection(pres)=pointsProj
			} else {
				pres=NULL
			}
		}
		if(!is.null(presSP)){ # for use in the sorting stage, where the spatial object was already made
			pres=presSP
		}
		# extra check, just in case fast sorting was used and cleaning wasn't performed
		pres=pres[!duplicated(coordinates(pres)),]

		#----------------------------------------------------------------
			#-- remove points outside the domain
		if( doClean & !is.null(pres) ){
			pres <- SpatialPoints(pres, proj4string=CRS(pointsProj))
			#pres=spTransform(pres,projection(env))
			keep=complete.cases(raster::extract(env,pres))
			pres=pres[keep,]

			if(doThin &&  length(pres)>thinCutoff){
				if(verbose>2) print(paste0(species,': ','  thinning presences'))
				#-- use regular thinning algorithm if <800 points
				if(length(pres)<800){
					#-- thin.algorithm wants degrees as units
					tmp.proj=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
					real.proj=projection(pres)
					tmp.pres=spTransform(pres,tmp.proj)
					tmp.pres2=data.frame(coordinates(tmp.pres), species=species)
					locs.thinned <- spThin::thin.algorithm(rec.df.orig = tmp.pres2, thin.par = 10, reps = 10)
					sample.sizes=which.max(unlist(lapply(locs.thinned,nrow)))[1]
					tmp1.pres=locs.thinned[[sample.sizes]][,c(3,1,2)]
					coordinates(tmp1.pres) =c(2,3)
					projection(tmp1.pres)=tmp.proj
					pres=spTransform(tmp1.pres,real.proj)
			  } else {
			  	#-- remove dups on a coarse grid if lots of points (thinning algorithm too slow)
			  	#coarse.grid=raster(paste0(dirs$miscDir,'coarse_grid_for_thinning.tif'))
			  	coarse.grid=raster(dirs$coarseGridFile)
			  	# Extract the unique cell numbers for the cells containing the points
					keep=!duplicated(cellFromXY(coarse.grid, pres))
					pres=pres[keep,]
			  }
			}
		} # end thin

		if(writeResult | overwrite){
		  if(verbose>0) print(paste0(species,': ','  writing presence csv'))
		  write.csv(pres,file=speciesCSV,row.names=FALSE)
		}

		#------------------------------------------------------------------
			#-- weights
		if( any(!is.null(stats$weightSet) & !stats$weightSet=='equalWeights') ){
		  if(verbose>0) print(paste0(species,': ','  calculating species weights'))
			d.w=utils::read.csv(dirs$otherDirs$spWeights,stringsAsFactors=F)
			d.w$sp=gsub(' ','_',d.w$sp)
			d.w$sp_neigh=gsub(' ','_',d.w$sp_neigh)
			d.w=subset(d.w,sp==species)
			d.w=d.w[,c(grep('sp_neigh',names(d.w)),grep('weight',names(d.w)))]
			weight.options=unique(stats$weightSet)
			toss.equal=grep('equalWeights',weight.options)
			if(length(toss.equal)>0) weight.options=weight.options[-grep('equalWeights',weight.options)]
			for(ii in weight.options){
				# can grab X most similar species or set a thresold at similarity Y
				keep=tail(order(d.w[,ii]),maxSimilarSpecies)
				related.sp=gsub(' ','_',d.w$sp_neigh[keep])
				all.sp.points=utils::read.csv(paste0(dirs$miscDir, '/allSpeciesPoints.csv'), stringsAsFactors=F)
				related.sp.pts=subset(all.sp.points,(sp %in% related.sp))
				related.sp.pts1=merge(related.sp.pts,d.w,by.x='sp',by.y='sp_neigh')
		 		pres2=data.frame(sp=species,lon=coordinates(pres)[,1], lat=coordinates(pres)[,2])
				aa=plyr::rbind.fill(pres2,related.sp.pts1)
				aa[aa$sp==species,grep('weight',names(aa))]=1
				coordinates(aa)=c('lon','lat')
				projection(aa)=CRS(proj4string(pres))
			}
			pres=aa
			pres$equalWeights=rep(1,length(pres)) #placeholder for no weights
			if(verbose>3) print(paste0(species,': ','     ',nrow(related.sp.pts1),' point from similar taxa'))

		} else {
			pres$sp=rep(species,length(pres))
			pres$equalWeights=rep(1,length(pres)) #placeholder for no weights
		}

		#------------------------------------------------------------------
		#-- manually remove a fake point in many data sets
		  #-- complicated due to some bullshit about machine precision that means == doesn't work
		toss=  which(abs(coordinates(pres)[,1]-(-1607582)) < 1 & abs(coordinates(pres)[,2]-(2625635)) < 1)
		if(length(toss)>0) pres=pres[-toss,]

		return(list(pres=pres))
	})
	return(out)
}

#======================================================================
#======================================================================
#' @title Find outlying occurrence data in geographic space
#'
#' @description Spatial  outliers
#'
#' @details
#' See Examples.
#'
#' @param pres a `SpatialPointsDataFrame`
#' @param pvalSet numeric; p-value used in Grubb's test for outlier (see package `outliers`)
#' @param checkPairs logical; check for a single pair of outliers. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye.
# @keywords
#' @export
#'
# @examples
# myPres=read.csv(system.file('ext/SampleData/Sp3Occurrence_v3.csv',
# package='occProfileR'))[,c(2,3)]
# myPres=myPres[complete.cases(myPres),]
# coordinates(myPres)=c(1,2)
# spOut=findSpatialOutliers(pres=myPres,pval=1e-5)
#'
#' @return Returns a vector of the IDs of outlying points
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

findSpatialOutliers=function(myPres,
                             pvalSet=1e-5,
                             checkPairs=T){
  #  for testing
  #  myPres=pres; pvalSet=1e-5
  #  myPres=mcp.pres; pvalSet=toDo$presences$pvalSet; checkPairs=T

  pres.inliers=myPres
  sp.toss.coord=NULL
  pval=0
  #toss singles
  while(pval<pvalSet){
  	dists=presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
    gt=outliers::grubbs.test(dists)
    #gt=dixon.test(dists)
  	#cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
    (pval=gt$p.value)
    # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
    if(gt$p.value<pvalSet){
      toss=which.max(dists)
      # IDs in the original data frame
      sp.toss.coord=rbind(sp.toss.coord,coordinates(pres.inliers)[toss,])
      pres.inliers=pres.inliers[-toss,]
      dists=dists[-toss]
    }
  }
  # toss pairs
  if(checkPairs){
  	if(length(pres.inliers)<31){
			pval=0
			# By turning off this loop, I'm ensuring that you can only toss 1 pair of outliers. with the loop, it tends to find lots of supposed outliers very confidently, but by eye, it tends to omit clusters
			#while(pval<pvalSet){
				dists=presPercentile(pres.inliers, percent=NULL)[[1]]$dist.from.centroid
				gt=outliers::grubbs.test(dists,type=20)

				#gt=dixon.test(dists)
				#cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
				(pval=gt$p.value)
				# conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
				if(gt$p.value<pvalSet){
					toss=tail(order(dists),2)
					# IDs in the original data frame
					sp.toss.coord=rbind(sp.toss.coord, coordinates(pres.inliers)[toss,])
					pres.inliers=pres.inliers[-toss,]
				}
			#}
		}
  }

  if(!is.null(sp.toss.coord)){
  	coor=coordinates(myPres)
  	sp.toss.id= apply(sp.toss.coord,1,function(x) which(x[1]==coor[,1] & x[2]==coor[,2]))
  } else {sp.toss.id=NULL}
  print(paste0(length(sp.toss.id),' geographic outliers found'))
  sp.toss.id
}

#======================================================================
#======================================================================
#' @title Find outlying occurrence data in environmental space
#'
#' @description Environmental outliers
#'
#' @details
#' See Examples.
#'
#' @param myPres a `SpatialPointsDataFrame`
#' @param env a `RasterStack` of env layers. If NULL, it is assumed that `myPres` is a data.frame of the environmental values (columns) at presence locations (rows)
#' @param pvalSet numeric; p-value used in Grubb's test for outlier (see package `outliers`)
# @keywords
#' @export
#'
# @examples
# myPres=read.csv(system.file('ext/SampleData/Sp3Occurrence_v3.csv',
#  package='occProfileR'))[,c(2,3)]
#  myPres=myPres[complete.cases(myPres),]
#  coordinates(myPres)=c(1,2)
# myEnv=raster::stack(system.file('ext/AllEnv.tif',package='occProfileR'))
#  envOut=findEnvOutliers(myPres=myPres,env=myEnv,pval=1e-5)
#'
#' @return Returns a list of SpatialPointsDataFrames with (1) good presence points (2) spatial outliers and (3) environmental outliers.
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

findEnvOutliers=function(myPres,
                         myEnv=NULL,
                         pvalSet=1e-5,
                         checkPairs=F){
  #  for testing
  #  myPres=presDF; pvalSet=1e-5; checkPairs=F; myEnv=NULL
  #  myEnv=env
  #  myPres=myPres; env=myEnv; pvalSet=1e-5
  if(!is.null(myEnv)){ p.env=raster::extract(myEnv,myPres)
  } else {p.env=myPres}
	p.env=scale(p.env)
	# remove variables that are the same for all observations
	f=which(apply(p.env,2,function(x) !all(is.nan(x))))
	p.env=p.env[,f]
  pres.inliers=p.env
  row.id=apply( p.env, 1 , paste , collapse = "-" )
  env.toss.id=NULL
  pval=0
  while(pval<pvalSet){
    dists=apply(p.env,1,function(x) sqrt(sum((x)^2)) )
    gt=outliers::grubbs.test(dists)
    pval=gt$p.value
    # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
    if(gt$p.value<pvalSet){
      toss=which.max(dists)
      # IDs in the original data frame
      thisID=paste(p.env[toss,],collapse='-')
      env.toss.id=c(env.toss.id,which(row.id == thisID))
      p.env=p.env[-toss,]
    }
  }
  if(checkPairs) print('checkPairs not yet implemented')
  env.toss.id=unique(env.toss.id)
  print(paste0(length(env.toss.id),' environmental outliers found'))
  env.toss.id
}

#======================================================================
#======================================================================
#' @title Find outlying occurrence data
#'
#' @description Spatial or environmental outliers
#'
#' @details
#' See Examples.
#'
#' @param pres a SpatialPointsDataFrame
#' @param spOutliers logical; perform spatial outlier analysis
#' @param envOutliers logical; perform environmental outlier analysis
#' @param doPlot FALSE
#' @param env NULL
#' @param plotFile NULL
#' @param species NULL
# @keywords
#' @export
#'
# @examples
# myPres=read.csv(system.file('ext/SampleData/Sp3Occurrence_v3.csv',
#                           package='occProfileR'))[,c(2,3)]
# myPres=myPres[complete.cases(myPres),]
# coordinates(myPres)=c(1,2)
# myEnv=raster::stack(system.file('ext/AllEnv.tif',package='occProfileR'))
# presOut=findOutlyingPoints(pres=myPres,
#                            spOutliers=T,
#                            envOutliers=T,
#                            doPlot=FALSE,
#                            env=NULL,
#                            plotFile=NULL,
#                            species=NULL,
#                            pval=1e-5)
#'
#' @return Returns a list of SpatialPointsDataFrames with (1) good presence points (2) spatial outliers and (3) environmental outliers.
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.


findOutlyingPoints=function(pres,
                            spOutliers=T,
                            envOutliers=T,
                            bestVar=NULL,
                            doPlot=FALSE,
                            env=NULL,
                            plotFile=NULL,
                            species=NULL,
                            pval=1e-5){

  #  for testing
  #  pres=myPres; env=myEnv;
  #  envOutliers=T; spOutliers=T; doPlot=T; plotFile=NULL; bestVar=best.var; pval=toDo$presences$pvalSet
  #  plotFile=paste0(dirs$sp.diag.path, '/OccurrenceOutliers.png')
  if(!is.null(species)){
    keep=pres$sp==species
    otherSp=pres[!keep,] # append this back on later
    if(length(otherSp)==0) otherSp=NULL; keep=1:length(pres)
  } else {otherSp=NULL; keep=1:length(pres)}
  mcp.pres=data.frame(coordinates(pres)[keep,])
  coordinates(mcp.pres)=c(1,2)

  if(spOutliers) { sp.toss.id=findSpatialOutliers(myPres=mcp.pres,pvalSet=pval)
	} else {sp.toss.id=NULL}
  if(envOutliers) {
  	presDF=pres[,bestVar]@data
  	env.toss.id=findEnvOutliers(myPres=presDF,pvalSet=pval)
	} else {env.toss.id=NULL}
  if(!is.null(c(sp.toss.id,env.toss.id))){
    spatialToss=pres[sp.toss.id,]
    envToss=pres[env.toss.id,]
    pres0=pres[-c(sp.toss.id,env.toss.id),]
  } else{ pres0=pres; spatialToss=NULL; envToss=NULL}

  if(!is.null(c(sp.toss.id,env.toss.id))){ # don't plot if there's no outliers
    if(doPlot){
      if(!is.null(plotFile)) {png(plotFile,h=1000,w=1000)}
      pe=extent(pres)
      ydif=pe[4]-pe[3] ; xdif=pe[2]-pe[1] # for plotting
      pe[1]=pe[1]-xdif; pe[2]=pe[2]+xdif; pe[3]=pe[3]-ydif; pe[4]=pe[4]+ydif
      env.bg=raster::crop(env[[1]],pe)
      plot(env.bg,col='grey50',legend=F)
      points(pres0,pch=3,cex=1.5,col='black')
      points(spatialToss,pch=16,cex=2.3,col='red3')
      points(envToss,pch=16,cex=1.5,col='steelblue3')
      legend('bottomleft',legend=c('good presence','spatial outlier','environmental outlier'),pch=c(3,16,16), col=c('black','red3','steelblue3'),bty='n',cex=2,pt.cex=3)
      dev.off()
    }
  }

  if(!is.null(otherSp)){
    pres1=sp::rbind.SpatialPointsDataFrame(pres0,otherSp)
  } else {pres1=pres0}

  return(list(pres=pres1,spatialToss=spatialToss,evnToss=envToss))

}



#===============================================================
# old way with all the functions wrapped together
#' @title Find outlying occurrence data
#'
#' @description Spatial or environmental outliers
#'
#' @details
#' See Examples.
#'
#' @param
# @keywords
#' @export
#'
# @examples
#'
#'
#' @return Returns a list of SpatialPointsDataFrames with (1) good presence points (2) spatial outliers and (3) environmental outliers.
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.


# findOutlyingPoints=function(pres,
# 														spOutliers=T,
# 														envOutliers=T,
# 														doPlot=FALSE,
# 														env=NULL,
# 														plotFile=NULL,
# 														species=NULL){
#
# 		#  for testing
# 		#  envOutliers=T; spOutliers=T; doPlot=T; plotFile=paste0(dirs$sp.diag.path, '/OccurrenceOutliers.png')
#
# 		keep=pres$sp==species
# 		otherSp=pres[!keep,] # append this back on later
#
# 		mcp.pres=data.frame(coordinates(pres)[keep,])
# 		coordinates(mcp.pres)=c(1,2)
# 		pe=extent(mcp.pres)
#
# 		pp=presPercentile(mcp.pres,percent=NULL)
# 		#dt=dixon.test(pp[[1]]$dist.from.centroid,opposite=F)
# 		#cot=chisq.out.test(pp[[1]]$dist.from.centroid,variance = var(pp[[1]]$dist.from.centroid),opposite = FALSE)
# 		pres.inliers=pp[[1]]$xy.t
# 		dists=pp[[1]]$dist.from.centroid
# 		sp.toss.id=NULL
# 		pval=0
# 		while(pval<1e-5){
# 			gt=outliers::grubbs.test(dists)
# 			pval=gt$p.value
# 		# conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
# 			if(gt$p.value<1e-5){
# 				toss=which.max(dists)
# 				# IDs in the original data frame
# 				sp.toss.id=c(sp.toss.id,which(pp[[1]]$dist.from.centroid == max(dists)))
# 				pres.inliers=pres.inliers[-toss,]
# 				dists=dists[-toss]
# 			}
# 		}
#
# 		# old way
# 			#== for larger sample sizes the crazy outliers are a smaller proportion. these cutoffs are just guesses.
#
# 		# 		if(length(mcp.pres)<10) { pres3=mcp.pres }
# 		# 		if(length(mcp.pres)>=25 && length(mcp.pres)<300){ pres3=presPercentile(mcp.pres,95)[[1]] }
# 		# 		if( length(mcp.pres)>=300 && length(mcp.pres)<1000 ){
# 		# 			pres3=presPercentile(mcp.pres,98)[[1]] }
# 		# 		if(length(mcp.pres)>=1000){ pres3=presPercentile(mcp.pres,99.8)[[1]] }
#
# 		env.toss.id=NULL
# 		if(envOutliers){
# 			p.env=raster::extract(env,pres)
# 			mean.p=apply(p.env,2,mean,na.rm=T)
#
# 			pres.inliers=pres
# 			dists=apply(p.env,1,function(x) sqrt(sum((x-mean.p)^2)) )
# 			env.toss.id=NULL
# 			pval=0
# 			while(pval<1e-5){
# 				gt=outliers::grubbs.test(dists)
# 				pval=gt$p.value
# 			# conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
# 				if(gt$p.value<1e-5){
# 					toss=which.max(dists)
# 					# IDs in the original data frame
# 					env.toss.id=c(env.toss.id,which(dists == max(dists)))
# 					pres.inliers=pres.inliers[-toss,]
# 					dists=dists[-toss]
# 				}
# 			}
# 		}
#
# 		if(!is.null(c(sp.toss.id,env.toss.id))){
# 			spatialToss=pres[sp.toss.id,]
# 			envToss=pres[env.toss.id,]
# 			pres0=pres[-c(sp.toss.id,env.toss.id),]
# 		} else{ pres0=pres; spatialToss=NULL; envToss=NULL}
#
# 		if(!is.null(c(sp.toss.id,env.toss.id))){ # don't plot if there's no outliers
# 			if(doPlot){
# 				if(!is.null(plotFile)) {png(plotFile,h=1000,w=1000)}
# 					ydif=pe[4]-pe[3] ; xdif=pe[2]-pe[1] # for plotting
# 					pe[1]=pe[1]-xdif; pe[2]=pe[2]+xdif; pe[3]=pe[3]-ydif; pe[4]=pe[4]+ydif
# 					env.bg=raster::crop(env[[1]],pe)
# 					plot(env.bg,col='grey50',legend=F)
# 					points(pres0,pch=3,cex=1.5,col='black')
# 					points(spatialToss,pch=16,cex=2.3,col='red3')
# 					points(envToss,pch=16,cex=1.5,col='steelblue3')
# 					legend('bottomleft',legend=c('good presence','spatial outlier','environmental outlier'),pch=c(3,16,16), col=c('black','red3','steelblue3'),bty='n',cex=2,pt.cex=3)
# 					dev.off()
# 			}
# 		}
#
# 		pres1=sp::rbind.SpatialPointsDataFrame(pres0,otherSp)
#
# 		return(list(pres=pres1,spatialToss=spatialToss,evnToss=envToss))
#
# 	}



####################################
#Segmentation functions (from http://rstudio-pubs-static.s3.amazonaws.com/10685_1f7266d60db7432486517a111c76ac8b.html)
###################################

#' Segmentation of SpatialLines object into muliple segments 
#'
#' Segment a SpatialLines object in segment of equal length (taken from: http://rstudio-pubs-static.s3.amazonaws.com/10685_1f7266d60db7432486517a111c76ac8b.html)
#' @param sl SpatialLines object
#' @param length Length of individual segment (in units of sl)
#' @param n.parts Alternatively, the number of segments to create
#' @param merge.last Whether the last segment (of a different length) should be merged with the previous segment
#' @keywords SegmentSpL CreateSegment CreateSegments MergeLast
#' @return A SpatialLines object
#' @export
#' @examples
#'  x <- c(1,5,4,8)
#'	y <- c(1,3,4,7)
#'	t1<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
#'	t2<-SegmentSpL(t1, n.parts=10, merge.last=F)
#' 	plot(t2, col = rep(c(1, 2), length.out = length(t2)), axes = T)
SegmentSpL <- function(sl, length = 0, n.parts = 0, merge.last = FALSE) {
  stopifnot((length > 0 || n.parts > 0))
  id <- 0
  newlines <- list()
  sl <- as(sl, "SpatialLines")
  for (lines in sl@lines) {
    for (line in lines@Lines) {
      crds <- line@coords
      # create segments
      segments <- CreateSegments(coords = crds, length, n.parts)
      if (merge.last && length(segments) > 1) {
        # in case there is only one segment, merging would result into error
        segments <- MergeLast(segments)
      }
      # transform segments to lineslist for SpatialLines object
      for (segment in segments) {
        newlines <- c(newlines, Lines(list(Line(unlist(segment))), ID = as.character(id)))
        id <- id + 1
      }
    }
  }
  return(SpatialLines(newlines))
}

CreateSegment <- function(coords, from, to) {
  distance <- 0
  coordsOut <- c()
  biggerThanFrom <- F
  for (i in 1:(nrow(coords) - 1)) {
    d <- sqrt((coords[i, 1] - coords[i + 1, 1])^2 + (coords[i, 2] - coords[i + 
                                                                             1, 2])^2)
    distance <- distance + d
    if (!biggerThanFrom && (distance > from)) {
      w <- 1 - (distance - from)/d
      x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
      y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
      coordsOut <- rbind(coordsOut, c(x, y))
      biggerThanFrom <- T
    }
    if (biggerThanFrom) {
      if (distance > to) {
        w <- 1 - (distance - to)/d
        x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
        y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
        coordsOut <- rbind(coordsOut, c(x, y))
        break
      }
      coordsOut <- rbind(coordsOut, c(coords[i + 1, 1], coords[i + 1, 
                                                               2]))
    }
  }
  return(coordsOut)
}

CreateSegments <- function(coords, length = 0, n.parts = 0) {
  stopifnot((length > 0 || n.parts > 0))
  # calculate total length line
  total_length <- 0
  for (i in 1:(nrow(coords) - 1)) {
    d <- sqrt((coords[i, 1] - coords[i + 1, 1])^2 + (coords[i, 2] - coords[i + 
                                                                             1, 2])^2)
    total_length <- total_length + d
  }
  
  # calculate stationing of segments
  if (length > 0) {
    stationing <- c(seq(from = 0, to = total_length, by = length), total_length)
  } else {
    stationing <- c(seq(from = 0, to = total_length, length.out = n.parts), 
                    total_length)
  }
  
  # calculate segments and store the in list
  newlines <- list()
  for (i in 1:(length(stationing) - 1)) {
    newlines[[i]] <- CreateSegment(coords, stationing[i], stationing[i + 
                                                                       1])
  }
  return(newlines)
}
MergeLast <- function(lst) {
  l <- length(lst)
  lst[[l - 1]] <- rbind(lst[[l - 1]], lst[[l]])
  lst <- lst[1:(l - 1)]
  return(lst)
}

####################################
# Intersect 
####################################

#' Intersection between individual trajectory and segmented line  
#'
#' Calculate summary statistics regarding frequency of intersections between an animal trajectory and a segmented SpatialLines object
#' @param traj An animal trajectory of class ltraj
#' @param corri A segmented SpatialLines as returned bt function SegmentSpL
#' @param per_day Specify if the number of crossing should be standardized over the time (in days) the animal was monitored or over the total number of locations, default=TRUE. Long gap in data should be identified as separate burst
#' @param plot Whether a plot should be returned (default=F)
#' @keywords SegmentSpL ltraj plotcorri_ind
#' @return A SpatialLinesDataFrame object
#' @export
#' @examples
#'	require(adehabitatLT)
#'  x <- c(0,0)
#'	y <- c(-6500000,-4500000)
#'	t1<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
#'	t2<-SegmentSpL(t1, n.parts=20, merge.last=F)
#'	data (albatross) #From package adehabitatLT
#'	t3<-corriIntersects(albatross[3], t2, plot=T)
#'  plot(ltraj2sldf(albatross[3]), add=T)
corriIntersects<-function(traj, corri, per_day=TRUE, plot=F) {
  data<-adehabitatLT::ld(traj)
  data2<-na.omit(cbind(data[,1],data[,2], data[,1]+data[,4],data[,2]+data[,5]))
  steps<-SpatialLines(apply(data2, 1, function(r) {
    Lines(list(Line(cbind(r[c(1,3)], r[c(2,4)]))), uuid::UUIDgenerate())
  }))
  inter<-rgeos::gIntersects(corri, steps, byid=T)
  pct<-colSums(inter)/(sum(na.omit(data)$dt)/3600/24)
  if(per_day==F) {pct<-colSums(inter)/nrow(inter)*100}
  out<-SpatialLinesDataFrame(corri, data.frame(count=colSums(inter), pct=pct, cross=ifelse(colSums(inter)>0, 1, 0)))
  if(plot==T) {plotcorri_ind(out, extent=raster::extent(steps))}
  return(out)  
}

#' Plot of density of crossing of an individual with linear features 
#'
#' Produce a color-coded plot of density of crossing of an individual with a segmented linear features. Used to display result of function corriIntersects
#' @param SpL A SpatialLinesDataFrame object returned by corriIntersects
#' @param nb_breaks The number of breaks to use in display (default =5), must be <=10
#' @param extent If not NULL, an extent object (from package raster) that specify boundaries of the display. 
#' @keywords corriIntersects
#' @return A plot
#' @export
#' @examples
#'	require(adehabitatLT)
#'  x <- c(0,0)
#'	y <- c(-6500000,-4500000)
#'	t1<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
#'	t2<-SegmentSpL(t1, n.parts=20, merge.last=F)
#'	data (albatross) #From package adehabitatLT
#'	t3<-corriIntersects(albatross[3], t2, plot=F)
#'	plotcorri_ind(t3, nb_breaks=4, extent=raster::extent(ltraj2spdf(albatross)))
plotcorri_ind<-function(SpL, nb_breaks=5, extent=NULL) {
  breaks<-seq(min(SpL$pct)-1e-16, max(SpL$pct), length.out=nb_breaks)
  colors <- c("yellowgreen", "yellow", "orange", "orangered", "orangered4", "red", "red4", "maroon4","slateblue4", "navy")
  liwd<-c(1,3,3,3,3,3,3,3,3,3)
  plot(SpL, col=colors[cut(SpL$pct, breaks)], lwd=liwd[cut(SpL$pct, breaks)])
  if(!is.null(extent)) { plot(SpL, col=colors[cut(SpL$pct, breaks)], lwd=liwd[cut(SpL$pct, breaks)], xlim=c(extent[1], extent[2]), ylim=c(extent[3], extent[4]) )}
  legend("topleft",legend=levels(cut(SpL$pct, breaks)), col=colors[1:nb_breaks], lwd=liwd[1:nb_breaks])
}
  
###################################
## Combine individuals 
###################################

#' Intersection between multiple trajectories and segmented line  
#'
#' Wrapper function that apply corriIntersects to all individual in a trajectory object
#' @param traj An animal trajectory of class ltraj containing multiple individuals
#' @param corri A segmented SpatialLines as returned bt function SegmentSpL
#' @keywords SegmentSpL ltraj plotcorri_ind
#' @return A list of SpatialLinesDataFrame object
#' @export
#' @examples
#'	require(adehabitatLT)
#'  x <- c(0,0)
#'	y <- c(-6500000,-4500000)
#'	t1<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
#'	t2<-SegmentSpL(t1, n.parts=20, merge.last=F)
#'	data (albatross) #From package adehabitatLT
#'	t3<-corriIntersects_All(albatross, t2) 
corriIntersects_All<-function(trajs, corri) {
  id<-unique(id(trajs))
  id2<-id(trajs)
  return(lapply(1:length(id), function(x) corriIntersects(trajs[which(id2==id[x])], corri)))
}

#' Averaged intersection between multiple trajectories and segmented line  
#'
#' Combined a list list of corriIntersects object to obtain population level averaged statistics regarding crossing count
#' @param SpLlst A list of animal trajectory of class ltraj containing multiple individuals
#' @keywords SegmentSpL ltraj plotcorri_ind
#' @return A  SpatialLinesDataFrame object
#' @export
#' @examples
#'	require(adehabitatLT)
#'  x <- c(0,0)
#'	y <- c(-6500000,-4500000)
#'	t1<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
#'	t2<-SegmentSpL(t1, n.parts=20, merge.last=F)
#'	data (albatross) #From package adehabitatLT
#'  t3<-corriIntersects_All(albatross, t2) 
#'	t4<-avg_inds(t3) 
avg_inds<-function(SpLlst) {
  count_mean<-rowMeans(sapply(SpLlst, function(x) x$count))
  count_sum<-rowSums(sapply(SpLlst, function(x) x$count))
  pct_mean_0<-rowMeans(sapply(SpLlst, function(x) x$pct))
  tt<-sapply(SpLlst, function(x) x$pct)
  tt<-ifelse(tt==0, NA, tt)
  pct_mean<-rowMeans(tt, na.rm=T)
  pct_mean[is.na(pct_mean)]<-0
  nb_ind<-rowSums(sapply(SpLlst, function(x) x$cross))
  tt2<-apply(sapply(SpLlst, function(x) x$cross), 1, function(x) which(x>0))
  tt3<-lapply(tt2, length)
  id<-ifelse(tt3==0, "None", "Multi")
  id[which(tt3==1)]<-unlist(tt2[which(tt3==1)])
  return(SpatialLinesDataFrame(SpLlst[[1]], data.frame(count_mean, count_sum, pct_mean_0, pct_mean, nb_ind, id), match.ID=F))
}



#' Plot of density of crossing of multiple individuals with linear features 
#'
#' Produce a color-coded plot of density of crossing of multiple individuals with a segmented linear features. Used to display result of function avg_inds
#' @param SpL A SpatialLinesDataFrame object returned by avg_inds
#' @param nb_breaks The number of breaks to use in display (default =5), must be <=10
#' @param var Variable to display (default =4) 1= Average count, 2=Total count, 3=Average count standardized by length of tracking, 4=Average count standardized by length of tracking and spatial sampling bias, 5= Number of different individuals  
#' @keywords corriIntersects_All avg_inds
#' @return A plot
#' @export
#' @examples
#'	require(adehabitatLT)
#'  x <- c(0,0)
#'	y <- c(-6500000,-4500000)
#'	t1<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
#'	t2<-SegmentSpL(t1, n.parts=20, merge.last=F)
#'	data (albatross) #From package adehabitatLT
#'	t3<-corriIntersects_All(albatross, t2)
#'	t4<-avg_inds(t3)  
#'  plotcorri_grp(t4, nb_breaks=5, var=4, main="Albatross")
plotcorri_grp<-function(SpL, nb_breaks=5, var=4, main="Default") {
  breaks=seq(min(SpL@data[,var])-1e-16, max(SpL@data[,var]), length.out=nb_breaks)
  colors <- c("yellowgreen", "yellow", "orange", "orangered", "orangered4", "red", "red4", "maroon4","slateblue4", "navy")
  liwd<-c(1,3,3,3,3,3,3,3,3,3)
  plot(SpL, col=colors[cut(SpL@data[,var], breaks)], lwd=liwd[cut(SpL@data[,var], breaks)], main=main)
  legend("topleft",legend=levels(cut(SpL@data[,var], breaks)), col=colors[1:nb_breaks], lwd=liwd[1:nb_breaks])
}



#####################################
## Homerange functions based on http://r-sig-geo.2731867.n2.nabble.com/Split-polygon-by-line-td7588625.html
#####################################

#' Intersection of a home-range (polygon) with linear features   
#'
#' Intersects a SpatialPolygons* object with a SpatialLines* object and return a divided polygon (when overlap) with summary metric regarding % of area left (bigger sectio), number of segments created, and total initial area of the polygon. Taken from: http://r-sig-geo.2731867.n2.nabble.com/Split-polygon-by-line-td7588625.html
#' @param pol A SpatialPolygons* object representing an animal home-range. 
#' @param line A linear feature of class SpatialLines
#' @keywords SpatialPolygons SpatialLines
#' @return A list containing the segmented polygon, the % area of the bigger segment, the number of segment, and polygon initial total area 
#' @export
#' @examples
#' pol<-SpatialPolygons(list(Polygons(list(Polygon(cbind(c(0,1,1,0,0),c(0,0,1,1,0)))), ID="polygon"))) 
#' line<-SpatialLines(list(Lines(list(Line(cbind(c(0,1),c(0.4,0.4)))), ID="line"))) 
#' splt<-hr_split(pol, line)
#' plot(splt[[1]])
#' splt[2:4] 
hr_split<-function(pol=hr, line=corri) {
    lpi <- rgeos::gIntersection(pol, line) # intersect line withthe polygon 
  blpi <- rgeos::gBuffer(lpi, width = 0.000001)  # create a very thin polygon 
  #buffer of the intersected line 
  dpi <- rgeos::gDifference(pol, blpi) # split using gDifference 
  area<-unlist(lapply(dpi@polygons[[1]]@Polygons, function(x) x@area))/1000/1000
  pct_left<-max(area)/sum(area)*100
  nb_seg<-length(area)
  out<-list(dpi, pct_left, nb_seg, area)
  names(out)<-c("Spatial object", "% of bigger fragment", "Nb of fragment", "Area of each fragment")
  return(out)
}

#######################################
### Misc functions 
#######################################

#' Range standardisation (0,1)  
#'
#' This function standardises a vector between 0 and 1
#' @param x A vector
#' @keywords range standardisation
#' @export
#' @examples
#' v<-c(1,2,2,3,4,4,5,6)
#' range01(v)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


#' Associate closest elements between two SpatialPoints features
#'
#' Evaluate the distance between a SpatialPoints/Lines* object and another SpatialPoints* object and return which features in first object that are the closest to each feature in second object
#' @param pts1 A SpatialLines* or SpatialPoints object. If a SpatialLines is provided, the object is first converted to point using getSpatialLinesMidPoints
#' @param pts2 A SpatialPoints* object 
#' @keywords match SpatialPoints
#' @export
#' @examples
#'  x <- c(1,5,4,8)
#'	y <- c(1,3,4,7)
#'	t1<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
#'	t2<-SegmentSpL(t1, n.parts=10, merge.last=F)
#'  Pts<-SpatialPoints(matrix(c(1,3,7,1,2,7), nrow=3, ncol=2))
#'  match_pts(t2, Pts)
match_pts<-function(pts1, pts2) {
  if(class(pts1)=="SpatialLinesDataFrame" | class(pts1)=="SpatialLines") {pts1<-getSpatialLinesMidPoints(pts1)}
  t1<-rgeos::gDistance(pts1, pts2, byid=T)
  t2<-unique(unlist(apply(t1, 1, function(x) which(x==min(x))[1])))
  return(t2)
}


######################################
### Optimization function
######################################

#' Optimization of wildlife crossing locations over a linear features
#' 
#' The function use linear programming to optimize the location of wildlife crossing over a linear feature. The function maximise the spatial spread of locations, 
#'    and the importance of specific location for animal crossing. The user can specify the number of crossing location desired, if some segment should be excluded, or if the location
#'     of some crossings are already decided. The user also need to specify the weight given to the spatial argument and the importance of crossing (default is equal importance to each).
#' 
#' @param corri A segmented SpatialLines object returned by avg_inds  
#' @param var Variable used to represent importance of segment for animal. Number refer to columns of dataframe produced by avg_inds (default = 4, average percent crossing) 
#' @param n Number of crossing to place (default =5) in addition to fixed points (i.e. if a SpatialPoints* object is provided to the add argument)
#' @param pct_keep Percentage of segment to consider in optimization. Removal is based on percentile of values of importance of segment for animal crossing
#' @param rm A SpatialPoints object indicating the location that most be excluded from the optimization (defaul is NULL)
#' @param add A SpatialPoints object indicating the location where a crossing is already present, or must be place to this location. Number of points included here will be added to n to give the total number of crossing selected. 
#' @param cost A vector of same length that the number of segment in corri giving the cost for each segment
#' @param weight Argument setting the weight given to the spatial spread relative to the crossing importance (default = 0.5 meaning equal importance)
#' @param ln Whether the natural logarithm of the variable value should be taken. Default=F
#' @param plot Whether a plot showing the crossing location should be returned (default=T)
#' @param ... additional arguments that can be specify to Rsymphony_solve_lp (time_limit, gap_limit, first_feasible)
#' @keywords SpatialPolygons SpatialLines
#' @return A list containing the segmented polygon, the % area of the bigger segment, the number of segment, and polygon initial total area 
#' @export
#' @examples
#'	require(adehabitatLT)
#'  x <- c(0,0)
#'	y <- c(-6500000,-4500000)
#'	t1<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
#'	t2<-SegmentSpL(t1, n.parts=20, merge.last=F)
#'	data (albatross) #From package adehabitatLT
#'  t3<-corriIntersects_All(albatross, t2) 
#'	t4<-avg_inds(t3)
#'  opti1<-optim_corri(t4, var=4, n=3, nb_ind=1, weight=0.5, plot=T) #Equal weight
#'  opti2<-optim_corri(t4, var=4, n=3, nb_ind=1, weight=0.95, plot=T) #More importance to spatial spread

optim_corri<-function(corri, var=4, n=5, pct_keep=1, nb_ind=1, ln=F, rm=NULL,add=NULL, weight=0.5, plot=T, time_limit=-1, gap_limit = -1, first_feasible=F, ...) {
  pts<-getSpatialLinesMidPoints(corri) 
  dmat<-weight*range01(as.matrix(dist(coordinates(pts))))
  diag(dmat)<-(1-weight)*range01(corri@data[,var])*2*n #Multiply by 2 bc should be twice the rest of spatial bc spatial has spread and distance 
  #tt<-which(diag(dmat) == 0)
  tt<-which(corri@data[,5]<nb_ind)
  #Condition if nb_ind = 1 then keep 1 segment per ind 
  if(nb_ind==1) { 
    sub<-which(corri@data[,5]==1)
    id<-corri@data[sub,6]
    val<-corri@data[sub,var]
    df<-data.frame(sub, id, val)
    df<-df[order(df$val),]
    df<-df[!duplicated(df$id, fromLast=T), ]
    sub2<-sub[!sub %in% df$sub]
    tt<-unique(c(tt, sub2))
    }
   if (pct_keep !=1)   {tt2<-which(diag(dmat)<= quantile(diag(dmat)[diag(dmat)>0], probs=pct_keep))
                    tt<- unique(c(tt, tt2))} 
  #Remove unavailable pts 
  if(!is.null(rm)) {
    pts_rm<-match_pts(corri, rm)
    tt<-unique(c(tt, pts_rm))
  }
  #Find fixed segments 
  lgth=0
  fixed<-rep(0, nrow(dmat))[-tt]
  if(!is.null(add)) {
    pts_fixed<-match_pts(corri, add)
    #if (sum(tt %in% pts_fixed)>=1) {tt<-tt[-which(tt %in% pts_fixed)]} # Make sure the one that are fixed are kept in calculation 
   fixed<-rep(0, nrow(dmat))#Vector of 1 where there are fixed otherwise zero 
    fixed[pts_fixed]<-1
    fixed<-fixed[-tt]
    lgth<-sum(fixed)
  }
  n<-n+lgth
  #Remove from dmat
  kept<-c(1:length(pts))[-tt]
  dmat<-dmat[-tt,-tt]
  if(ln==T) {diag(dmat)<-(1-weight)*range01(log(corri@data[,var][-tt]))*2*n}
  
  #Matrix of constraints 
  dmat2<-matrix(0, nrow=nrow(dmat), ncol=ncol(dmat))
  diag(dmat2)<-1
  lg<-sum(upper.tri(dmat, diag=T))
  r1<-c(dmat2[upper.tri(dmat2, diag=T)], rep(0, nrow(dmat)*1)) # diagonal only constraints 
  r2<-c(rep(1, lg), rep(0,nrow(dmat)*1)) #Maximum number of pixels occupied = n*n
  r3<-c(rep(0, lg),  rep(1, nrow(dmat))) # avg dist constraints = n 
 r4<-c(rep(0, lg),  fixed) # For the one that are fixed
  dmat3<-matrix(0, nrow=nrow(dmat), ncol=ncol(dmat))
  dmat4<-matrix(0, nrow=nrow(dmat),ncol=sum(upper.tri(dmat, diag=T)))
  for (i in 1:nrow(dmat)) {
    dmat5<-dmat3
    dmat5[i,]<-1
    dmat5[,i]<-1
    dmat4[i,]<-dmat5[upper.tri(dmat5, diag=T)]
   }
  dmat5<-matrix(0, nrow=nrow(dmat), ncol=ncol(dmat)) # Y side of constraints for 
  diag(dmat5)<--n
  mat<-rbind(r1,r2,r3,r4, cbind(dmat4, dmat5)) # cbind(dmat4, dmat5), cbind(dmat3, dmat6, dmat5), cbind(dmat4, dmat6, dmat5))
  #Direction of equalities/inequalities
  dir<-c("==", "==", "<=","==", rep("==", nrow(dmat4)*1))
  #Right hand side 
  rhs<-c(n,sum(upper.tri(matrix(0, nrow=n, ncol=n), diag=T)), n, lgth, rep(0, nrow(dmat4)*1))
  #Objective function
  obj<-c(as.numeric(dmat[upper.tri(dmat, diag=T)]),  n*(weight-(colMeans(dmat))))  #Step distance divided by 2 bc added twice, avg distance multiplied by number of locs,  
  #Linear prog solver 
  opti<-Rsymphony::Rsymphony_solve_LP(obj, mat, dir, rhs, types="B", max=T, time_limit = time_limit, gap_limit=gap_limit, first_feasible=first_feasible)
  tt2<-tail(opti$solution, ncol(dmat))
  t1<-kept[which(tt2>0)]
  if(!is.null(add)) {t2<-kept[which(fixed>0)]}
  if(plot==T) {
    plotcorri_grp(corri, var=var)
    plot(pts[t1], col="blue", add=T)
    if(!is.null(add)) {plot(pts[t2], col="red", add=T)}
  }
  if(!is.null(add)) {return(list(opti, pts[t1], pts[t2]))}
  if(is.null(add)) {return(list(opti, pts[t1]))}
}

#' Plot of selected crossing structures following optimization using ILP
#'
#' Produce a color-coded plot of density of crossing of multiple individualsand optimal crossing structure. Used to display result of function optim_corri
#' @param corri A SpatialLinesDataFrame object returned by avg_inds
#' @param var Variable to display (default =4) 1= Average count, 2=Total count, 3=Average count standardized by length of tracking, 4=Average count standardized by length of tracking and spatial sampling bias, 5= Number of different individuals  
#' @param optim Output of the optimization function optim_corri
#' @keywords corriIntersects_All avg_inds optim_corri
#' @return A plot
#' @export
#' @examples
#'	require(adehabitatLT)
#'  x <- c(0,0)
#'	y <- c(-6500000,-4500000)
#'	t1<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
#'	t2<-SegmentSpL(t1, n.parts=20, merge.last=F)
#'	data (albatross) #From package adehabitatLT
#'	t3<-corriIntersects_All(albatross, t2)
#'	t4<-avg_inds(t3) 
#'  opti1<-optim_corri(t4, var=4, n=3, nb_ind=1, weight=0.5, plot=T) #Equal weight 
#'  plot_optim(t4,  var=4, opti1, main="Crossings") 
plot_optim<-function(corri, var=4, optim, main="Default", ...){
  plotcorri_grp(corri, var, main=main)
  plot(optim[[2]], col="blue", add=T,pch=4, ps=2)
  if(length(optim)==3) {plot(optim[[3]], col="red",pch=4, add=T)}
}





###### NEW OPTIMIZATION FUNCTION 
#' Optimization of wildlife crossing locations over a linear features using the maximum coverage location problem
#' 
#' The function use linear programming to optimize the location of wildlife crossing over a linear feature. The function maximise the number of individuals the selected features will assist.
#'  The user can specify the number of crossing location desired, a coverage distance,  if some segment should be excluded, or if the location
#'     of some crossings are already decided. 
#' 
#' @param corri_ls A segmented SpatialLines list returned by corriIntersects_All  
#' @param n Number of crossing to place (default =5) in addition to fixed points (i.e. if a SpatialPoints* object is provided to the add argument)
#' @param dist Distance used to considered a given segment as covered (ie radius)
#' @param rm A SpatialPoints object indicating the location that most be excluded from the optimization (defaul is NULL)
#' @param add A SpatialPoints object indicating the location where a crossing is already present, or must be place to this location. Number of points included here will be added to n to give the total number of crossing selected. 
#' @param plot Whether a plot showing the crossing location should be returned (default=T)
#' @param ... additional arguments that can be specify to Rsymphony_solve_lp (time_limit, gap_limit, first_feasible)
#' @keywords SpatialPolygons SpatialLines
#' @return A list containing the segmented polygon, the % area of the bigger segment, the number of segment, and polygon initial total area 
#' @export
#' @examples
#'  require(adehabitatLT)
#'  x <- c(0,0)
#'  y <- c(-6500000,-4500000)
#'  t1<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
#'  t2<-SegmentSpL(t1, n.parts=20, merge.last=F)
#'  data (albatross) #From package adehabitatLT
#'  t3<-corriIntersects_All(albatross, t2) 
#'  test<-optim_mclp(t3, n=2,dist=10*1000, plot=T)


optim_mclp<-function(corri_ls, n=5, dist=10*1000,rm=NULL,add=NULL, plot=T, time_limit=-1, gap_limit = -1, first_feasible=F, ... ) {
t4<-avg_inds(corri_ls)  
t5<-getSpatialLinesMidPoints(t4)
#List of not empty segments
seg<-which(t4$nb_ind>0)
#Segment to exclude 
if(!is.null(rm)) {
  pts_rm<-match_pts(t4, rm)
tt<-which(seg %in% pts_rm)
if (length(tt)>0) {seg<-seg[-tt]}
  }
#Matrix of constraints
fct<-function(x) { 
  tt<-vector()
  for (i in 1:length(seg)) {
    tt[i]<-ifelse(sum(spDistsN1(x, t5[seg[i]])<dist)>0,1,0)
  }
  return(tt)
}

mat1<-matrix(unlist(lapply(corri_ls, function(x) fct(getSpatialLinesMidPoints(x)[which(x$cross==1)]))), nrow=length(corri_ls), ncol=length(seg), byrow=T)

#Segment fixed 
rw<-rep(0,ncol(mat1))
l<-0
if(!is.null(add)) {
  pts_add<-match_pts(t4, add)
  tt2<-which(seg %in% pts_add)
  n<-n+length(tt2)
  rw[tt2]<-1
  l<-length(tt2)
}

mat2<-matrix(0, nrow=nrow(mat1), ncol=nrow(mat1))
diag(mat2)<--1
mat3<-cbind(mat1,mat2)
#mat2<-rbind(rep(1, ncol(mat1)),rw, mat1)
#mat3<-cbind(mat2, c(0,0, rep(-1, nrow(mat1))))
mat4<-rbind(c(rep(1, ncol(mat1)),rep(0, nrow(mat2))),c(rw, rep(0, nrow(mat2))), mat3)
#Other parameters
dir<-c("==", "==", rep(">=", nrow(mat1)))
rhs<-c(n, l, rep(0, nrow(mat1)))
#obj<-c(colSums(mat1),0)
obj<-c(rep(0, ncol(mat1)), rep(1, nrow(mat1)))
  
  
#Optimization
opti<-Rsymphony::Rsymphony_solve_LP(obj, mat4, dir, rhs, types="B", max=T, time_limit = time_limit, gap_limit=gap_limit, first_feasible=first_feasible)

#Output results 
out<-opti$solution
t1<-seg[which(out[1:ncol(mat1)]>0)]
if(!is.null(add)) {t2<-seg[which(rw>0)]}
if(plot==T) {
  plotcorri_grp(t4, var=5)
  plot(t5[t1], col="blue", add=T)
  if(!is.null(add)) {plot(t5[t2], col="red", add=T)}
}
if(!is.null(add)) {return(list(opti, t5[t1], t5[t2]))}
if(is.null(add)) {return(list(opti, t5[t1]))}
}



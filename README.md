## wildxing  ##

**wildxing** is an R package for optimal positioning of wildlife crossing structures using GPS telemetry.

This is the development area for the package `wildxing`, which provides a series of functions to analyze movement data to optimally locate wildlife crossing structures. 

*References*: To be submitted. 

For questions: gbr |at| colostate.edu

## Installation of the development version  ##

You need to use the package `devtools`from Hadley Wickham. 
    
    library(devtools)
    install_github("BastilleRousseau/wildxing")


## Getting started ##

The package main functions are `corriIntersects`, `hr_split`, `avg_inds` and `optin_corri`.  For a list of documented functions see the Reference manual. 
For examples of how to use the main functions, you can also look at the vignette. 

Alternatively, here is a quick example to get you going: 
  
    require(adehabitatLT)
    x <- c(0,0)
    y <- c(-6500000,-4500000)
    t1<-SpatialLines(list(Lines(Line(cbind(x,y)), ID="a")))
    t2<-SegmentSpL(t1, n.parts=20, merge.last=F)
    data (albatross) #From package adehabitatLT
    t3<-corriIntersects_All(albatross, t2) 
    t4<-avg_inds(t3)
    opti1<-optim_corri(t4, var=4, n=3, nb_ind=1, weight=0.5, plot=T) #Equal weight
    opti2<-optim_corri(t4, var=4, n=3, nb_ind=1, weight=0.95, plot=T) #More importance to spatial spread

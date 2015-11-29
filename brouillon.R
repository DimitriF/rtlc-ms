setwd('/home/dimitri/rMS/www')

# library("BiocInstaller")
# biocLite("RforProteomics", dependencies = TRUE)
library(sqldf)
library(enviPick)
library(tidyr)
library(abind)

# library("mzR")


system.time(df <- mass2R(c('propolis pl1 tr1 RP.mzXML','propolis pl1 tr1 RP.mzXML'),900))
data <- mass2R(c('propolis pl1 tr1 RP.mzXML'),900)
for(i in seq(length(data))){
  data[[i]] <- data[[i]][,as.numeric(c(163,255,285,177))]
}
Smoot <- list(window.size = 9,poly.order=1,diff.order=0)
data.1 <- lapply(data,function(x){
  t(f.preprocess(t(x),preprocess.order=c('Smoothing'),
                 preprocess.option=list(Smoothing=Smoot),
                 training.data=t(x)))}
)
as.data.frame(lapply(data.1,function(x){apply(x,2,sum)}))

## Usable to peak detection, index are available , look for area
baseline.peakDetection(t(data[[1]]),5,50,5,5,100)$peaks

# library("mzID")
# library("MSnID")
# library("MSnbase")
# library("rpx")
# library("MLInterfaces")
# library("pRoloc")
# library("pRolocdata")
# library("MSGFplus")
# library("rols")
# library("hpar")

# work but seems limited
library(mzR)
aa <- openMSfile('propolis pl1 tr1 RP.mzXML','Ramp')

pl <- peaks(aa,1)
head(pl,10)
plot(pl[,1],pl[,2],type='h')

col <- runInfo(aa)[[1]]
eic <- c()
range <- c(0,1000)
for(i in seq(col)){
  truc <- peaks(aa,i)
  truc <- sum(truc[truc[,1] < range[2] & truc[,1] > range[1],2])
  eic <- c(eic,truc)
}
plot(eic,type='l')


# produce the list
aa <- openMSfile('propolis pl1 tr1 RP.mzXML','Ramp')
col <- runInfo(aa)[[1]]
ls.1 <- list()
for(i in seq(col)){
  ls.1[[i]] <- peaks(aa,i)
}
# get the min and max
mz.min <- round(runInfo(aa)[[2]])
mz.max <- round(runInfo(aa)[[3]])
# round the mz
ls.2 <- lapply(ls.1,function(x){x[,1]<-round(x[,1]);return(x)})
# change to vector
ls.3 <- lapply(ls.2,
               function(x){
                 vec <- rep(0,mz.max)
                 for(i in seq(nrow(x))){
                   vec[x[i,1]]<-x[i,2]
                 }
                 return(vec)
               }
)
# create the df
df <- as.matrix(as.data.frame(ls.3))

dim(df)
df <- df[1:300,seq(1,945,by=3)]
df.2 <- round(df/10000)

library(rgl)
require("RColorBrewer")
# nbcol = 100
# color = rev(heat.colors(nbcol,alpha=1))
# zcol  = cut(df, nbcol)
persp3d(y = seq(from=1,to=945,len=ncol(df)),
        x = seq(mz.min, 400, len = nrow(df)),
        z=df,xlab="mz",ylab="spectrum number",zlab="",
        col='grey',back="lines")


# RforProteomics
ms <- openMSfile('propolis pl1 tr1 RP.mzXML','Ramp')
hd <- header(ms)
## a set of spectra of interest: MS1 spectra eluted
## between 30 and 35 minutes retention time
ms1 <- which(hd$msLevel == 1)
rtsel <- hd$retentionTime[ms1] / 60 > 0 &
  hd$retentionTime[ms1] / 60 < 10

## the map
M <- MSmap(ms, ms1[rtsel], 100, 700, .1, hd)

## 1

plot(M, aspect = 1, allTicks = FALSE)
plot3D(M)


# don't work
library(affyio)
data <- read.cdffile.list('MaryCT82ng.cdf')

# test
library(xcms)
xr <- xcmsRaw("propolis pl1 tr1 RP.mzXML")
plot(xr@tic,type='l')
getPeaks(xr,peakrange=data.frame(mzmin=300,mzmax=400,rtmin=min(xr@scantime),rtmax=max(xr@scantime)),step=0.1)


# test: working but strangely time consuming
library("readMzXmlData")
exampleDirectory <- system.file("Examples", package="readMzXmlData")
spec <- readMzXmlDir(exampleDirectory)
plot(spec[[1]]$spectrum$mass, spec[[1]]$spectrum$intensity, type="n")
l <- length(spec)
legendStr <- character(l)
for (i in seq(along=spec)) {
  lines(spec[[i]]$spectrum$mass, spec[[i]]$spectrum$intensity, type="l",
        col=rainbow(l)[i])
  legendStr[i] <- basename(spec[[i]]$metaData$file)
}
spec <- readMzXmlFile('propolis pl1 tr1 RP.mzXML')
chrom <- lapply(spec,function(x){as.data.frame(x[[1]])})

# enviPick
library(enviPick)
df <- readMSdata('propolis pl1 tr1 RP.mzXML')$Scans[[2]]
df[,1] <- round(df[,1])
df <- as.data.frame(df[,1:3])
colnames(df) <- c('mz','intensity','rt')
## sql: super
df.2 <- sqldf('select mz, max(intensity), rt from df group by mz, rt')
colnames(df.2) <- c('mz','intensity','rt')
# library(plyr)
# truc <- daply(df, .(mz, rt), function(x){max(x$intensite)})
library(tidyr)
truc <- spread(df.2, mz, intensity,fill=0)

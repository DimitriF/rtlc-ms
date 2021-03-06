f.plot.array<-function(data,id,label,hauteur,Zf,dist.bas,reconstruct=T,xlim=c(-dist.bas/(Zf-dist.bas),(hauteur-dist.bas)/(Zf-dist.bas)),inverse=F,ylim.raster=1.3){
  if(reconstruct==T){
    if(is.null(label)){
      plot(x=seq((hauteur-dist.bas)/(Zf-dist.bas),-dist.bas/(Zf-dist.bas),length.out=dim(data)[2]),y=as.vector(data[id,,1]),
           ylim=c(0,ylim.raster),xlim=xlim,xlab='',ylab="",
           type="n")
      par(new=T)
    }else{
      plot(x=seq((hauteur-dist.bas)/(Zf-dist.bas),-dist.bas/(Zf-dist.bas),length.out=dim(data)[2]),y=as.vector(data[id,,1]),
           ylim=c(0,ylim.raster),xlim=xlim,main=label[id],xlab='',ylab="",
           type="l",col="red")
      par(new=T)
    }
    plot(x=seq((hauteur-dist.bas)/(Zf-dist.bas),-dist.bas/(Zf-dist.bas),length.out=dim(data)[2]),y=as.vector(data[id,,1]),
         ylim=c(0,ylim.raster),xlim=xlim,xlab='',ylab="",
         type="l",col="red")
    par(new=T)
    plot(x=seq((hauteur-dist.bas)/(Zf-dist.bas),-dist.bas/(Zf-dist.bas),length.out=dim(data)[2]),y=as.vector(data[id,,2]),
         ylim=c(0,ylim.raster),xlim=xlim,xlab='',ylab='',
         type="l",col="green")
    par(new=T)
    plot(x=seq((hauteur-dist.bas)/(Zf-dist.bas),-dist.bas/(Zf-dist.bas),length.out=dim(data)[2]),y=as.vector(data[id,,3]),
         ylim=c(0,ylim.raster),xlim=xlim,xlab='',ylab='',
         type="l",col="blue")
    par(new=T)
    plot(x=seq((hauteur-dist.bas)/(Zf-dist.bas),-dist.bas/(Zf-dist.bas),length.out=dim(data)[2]),y=as.vector(data[id,,4]),
         ylim=c(0,ylim.raster),xlim=xlim,xlab='',ylab='',
         type="l",col="black")
    if(inverse==F){
      data.plot<-round(array(data[id,,c(1,2,3)],dim=c(1,dim(data)[2],3))*256)/256
      rasterImage(data.plot,(hauteur-dist.bas)/(Zf-dist.bas) , 1.1, -dist.bas/(Zf-dist.bas), 1.3)
    }else{
      data.plot<-round(array(data[id,dim(data)[2]:1,c(1,2,3)],dim=c(1,dim(data)[2],3))*256)/256
      rasterImage(data.plot, (hauteur-dist.bas)/(Zf-dist.bas) , 1.1, -dist.bas/(Zf-dist.bas), 1.3)
    }
  }else{
    plot(x=seq((hauteur-dist.bas)/(Zf-dist.bas),-dist.bas/(Zf-dist.bas),length.out=dim(data)[2]),y=as.vector(data[id,,1]),
         ylim=c(min(data),max(data)),xlim=xlim,main=label[id],xlab='',ylab='',
         type="l",col="red")
    par(new=T)
    plot(x=seq((hauteur-dist.bas)/(Zf-dist.bas),-dist.bas/(Zf-dist.bas),length.out=dim(data)[2]),y=as.vector(data[id,,2]),
         ylim=c(min(data),max(data)),xlim=xlim,xlab='',ylab='',
         type="l",col="green")
    par(new=T)
    plot(x=seq((hauteur-dist.bas)/(Zf-dist.bas),-dist.bas/(Zf-dist.bas),length.out=dim(data)[2]),y=as.vector(data[id,,3]),
         ylim=c(min(data),max(data)),xlim=xlim,xlab='',ylab='',
         type="l",col="blue")
    par(new=T)
    plot(x=seq((hauteur-dist.bas)/(Zf-dist.bas),-dist.bas/(Zf-dist.bas),length.out=dim(data)[2]),y=as.vector(data[id,,4]),
         ylim=c(min(data),max(data)),xlim=xlim,xlab='',ylab='',
         type="l",col="black")
  }
  mtext(side = 1, expression("R"['F']), line = 2)
  mtext(side = 2, "Intensity", line = 2)
}
mass2R <- function(files,cut.off=900){
  ms.list <- list()
  for(j in seq(length(files))){
    df <- readMSdata(files[j])$Scans[[2]]
    df[,1] <- round(df[,1])
    df <- as.data.frame(df[,1:3])
    colnames(df) <- c('mz','intensity','rt')
    df.2 <- sqldf('select mz, max(intensity), rt from df group by mz, rt')
    colnames(df.2) <- c('mz','intensity','rt')
    truc <- spread(df.2, mz, intensity,fill=0)
    truc <- truc[1:cut.off,-1]
    truc <- cbind(matrix(rep(0,cut.off*99),nrow=cut.off,ncol=99),truc)
    ms.list[[j]] <- as.matrix(truc)
  }
  return(ms.list)
}

# mass2R <- function(files,type=c("Ramp", "pwiz", "netCDF"),cut.off=900){
#   ms.list <- list()
#   for(j in seq(length(files))){
#     aa <- openMSfile(files[j],backend=type)
#     ls.1 <- list()
#     for(i in seq(cut.off)){
#       ls.1[[i]] <- peaks(aa,i)
#     }
#     # get the min and max
#     mz.min <- round(runInfo(aa)[[2]])
#     mz.max <- round(runInfo(aa)[[3]])
#     # round the mz
#     ls.2 <- lapply(ls.1,function(x){
#       x[,1]<-round(x[,1])
# #       dupli <- duplicated(x[,1])
# #       dupli[1:100] <- F
# #       x <- x[!dupli,]
# #       colnames(x) <- c('mz','int')
# #       do.call("rbind",lapply(split(x, x[,'mz']), function(y){y[which.max(y[,2]),]}))
# #       x.agg <- aggregate(int ~ mz, x, max)
# #       x <- merge(x.agg, x)
#       return(x)
#       })
#     # change to vector
#     ls.3 <- lapply(ls.2,
#                    function(x){
#                      vec <- rep(0,mz.max)
#                      for(i in seq(nrow(x))){
#                        vec[x[i,1]]<-x[i,2]
#                      }
#                      return(vec)
#                    }
#     )
#     ms.list[[j]] <- t(as.matrix(as.data.frame(ls.3)))
#   }
#   return(ms.list)
# }
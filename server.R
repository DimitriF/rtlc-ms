
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(enviPick)
library(tidyr)
library(sqldf)
require("baseline")
require("prospectr")
library(plyr)
library(abind)
options(shiny.maxRequestSize=1000*1024^2)

source('mass2R.R')
source('f.plot.array.R')
source("PreProcess.function.R")

shinyServer(function(input, output) {
  file.Rdata <- reactive({
    validate(
      need(input$file.Rdata != "", "Please upload Rdata file")
    )
    load(input$file.Rdata$datapath)
    return(data)
  })
  Rdata.chrom <- reactive({
    file.Rdata()[[1]]
  })
  Rdata.batch <- reactive({
    file.Rdata()[[2]]
  })
  output$batch.simple <- renderDataTable({
    Rdata.batch()
  })
  
  Rdata.dim.vertical <- reactive({
    file.Rdata()[[3]]
  })
  
  
  file.ms <- reactive({
    validate(
      need(input$file.ms != "", "Please upload MS file(s)")
    )
    inFile <- input$file.ms
  })
  ms.array.0 <- reactive({
    inFile <- file.ms()
    mass2R(inFile[,4],input$cut.off)
  })
  ms.array.0.5 <- reactive({
    ms.array.0()[[as.numeric(input$sample.select)]]
  })
  ms.array.1 <- reactive({ # selection des mz d'interet
    data <- ms.array.0()
    for(i in seq(length(data))){
      data[[i]] <- data[[i]][,as.numeric(input$integrate.mz.choice)]
    }
    return(data)
  })
  ms.array.2 <- reactive({ # apply le preprocess
    lapply(ms.array.1(),function(x){
      t(f.preprocess(t(x),preprocess.order=input$integrate.mz.preprocessing.order,
                     preprocess.option=Preprocess.options(),
           training.data=t(x)))}
    )
  })
  output$ms.array.3 <- renderTable({
    ms.array.3()
  })
  ms.array.3 <- reactive({ # apply sum on each chromato
    data <- as.data.frame(lapply(ms.array.2(),function(x){apply(x,2,sum)}))
    colnames(data) <- file.ms()[,1]
    data <- t(data)
    colnames(data) <- input$integrate.mz.choice
    data
  })

  
  output$sample.select <- renderUI({
    choices=seq(length(file.ms()[,1]))
    names(choices) <- file.ms()[,1]
    selectizeInput('sample.select','Sample to explore',choices=choices)
  })

  
  #   output$myWebGL <- renderWebGL({
  # #     df <- ms.array.0.5()[seq(1,input$cut.off,by=3),100:400]
  # #     persp3d(y = seq(from=1,to=input$cut.off,len=ncol(df)),
  # #             x = seq(100, 400, len = nrow(df)),
  # #             z=df,xlab="mz",ylab="spectrum number",zlab="",
  # #             col='grey',back="lines")
  #     points3d(1:10, 1:10, 1:10)
  #     axes3d()
  #   })
  output$Inter.f.plot.array <- renderPlot({
    data <- Rdata.chrom()
    n.band <- as.numeric(input$sample.select)
    f.plot.array(data,n.band,file.ms()[,1],Rdata.dim.vertical()[2],Rdata.dim.vertical()[3],Rdata.dim.vertical()[4],inverse=F)
  })
  
  output$Inter.tic <- renderPlot({
    df <- ms.array.0.5()
    plot(apply(df,1,sum),type='l',main='TIC \nClick to view one spectrum')
  })
  output$Inter.eim <- renderPlot({
    df <- ms.array.0.5()
    brush <- input$Inter.ms_brush
    if(is.null(brush)){x1 <- 0;x2 <- 700}else{x1 <- round(as.numeric(brush[1]));x2<-round(as.numeric(brush[2]))}
    plot(apply(df[,x1:x2],1,sum),type='l',main=paste0('EIC: m/z = ', x1,' to ',x2))
  })
  output$Inter.sim <- renderPlot({
    df <- ms.array.0.5()
    click <- input$Inter.ms_click
    click.tic <- input$Inter.tic_click
    if(is.null(click.tic)){x1 <- 200}else{
      if(is.null(click)){x1 <- which.max(df[round(as.numeric(click.tic[1])),])}else{x1 <- round(as.numeric(click[1]))}
    }
    plot(df[,x1],type='l',main=paste0('SIM chromatogram: m/z = ', x1,'\nAt Tic click, change to spectrum max'))
  })
  output$Inter.ms.moyen <- renderPlot({
    df <- ms.array.0.5()
    plot(apply(df,2,sum),type='h' ,main='Mass Spectrum average')
  })
  output$Inter.ms <- renderPlot({
    df <- ms.array.0.5()
    click <- input$Inter.tic_click
    if(is.null(click)){x1 <- 1}else{x1 <- round(as.numeric(click[1]))}
    plot(df[x1,],type='h' ,main=paste0('Mass Spectrum for the spectrum ',x1,'\nBrush to zoom and produce Extracted Ion Chromatogram'),xlim=c(0,700))
  })
  output$Inter.ms.zoom <- renderPlot({
    df <- ms.array.0.5()
    click <- input$Inter.tic_click
    if(is.null(click)){x1 <- 1}else{x1 <- round(as.numeric(click[1]))}
    brush <- input$Inter.ms_brush
    if(is.null(brush)){zoom1 <- 0;zoom2 <- 700}else{zoom1 <- round(as.numeric(brush[1]));zoom2<-round(as.numeric(brush[2]))}
    plot(df[x1,],type='h'  ,main=paste0('Zomm Mass Spectrum for the spectrum ',x1,'\nClick to produce SIM chromatograme'),xlim=c(zoom1,zoom2))
  })
  
  output$download.csv.zip <- downloadHandler(
    filename = "mass_export.zip",
    content = function(file) {
      fs <- c()
      fichiers <- file.ms()[,1]
      for(i in seq(length(fichiers))){
        path <- paste0(fichiers[i],'.csv')
        fs <- c(fs,path)
        write.table(ms.array.0()[[i]],file=path,row.names = F,col.names = F,sep=';')
      }
      tempFile <- tempfile(fileext = ".zip")
      zip(zipfile=tempFile, files=fs)
      file.rename(tempFile, file)
    },
    contentType = "application/zip"
  )
  
  Preprocess.options <- reactive({
    Smoothing <- list(window.size = input$window.size,poly.order=input$poly.order,diff.order=input$diff.order)
    if(input$baseline == "als"){Baseline <- list(method=input$baseline,lambda.1=input$lambda.1,p=input$p,maxit.1=input$maxit.1)}
    if(input$baseline == "fillPeaks"){Baseline <- list(method=input$baseline,lambda.2=input$lambda.2,hwi=input$hwi,it=input$it,int=input$int)}
    if(input$baseline == "irls"){Baseline <- list(method=input$baseline,lambda1=input$lambda1,lambda2=input$lambda2,maxit.2=input$maxit.2,wi=input$wi)}
    if(input$baseline == "lowpass"){Baseline <- list(method=input$baseline,steep=input$steep,half=input$half)}
    if(input$baseline == "medianWindow"){Baseline <- list(method=input$baseline,hwm=input$hwm,hws=input$hws,end=input$end)}
    if(input$baseline == "modpolyfit"){Baseline <- list(method=input$baseline,degree=input$degree,tol=input$tol,rep=input$rep)}
    if(input$baseline == "peakDetection"){Baseline <- list(method=input$baseline,left=input$left,right=input$right,lwin=input$lwin,rwin=input$rwin)}
    if(input$baseline == "rfBaseline"){Baseline <- list(method=input$baseline)}
    if(input$baseline == "rollingBall"){Baseline <- list(method=input$baseline,wm=input$wm,ws=input$ws)}
    return(list(Smoothing=Smoothing,Baseline.correction=Baseline))
  })
})

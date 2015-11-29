#### License ####
#Copyright (C) {2015}  {Fichou Dimitri} 
#{dimitrifichou@laposte.net}

#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 2 of the License, or
# any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along
#with this program; if not, write to the Free Software Foundation, Inc.,
#51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#### rMS ######

library(shiny)

shinyUI(navbarPage("rMS",
                   tabPanel("Data input",
                            fluidRow(
                                     wellPanel(
                                       tabsetPanel(
                                         tabPanel('Load',
                                                  p('Originaly, it was possible to choose the format but because we focus on DART data, only mzXML is supported'),
                                                  fileInput('file.ms', 'Upload mzXML files',multiple=T),
                                                  fileInput('file.Rdata','Upload a Rdata file from rTLC'),
                                                  numericInput('cut.off','Cut off to take for the spectrum',900),
                                                  p('We need here to add some informations about the dart acquisition for a "click on picture and see the spectrum",
                                                    I put the cut-off for the dart so there will only be 900 spectrum by default but we need to know what is the corresponfing Rf on the plate
                                                    '),
                                                  downloadButton('download.csv.zip','Save zip file with csv')
                                                  ),
                                         tabPanel("batch",
                                                  dataTableOutput('batch.simple')
                                         ),
                                         tabPanel('Interactive',
                                                  uiOutput('sample.select'),
                                                  column(4,
                                                         plotOutput('Inter.f.plot.array'),
                                                         plotOutput('Inter.tic',
                                                                    click = "Inter.tic_click")
                                                  ),
                                                  column(4,
                                                         plotOutput('Inter.ms',
                                                                    brush = "Inter.ms_brush"),
                                                         plotOutput('Inter.ms.zoom',
                                                                    click = "Inter.ms_click")
                                                  ),
                                                  column(4,
                                                         plotOutput('Inter.eim'),
                                                         plotOutput('Inter.sim')
                                                  )
                                         ),
                                         tabPanel('Integration Options',
                                                  column(4,
                                                         h3('mz to integrate'),
                                                         selectizeInput('integrate.mz.choice','Choice',choices=c(100:700),multiple=T),
                                                         h5('Note that it\'s not the area under the peaks that are returned but the sum of the whole SIM chromatograms after preprocessing'),
                                                         h3('Preprocessing'),
                                                         selectizeInput('integrate.mz.preprocessing.order','Order to preprocess',choices=c('Smoothing','Baseline.correction'),multiple=T)
                                                         ),
                                                  column(4,
                                                         h4("Smoothing"),
                                                         helpText(   a("Click Here for help with this smoothing feature",target="_blank",     
                                                                       href="http://www.inside-r.org/node/206625")
                                                         ),
                                                         helpText(   a("Wikipedia link",target="_blank",     
                                                                       href="https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter")
                                                         ),
                                                         numericInput("window.size","size of the windows",3,min=3,max=NA,step=2),
                                                         numericInput("poly.order","polynomial order",1),
                                                         numericInput("diff.order","differentiation order",0)
                                                  ),
                                                  column(4,
                                                         h4("Baseline"),
                                                         helpText(   a("Click Here for help with the Baseline feature",target="_blank",     
                                                                       href="http://cran.r-project.org/web/packages/baseline/baseline.pdf")
                                                         ),
                                                         selectizeInput("baseline", "type of baseline", choices=c("als","fillPeaks","irls","lowpass","medianWindow","modpolyfit","peakDetection","rfbaseline","rollingBall"),select=NULL),
                                                         conditionalPanel(condition="input.baseline=='als'",
                                                                          numericInput("lambda.1","lambda : 2nd derivative constraint",5),
                                                                          numericInput("p","p : weighting of positive residuals",0.05),
                                                                          numericInput("maxit.1","maxit : maximum number of iterations",20)
                                                         ),
                                                         conditionalPanel(condition="input.baseline=='fillPeaks'",
                                                                          numericInput("lambda.2","lambda : 2nd derivative constraint for primary smoothing",6),
                                                                          numericInput("hwi","hwi : half width of local windows",100),
                                                                          numericInput("it","it : number of iterations in suppression loop",10),
                                                                          numericInput("int","int : number of buckets to divide spectra into",200)
                                                         ),
                                                         conditionalPanel(condition="input.baseline=='irls'",
                                                                          numericInput("lambda1","lambda1 : 2nd derivative constraint for primary smoothing",5),
                                                                          numericInput("lambda2","lambda2 : 2nd derivative constraint for secondary smoothing",9),
                                                                          numericInput("maxit.2","maxit : maximum number of iterations",200),
                                                                          numericInput("wi","wi : weighting of positive residuals",0.05)
                                                         ),
                                                         conditionalPanel(condition="input.baseline=='lowpass'",
                                                                          numericInput("steep","steep : steepness of filter curve",2),
                                                                          numericInput("half","half : half way point of filter curve",5)
                                                         ),
                                                         conditionalPanel(condition="input.baseline=='medianWindow'",
                                                                          numericInput("hwm","hwm : window half width for local medians",300),
                                                                          numericInput("hws","hws : window half width for local smoothing",5),
                                                                          checkboxInput("end","end : original endpoint handling",F)
                                                         ),
                                                         conditionalPanel(condition="input.baseline=='modpolyfit'",
                                                                          numericInput("degree","degree : degree of polynomial",4),
                                                                          numericInput("tol","tol : tolerance of difference between iterations",0.001),
                                                                          numericInput("rep","rep : maximum number of iterations",100)
                                                         ),
                                                         conditionalPanel(condition="input.baseline=='peakDetection'",
                                                                          numericInput("left","left : smallest window size for peak widths",30),
                                                                          numericInput("right","right : largest window size for peak widths",300),
                                                                          numericInput("lwin","lwin : Smallest window size for minimums and medians in peak removed spectra",50),
                                                                          numericInput("rwin","rwin : Largest window size for minimums and medians in peak removed spectra",50),
                                                                          numericInput("snminimum","snminimum : Minimum signal to noise ratio for accepting peaks",10)
                                                         ),
                                                         conditionalPanel(condition="input.baseline=='rollingBall'",
                                                                          numericInput("wm","wm : Width of local window for minimization/maximization",200),
                                                                          numericInput("ws","ws : Width of local window for smoothing",200)
                                                         )
                                                         )
                                                  ),
                                         tabPanel('Data integrated',
                                                  tableOutput('ms.array.3')
                                                  )
                                       )
                                     )
                              )
                   ),
                   tags$head(tags$style(type="text/css", "tfoot {display: table-header-group}")),
                   tags$head(tags$style(HTML(".shiny-output-error-validation {color: red;font-size: 30px}"))),
                   tags$head(tags$style(type="text/css", ".shiny-progress .progress {position: absolute;width: 100%;top: 100px;height: 10px;margin: 0px;}")),
                   tags$head(tags$style(type="text/css", ".shiny-progress .progress-text {position: absolute;border-style: solid;
                                        border-width: 2px;right: 10px;height: 36px;width: 50%;background-color: #EEF8FF;margin: 0px;padding: 2px 3px;opacity: 1;}"))
                   
                   ))

################################################################################
#Libraries
################################################################################

library(shiny)
library(readxl)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(EnhancedVolcano)
library(Hmisc)
library(plotly)
library(biomaRt)
library(shinycssloaders)
library(BiocManager)

options(repos = BiocManager::repositories())
options(shiny.maxRequestSize = 30*1024^2)

## Versions

#v.2.0 (September 2023): Input (xlsx, csv, and txt)
                        #Update on geometric mean for duplicated values
                        #Users can choose groups to be plotted in Boxplots

################################################################################
#Required functions
################################################################################

#Normfinder 

Normfinder=function(dat0,Groups=TRUE,ctVal=TRUE,pStabLim=0.25){
  #
  #
  ntotal=dim(dat0)[2] # number of samples
  k0=dim(dat0)[1] # number of rows
  #
  if (Groups){
    ngenes=k0-1 # number of genes
    genenames=rownames(dat0)[-k0]
    grId=dat0[k0,]
    dat0=dat0[-k0,]
  } else {
    ngenes=k0 # number of genes
    genenames=rownames(dat0)
    grId=rep(1,ntotal)
  }
  #
  dat=matrix(as.numeric(unlist(dat0)),ngenes,ntotal) # matrix instead of list
  #
  if (!ctVal){dat=log2(dat)} # transform to log2 values
  #
  samplenames=colnames(dat0)
  grId=factor(unlist(grId))  # group identifier
  groupnames=levels(grId)  # group names
  ngr=length(levels(grId)) # number of groups
  # Number of samples in each group:
  nsamples=rep(0,ngr)
  for (group in 1:ngr){nsamples[group]=sum(grId==groupnames[group])}
  #
  #
  MakeStab=function(da){
    ngenes=dim(da)[1]
    # Sample averages
    sampleavg=apply(da,2,mean)
    # Gene averages within group
    genegroupavg=matrix(0,ngenes,ngr)
    for (group in 1:ngr){
      genegroupavg[,group]=apply(da[,grId==groupnames[group]],1,mean)}
    # Group averages
    groupavg=rep(0,ngr)
    for (group in 1:ngr){groupavg[group]=mean(da[,grId==groupnames[group]])}
    
    # Variances 
    GGvar=matrix(0,ngenes,ngr)
    for (group in 1:ngr){
      grset=(grId==groupnames[group])
      a=rep(0,ngenes)
      for (gene in 1:ngenes){
        a[gene]=sum((da[gene,grset]-genegroupavg[gene,group]-
                       sampleavg[grset]+groupavg[group])^2)/(nsamples[group]-1)
      }
      GGvar[,group]=(a-sum(a)/(ngenes*ngenes-ngenes))/(1-2/ngenes)
    }
    #
    # Change possible negative values
    genegroupMinvar=matrix(0,ngenes,ngr)
    for (group in 1:ngr){
      grset=(grId==groupnames[group])
      z=da[,grset]
      for (gene in 1:ngenes){
        varpair=rep(0,ngenes)
        for (gene1 in 1:ngenes){varpair[gene1]=var(z[gene,]-z[gene1,])}
        genegroupMinvar[gene,group]=min(varpair[-gene])/4
      }
    }
    #
    # Final variances
    GGvar=ifelse(GGvar<0,genegroupMinvar,GGvar)
    #
    # Old stability measure for each gene is calculated:
    #
    dif=genegroupavg
    difgeneavg=apply(dif,1,mean)
    difgroupavg=apply(dif,2,mean)
    difavg=mean(dif)
    for (gene in 1:ngenes){
      for (group in 1:ngr){
        dif[gene,group]=dif[gene,group]-difgeneavg[gene]-difgroupavg[group]+difavg
      }
    }
    #
    nsampMatrix=matrix(rep(nsamples,ngenes),ngenes,ngr,byrow=T)
    vardif=GGvar/nsampMatrix
    gamma=sum(dif*dif)/((ngr-1)*(ngenes-1))-sum(vardif)/(ngenes*ngr)
    gamma=ifelse(gamma<0,0,gamma)
    #
    difnew=dif*gamma/(gamma+vardif)
    varnew=vardif+gamma*vardif/(gamma+vardif)
    Ostab0=abs(difnew)+sqrt(varnew)
    Ostab=apply(Ostab0,1,mean)
    #
    # Measure of group differences:
    mud=rep(0,ngenes)
    for (gene in 1:ngenes){
      mud[gene]=2*max(abs(dif[gene,]))
    }
    # Common variance:
    genevar=rep(0,ngenes)
    for (gene in 1:ngenes){
      genevar[gene]=sum((nsamples-1)*GGvar[gene,])/(sum(nsamples)-ngr)
    }
    Gsd=sqrt(genevar)
    #
    # Return results:
    #
    return(cbind(mud,Gsd,Ostab,rep(gamma,ngenes),GGvar,dif))
  }    # End of function MakeStab
  #
  #
  MakeComb2=function(g1,g2,res){
    gam=res[1,4]
    d1=res[g1,(4+ngr+1):(4+ngr+ngr)]; d2=res[g2,(4+ngr+1):(4+ngr+ngr)]
    s1=res[g1,(4+1):(4+ngr)]; s2=res[g2,(4+1):(4+ngr)]
    rho=abs(gam*d1/(gam+s1/nsamples)+gam*d2/(gam+s2/nsamples))*
      sqrt(ngenes/(ngenes-2))/2
    rho=rho+sqrt(s1/nsamples+gam*s1/(nsamples*gam+s1)+
                   s2/nsamples+gam*s2/(nsamples*gam+s2))/2
    return(mean(rho))
  }
  #
  #
  MakeStabOne=function(da){
    ngenes=dim(da)[1]
    # Sample averages
    sampleavg=apply(da,2,mean)
    # Gene averages
    geneavg=apply(da,1,mean)
    totalavg=mean(da)
    #
    # Variances 
    genevar0=rep(0,ngenes)
    for (gene in 1:ngenes){
      genevar0[gene]=
        sum((dat[gene,]-geneavg[gene]-sampleavg+totalavg)^2)/
        ((ntotal-1)*(1-2/ngenes))
    }
    genevar=genevar0-sum(genevar0)/(ngenes*ngenes-ngenes)
    #
    # Change possible negative values
    geneMinvar=rep(0,ngenes)
    z=da
    for (gene in 1:ngenes){
      varpair=rep(0,ngenes)
      for (gene1 in 1:ngenes){varpair[gene1]=var(z[gene,]-z[gene1,])}
      geneMinvar[gene]=min(varpair[-gene])/4
    }
    # Final variances
    genevar=ifelse(genevar<0,geneMinvar,genevar)
    #
    return(genevar)
  }
  #     End of function MakeStabOne
  #
  #################################################
  #
  # Main part
  #
  if (ngr>1){   # More than one group.
    #
    res=MakeStab(dat)
    #
    gcand=c(1:ngenes)[res[,3]<pStabLim]
    ncand=length(gcand)
    if (ncand<4){
      if (ngenes>3){
        li=sort(res[,3])[4]
        gcand=c(1:ngenes)[res[,3]<=li]
        ncand=length(gcand)
      } else {
        gcand=c(1:ngenes)
        ncand=length(gcand)
      }
    }
    #
    vv2=c()
    #
    for (g1 in 1:(ncand-1)){
      for (g2 in (g1+1):ncand){
        qmeas=MakeComb2(gcand[g1],gcand[g2],res)
        vv2=rbind(vv2,c(gcand[g1],gcand[g2],qmeas))
      }}
    #
    ord=order(res[,3])
    FinalRes=list(Ordered=
                    data.frame("GroupDif"=round(res[ord,1],2),"GroupSD"=round(res[ord,2],2),
                               "Stability"=round(res[ord,3],2),row.names=genenames[ord]),
                  UnOrdered=
                    data.frame("GroupDif"=round(res[,1],2),"GroupSD"=round(res[,2],2),
                               "Stability"=round(res[,3],2),
                               "IGroupSD"=round(sqrt(res[,(4+1):(4+ngr)]),2),
                               "IGroupDif"=round(res[,(4+ngr+1):(4+ngr+ngr)],2),
                               row.names=genenames),
                  PairOfGenes=
                    data.frame("Gene1"=genenames[vv2[,1]],"Gene2"=genenames[vv2[,2]],
                               "Stability"=round(vv2[,3],2)))
    #
    return(FinalRes)
    #
  } else {    # End of more than one group: next is for one group only.
    #
    #
    sigma=sqrt(MakeStabOne(dat))
    #
    siglim=(min(sigma)+0.1)
    gcand=c(1:ngenes)[sigma<siglim]
    ncand=length(gcand)
    #
    if ((ncand>=2)&(ngenes>3)){
      #
      vv2=c()
      #
      for (g1 in 1:(ncand-1)){
        for (g2 in (g1+1):ncand){
          dat1=rbind(dat[-c(gcand[g1],gcand[g2]),],
                     apply(dat[c(gcand[g1],gcand[g2]),],2,mean))
          qmeas=sqrt(MakeStabOne(dat1))
          vv2=rbind(vv2,c(gcand[g1],gcand[g2],qmeas[ngenes-1]))
        }}
      ord=order(sigma)
      FinalRes=list(Ordered=
                      data.frame("GroupSD"=round(sigma[ord],2),row.names=genenames[ord]),
                    PairOfGenes=
                      data.frame("Gene1"=genenames[vv2[,1]],"Gene2"=genenames[vv2[,2]],
                                 "GroupSD"=round(vv2[,3],2)))
    } else { # No combined genes to consider
      ord=order(sigma)
      FinalRes=list(Ordered=
                      data.frame("GroupSD"=round(sigma[ord],2),row.names=genenames[ord]))
    } # End ncand<2 or ngenes<=3
    #
    return(FinalRes)
    #
  }  # End one group only
  #
} # End of main function

# Version history:
# Version dated 17/12-2013.
# Version dated 18/06-2014: Correction in line 121: s2 is squared.
# Version dated 23/08-2014: Correction to above correction in line 121: 
#   s1 and s2 are not squared. 
#   We thank John R. Young for pointing to the error in line 121.
# Version dated 05/01-2015: sum()/2 changed to mean() in line 126.
#   This misprint appears in the Supplementary Information as well.
#   We thank Mark Birtwistle for pointing out this error.

#Color verifying

verificando_cores<-function(x){
  try(plotrix::color.id(x)[1],silent=T)->x
  if(x %in% colors()){
    TRUE
  }else{
    FALSE
  }
}

#Correlat table
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  out<-data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
  return(out)
}

################################################################################
#UIs
################################################################################

Data <- tabPanel(
  title = HTML(paste("","Input",sep = "<br/>")),
  titlePanel("Data Input"),
  sidebarLayout(
    sidebarPanel(
      title = "Inputs",
      selectInput("data_format","Choose the format of your data",
                  c("One column per gene","One column with gene names")),  
      fileInput("xlsx_input","Upload the file with Ct data (xls, xlsx, txt, or csv formats are allowed)",
                accept=c(".xlsx",".txt",".csv")),
      actionButton("run_button","Start!",icon=icon("play")),
      br(),
      br(),
      uiOutput("id0"),
      uiOutput("id3"),
      uiOutput("id4"),
      uiOutput("id5"),
      br(),
      br(),
      uiOutput("idgene"),
      uiOutput("id6"),
      uiOutput("id7"),
      uiOutput("id8"),
      br(),
      br(),
      uiOutput("id9")
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "Data verification",
          uiOutput("id2"),
          textOutput("verif1"),
          textOutput("verif2"),
          tags$head(tags$style("#verif1{color: red;
                                 font-size: 18px;
                                 font-style: italic;
                                 }"
          )
          )
        #   ),
        # tabPanel(
        #   title="Color Palette",
        #   plotOutput("palette_plot1",height = "150px"),
        #   htmlOutput("palette_text1")
        )
      )
    )
  )
)

outliers_normal<-tabPanel(
  title=HTML(paste("    Data    ","Normalization",sep = "<br/>")),
  titlePanel("Data Normalization"),
  sidebarLayout(
    sidebarPanel(
      title="Inputs",
      # selectInput("gene_norm","Select one gene from the list",c(NULL)),
      # checkboxInput("out_norm","Do you want to remove outliers? Only for this plot. Then you can choose if remove outliers for further analysis",FALSE),
      # actionButton("run_button_Norm","Plot this gene",icon=icon("play")),
      # br(),
      # br(),
      selectInput("house_analysis","Do you want to indicate housekeeping gene candidates?",c("No","Yes")),
      uiOutput("out_normal1"),
      actionButton("house_Norm","Evaluate housekeeping gene(s)",icon=icon("edge")),
      br(),
      br(),
      uiOutput("out_normal2"),
      uiOutput("out_normal3"),
      uiOutput("out_normal4"),
      br(),
      uiOutput("out_normal5"),
      uiOutput("out_normal6"),width=3
      
    ),
    mainPanel(
      tabsetPanel(
        # tabPanel(
        # title="-??Ct plot",
        # plotOutput("norm_plot",height = "500px"),
        # textOutput("norm_text"),
        # tags$head(tags$style("#norm_text{color: blue;
        #                          font-size: 18px;
        #                          font-style: italic;
        #                          }"
        # )
        # )
        # ),
        tabPanel(
          title="Housekeeping genes",
          textOutput("house_verif1"),
          tableOutput("house_table"),
          tableOutput("house_table2"),
          textOutput("norm_text1"),
          tags$head(tags$style("#norm_text1{color: blue;
                                 font-size: 18px;
                                 font-style: italic;
                                 }"
          )),
          textOutput("norm_text2"),
          tags$head(tags$style("#norm_text2{color: red;
                                 font-size: 18px;
                                 font-style: italic;
                                 }"
          )
          ),
          uiOutput("out_resul1")
        )
      )
    )
  )
)
  
  
  
  
  
expressionAna <- tabPanel(
  title = HTML(paste("Expression","Analysis",sep = "<br/>")),
  titlePanel("Gene Expression"),
  sidebarLayout(
    sidebarPanel(
      title = "Inputs",
      selectInput("ref.group0","Select the variable to be compared",c(NULL)),
      checkboxGroupInput("ref.group0_5","Select groups to be compared",c(NULL)),
      selectInput("ref.group","Select the reference group (for statistical comparisons)",c(NULL)),  
      selectInput("gene","Select one gene from the list",c(NULL)),
      # checkboxInput("specific_expr","Do you want to focus on a subcohort based on gene expression?"),
      uiOutput("UIinputexpr1"),
      uiOutput("UIinputexpr2"),
      selectInput("method_exp","Choose the test to compare your data:",c("Depending on their distribution (see Normality tab)",
                                                                         "Non-parametric test (Mann-Whitney)",
                                                                         "Parametric test (T-test)")),
      actionButton("run_button_Exp","Run analysis",icon=icon("play")),width=3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "-\u25B3Ct Plot",
          plotOutput("express2Plot",height = "500px"),
          uiOutput("UIexpr1")
        ),
        tabPanel(
          title = "2^-\u25B3\u25B3Ct Plot",
          plotOutput("expressPlot",height = "500px"),
          uiOutput("UIexpr2")
        ),
        tabPanel(
          title = "Normality Analysis",
          tableOutput("NormTable")
        ),
        tabPanel(
          title = "Results",
          tableOutput("expressTable"),
          uiOutput("idExp1")
        )
        
      )
    )
  )
)

volcano_plots <-tabPanel(
  title=HTML(paste("Volcano","Plots",sep = "<br/>")),
  titlePanel("Volcano plots"),
  sidebarLayout(
    sidebarPanel(
      selectInput("group0","Select the variable to be compared",c(NULL)),
      selectInput("group1","Select the first group for comparison",c(NULL)),
      selectInput("group2","Select the second group for comparison",c(NULL)),
      selectInput("method_volc","Choose the test to compare your data:",c("Non-parametric test (Mann-Whitney)",
                                                                         "Parametric test (T-test)")),
      selectInput("p_plotted","Choose the p-value for the plot",c("adjusted p-value (BH)","non adjusted p-value")),
      textInput("p_val_volc","Please set the p value cutoff","1e-4"),
      textInput("FC_volc","Please set the Fold change cutoff","2"),
      # checkboxInput("specific_volc","Do you want to focus on a subcohort based on gene expression?"),
      uiOutput("UIinputvolc1"),
      uiOutput("UIinputvolc2"),
      actionButton("run_button_volc","Run analysis",icon=icon("play")),width=3),
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "Graphics",
          plotOutput("volcanoPlot",height = "600px"),
          uiOutput("UIvolcano1")
        ),
        tabPanel(
          title = "Table",
          uiOutput("UIvolcano2"),
          tableOutput("volcanoTable")
          
        )
      )
    )
  )
)

correlat_plots <- tabPanel(
  title = HTML(paste("Correlation","Analysis",sep = "<br/>")), 
  titlePanel("Correlation"),
  sidebarLayout(
    sidebarPanel(
      selectInput("Treatment0","Select the variable to be compared",c(NULL)),
      selectInput("Treatment2","Choose an specific group",c(NULL)),
      numericInput("valor_p","Choose a value for alpha for correlation evaluation (1%:p<0.01)",value=1),
      actionButton("run_button_correlat","Run analysis",icon=icon("play")),width=3),
    mainPanel(
      plotlyOutput("correlatPlot")
    )
  )   
)

scatter_plots <- tabPanel(
  title = HTML(paste("Scatter","Plots",sep = "<br/>")),
  titlePanel("Scatter"),
  sidebarLayout(
    sidebarPanel(
      selectInput("sc_group0","Select the variable to be compared",c(NULL)),
      selectInput("sc_group1","Choose the first group to be plotted",c(NULL)),
      selectInput("sc_group2","Choose the second group to be plotted",c(NULL)),
      selectInput("sc_gene1","Select the first gene for correlation analysis",c( NULL)),
      selectInput("sc_gene2","Select the second gene for correlation analysis",c(NULL)),
      actionButton("run_button_sc","Run analysis",icon=icon("play")),width=3),
    mainPanel(
      plotOutput("scatterPlot"),
      uiOutput("UIscatter1")
    )
  )   
)

glossary<-tabPanel(
  title = HTML(paste("","Glossary",sep = "<br/>")),
  titlePanel("Glossary"),
  sidebarLayout(
    sidebarPanel(
      selectInput("gloss_specie","Select the specie you are working with",c(NULL)),
      selectInput("gloss_type","Choose the format of your gene list",c(NULL)),
      actionButton("run_gloss","Find genes!",icon=icon("play")),width=3),
  mainPanel(
    uiOutput("gloss_id0"),
    dataTableOutput("gloss")
  ))
)


user_guide<-tabPanel(
  title=HTML(paste("","User Guide",sep = "<br/>")),
  titlePanel("Guide"),
  sidebarLayout(
    sidebarPanel(
      p("You can check our guide for the deltaXpress (\u25B3Xpress) app at "),
      a("this link", 
        href = "https://drive.google.com/file/d/1RfwtaTxasq_WgTXkZW7QXVB4SIgYFOca/view?usp=share_link",
        target="_blank"),
      br(),
      br(),
      h4("Contact"),
      p("In case of comments or suggestions, please send an e-mail to Alexis Murillo (agmurilloc@usp.br) or Roger Chammas (rchammas@usp.br)")
    ),
    mainPanel()
  )
)


ui <- navbarPage(
  tags$head(
    tags$style(HTML('.navbar-nav > li > a, .navbar-brand {
                            padding-top:10 px !important; 
                            padding-bottom:2 px !important;
                            height: 60px;
                            }
                           .navbar {min-height:60px !important;}
                    '))
  ),
  title = HTML(paste("deltaXpress: A tool for mapping", "differentially correlated genes",sep = "<br/>")),  #&#916
  Data,
  outliers_normal,
  expressionAna,
  volcano_plots,
  correlat_plots,
  scatter_plots,
  glossary,
  user_guide
)


################################################################################
#Server
################################################################################

server <- function(input, output,session){
  
  datos <- reactiveValues(x = NULL,ref=NULL,order_color=NULL,
                          final=NULL,o_var=NULL)
  
  observe({
    req(input$xlsx_input)

    
    tail(unlist(strsplit(input$xlsx_input$datapath,"[.]")),1)->formato
    if(formato=="xlsx"){
      x <-  read_excel(input$xlsx_input$datapath)
    }else if(formato=="txt"){
      x <- read.table(input$xlsx_input$datapath,sep="\t",
                      header = TRUE, 
                      stringsAsFactors = FALSE)
      
    }else if (formato=="csv"){
      x <- read.csv(input$xlsx_input$datapath)
    }
    
    gm_mean<-function(w){
      w<-na.omit(w)
      exp(mean(log(as.numeric(w))))
    }
    
    
    tryCatch(
      {
        
        
        if(input$data_format=="One column per gene"){
                

                if(length(class(as.data.frame(x)[,3]))==1 && class(as.data.frame(x)[,3])=="numeric"){
                
                  colnames(x)<-c("Sample","Group",colnames(x[,-c(1:2)]))
                
                  aggregate(x=as.data.frame(x),
                          by=list(x$Sample,
                                  x$Group),
                          FUN=gm_mean )->x
                
                x[,-c(3:4)]->x
                
                colnames(x)<-c("Sample","Group",colnames(x[,-c(1:2)]))
                
                x<-x%>%pivot_longer(!c("Sample","Group"),names_to="Name",values_to="Value")
                
                renderText({"File correctly loaded. Please, click on Start!"})->output$verif1
                renderTable({NULL})->output$verif
                renderText({""})->output$verif2
                
                } else {
                  x<-NULL
                  renderText({"Please, verify the format of the input file."})->output$verif1
                  renderTable({NULL})->output$verif
                  renderText({""})->output$verif2
                }

                
                
        } else {
          
          
          if(ncol(as.data.frame(x))==4){
            
            colnames(x)<-c("Sample","Group","Name","Value")
            
            x<-x%>%group_by(Sample,Group,Name)%>%
              mutate(Value=gm_mean(Value))%>%distinct()
            
            renderText({"File correctly loaded. Please, click on Start!"})->output$verif1
            renderTable({NULL})->output$verif
            renderText({""})->output$verif2
            
            }else{
            
            x<-NULL
            renderText({"Please, verify the format of the input file."})->output$verif1
            renderTable({NULL})->output$verif
            renderText({""})->output$verif2
          }
          
          # if(length(excel_sheets(input$xlsx_input$datapath))>1){
          #   order_color<-read_excel(input$xlsx_input$datapath,sheet = 2)
          #   grupos<-as.matrix(unique(x[,2]))
          #   if(length(grupos[grupos %in% as.character(as.matrix(order_color[1:nrow(grupos),1]))])==length(grupos) &
          #      is.numeric(as.matrix(order_color[1:length(grupos),2]))){
          #     if(length(unique(sapply(as.matrix(order_color[1:length(grupos),3]),verificando_cores)))==1 &
          #        unique(sapply(as.matrix(order_color[1:length(grupos),3]),verificando_cores))== TRUE){
          #       
          #       datos$order_color<-order_color[1:length(grupos),]
          #       
          #     }
          #   }
          #   
          # }
        }
        
        },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      })
    # Reactive values updated from x
    
    
    datos$x <- as.data.frame(x)
    })

  observeEvent(input$run_button,{
    
    req(datos$x)

    output$id2 <-
        renderUI({
          shinycssloaders::withSpinner(
            tableOutput("verif")
          )
        })
      
    if(is.null(datos$order_color)){
      hue_pal()(length(unique(datos$x[,2])))->col.palette
      as.matrix(unique(datos$x[order(datos$x[,2]),][,2]))->lista_grupos
      
    }else{
      as.character(as.matrix(datos$order_color[order(datos$order_color[,2]),][,3]))->col.palette
      as.matrix(datos$order_color[order(datos$order_color[,2]),][,1])->lista_grupos
    }

    data.frame(Parameters=c("Sample information (column name)",
                            "Group information (column name)",
                            "Number of samples",
                            "Number of groups",
                            "Groups of Samples",
                            "Number of genes",
                            "List of genes"),
               Values=c(colnames(datos$x[,1:2])[1],
                        colnames(datos$x[,1:2])[2],
                        length(unique(datos$x[,1])),
                        length(unique(datos$x[,2])),
                        paste0(t(unique(datos$x[order(datos$x[,2]),][,2])),collapse=", "),
                        length(unique(datos$x[,3])),
                        gsub("(.{35})","\\1\n\\2",paste0(t(unique(datos$x[order(datos$x[,3]),][,3])),collapse=", ")))
               )->resumen
   
    
    renderTable({resumen})->output$verif
    renderText({""})->output$verif1
    renderText({"If the above results agree with the expected ones, continue to the next section. If not, revise the uploaded file."})->output$verif2
    renderPlot({show_col(col.palette,labels=F,ncol=length(unique(datos$x[,2])))},height = 400)->output$palette_plot1
    renderUI({HTML(paste0("From left to right: <br/>",
                            paste0(as.character(lista_grupos),collapse=" <br/> ")))})->output$palette_text1
    
    output$id0 <-
      renderUI({checkboxInput("change_group","Do you want to change the name of some group?",value=FALSE)
      })
    
    output$idgene <-
      renderUI({checkboxInput("change_gene","Do you want to change the name of some gene?",value=FALSE)
      })
    
    output$id9 <-
      renderUI({
        actionButton("go_next","Go to the next section",icon=icon("check",lib = "font-awesome"))
      })
    
    })
  
  
  output$id3 <-
    renderUI({
      req(input$change_group)
      selectInput("group_to_change","Select the group",paste0(t(unique(datos$x[order(datos$x[,2]),][,2]))))
    })
  
  output$id4 <-
    renderUI({
      req(input$change_group)
      textInput("new_group_name","Please digit the new name for the above group")
    })
  
  output$id5 <-
    renderUI({
      req(input$change_group)
      actionButton("change_the_name","Change the group name",icon=icon("right-left",lib = "font-awesome"))
    })
  
  output$id6 <-
    renderUI({
      req(input$change_gene)
      selectInput("gene_to_change","Select the gene",paste0(t(unique(datos$x[order(datos$x[,3]),][,3]))))
    })
  
  output$id7 <-
    renderUI({
      req(input$change_gene)
      textInput("new_gene_name","Please digit the new name for the above gene")
    })
  
  output$id8 <-
    renderUI({
      req(input$change_gene)
      actionButton("change_the_gene_name","Change the gene name",icon=icon("right-left",lib = "font-awesome"))
    })
  
  
  output$out_normal1<-
    renderUI({
      req(input$house_analysis=="Yes")
      textInput("house_genes","Please input the housekeeping genes of your experiment (more than 2). Separate genes by space. i.e. ACTB GAPDH B2M")
    })
  
  observeEvent(input$change_the_name,{
    datos$x[which(datos$x[,2]==input$group_to_change),2]<-input$new_group_name
    
    
    if(is.null(datos$order_color)){
      hue_pal()(length(unique(datos$x[,2])))->col.palette
      as.matrix(unique(datos$x[order(datos$x[,2]),][,2]))->lista_grupos
    #   
    # }else{
    #   datos$order_color[which(datos$order_color[,1]==input$group_to_change),1]<-input$new_group_name
    #   as.character(as.matrix(datos$order_color[order(datos$order_color[,2]),][,3]))->col.palette
    #   as.matrix(datos$order_color[order(datos$order_color[,2]),][,1])->lista_grupos
    }
    
    data.frame(Parameters=c("Sample information (column name)",
                            "Group information (column name)",
                            "Number of samples",
                            "Number of groups",
                            "Groups of Samples",
                            "Number of genes",
                            "List of genes"),
               Values=c(colnames(datos$x[,1:2])[1],
                        colnames(datos$x[,1:2])[2],
                        length(unique(datos$x[,1])),
                        length(unique(datos$x[,2])),
                        paste0(t(unique(datos$x[order(datos$x[,2]),][,2])),collapse=", "),
                        length(unique(datos$x[,3])),
                        gsub("(.{35})","\\1\n\\2",paste0(t(unique(datos$x[order(datos$x[,3]),][,3])),collapse=", ")))
    )->resumen
    
    
    renderTable({resumen})->output$verif
    renderText({""})->output$verif1
    renderText({"If the above results agree with the expected ones, continue to the next section. If not, revise the uploaded file."})->output$verif2
    renderPlot({show_col(col.palette,labels=F,ncol=length(unique(datos$x[,2])))},height = 400)->output$palette_plot1
    renderUI({HTML(paste0("From left to right: <br/>",
                          paste0(as.character(lista_grupos),collapse=" <br/> ")))})->output$palette_text1
    
    output$id3 <-
      renderUI({
        req(input$change_group)
        selectInput("group_to_change","Select the group",paste0(t(unique(datos$x[order(datos$x[,2]),2]))))
      })
    
    
    })
  
  observeEvent(input$change_the_gene_name,{
    

    
    datos$x[which(datos$x[,3]==input$gene_to_change),3]<-input$new_gene_name
    data.frame(Parameters=c("Sample information (column name)",
                            "Group information (column name)",
                            "Number of samples",
                            "Number of groups",
                            "Groups of Samples",
                            "Number of genes",
                            "List of genes"),
               Values=c(colnames(datos$x[,1:2])[1],
                        colnames(datos$x[,1:2])[2],
                        length(unique(datos$x[,1])),
                        length(unique(datos$x[,2])),
                        paste0(t(unique(datos$x[order(datos$x[,2]),][,2])),collapse=","),
                        length(unique(datos$x[,3])),
                        gsub("(.{35})","\\1\n\\2",paste0(t(unique(datos$x[order(datos$x[,3]),][,3])),collapse=", ")))
    )->resumen
    
    
    renderTable({resumen})->output$verif
    
    output$id6 <-
      renderUI({
        req(input$change_gene)
        selectInput("gene_to_change","Select the gene",paste0(t(unique(datos$x[,3]))))
      })
    
    
  })
  
  observeEvent(input$go_next,{
    
    updateSelectInput(session,"gene_norm","Select one gene from the list",
                      paste0(t(unique(datos$x[order(datos$x[,3]),][,3]))))
    
  })
  
  
  # observeEvent(input$run_button_Norm,{
  #   
  #   if(!is.null(input$gene_norm) & !is.null(datos$x)){
  #     input$gene_norm->gene_norm
  #     dados_norm<-datos$x%>%group_by(Name)%>%
  #       mutate(row=row_number())%>%
  #       tidyr::pivot_wider(names_from = Name, values_from = Value)
  #     dados_norm<-dados_norm[,-3]
  #     
  #     if(length(unique(table(dados_norm[,1])))>1){
  #       dados_norm<-datos$x%>%group_by(Name)%>%
  #         tidyr::pivot_wider(names_from = Name, values_from = Value)
  #     }
  #     
  # 
  #     dados_norm<-unnest(dados_norm,colnames(dados_norm[,3:ncol(dados_norm)]))
  #     
  #     unique(dados_norm[,2])->nome_grupo
  #     which(colnames(dados_norm)==gene_norm)->coluna_norm
  #     
  #     
  #     if(input$out_norm==TRUE){
  #       
  #         for(k in 1:nrow(nome_grupo)){
  #           
  #           quartiles <- quantile(dados_norm[which(dados_norm[,2]==as.character(nome_grupo[k,1])),coluna_norm], probs=c(.25, .75), na.rm = TRUE)
  #           IQR <- IQR(as.matrix(dados_norm[which(dados_norm[,2]==as.character(nome_grupo[k,1])),coluna_norm]), na.rm = TRUE)
  # 
  #           Lower <- quartiles[1] - 1.5*IQR
  #           Upper <- quartiles[2] + 1.5*IQR
  # 
  #           dados_norm[which(dados_norm[,2]==as.character(nome_grupo[k,1]) & dados_norm[,coluna_norm]<Lower),coluna_norm]<-NA
  #           dados_norm[which(dados_norm[,2]==as.character(nome_grupo[k,1]) & dados_norm[,coluna_norm]> Upper),coluna_norm]<-NA
  #         }
  #       
  #       
  #     }
  #     
  #     
  #     
  #     filt_norm<-dados_norm[which(!is.na(dados_norm[,coluna_norm])),c(2,coluna_norm)]
  #     
  #     p1_norm<-dados_norm[which(!is.na(dados_norm[,coluna_norm])),c(2,coluna_norm)]%>%
  #       ggplot(aes(Group,get(gene_norm)))+geom_boxplot()+geom_jitter(width=0.1,size=1,alpha=0.7)+
  #       theme_classic(base_size=18)+
  #       theme(axis.text.x = element_text(angle = 45,hjust=1))+
  #       labs(title=paste(gene_norm,"levels across groups"),x="Groups",y=bquote(italic(.(gene_norm)) ~ levels (Ct)))
  #     
  #     
  #     renderPlot({p1_norm}, height = 500, width = ifelse(length(unique(as.data.frame(dados_norm)[,2]))>2,
  #                                                        150*length(unique(as.data.frame(dados_norm)[,2])),
  #                                                        400) )->output$norm_plot
  #     
  #     
  #   }
  # })

  
  observeEvent(input$house_Norm,{
    
    if(!is.null(datos$x)){
      dados_norm<-datos$x%>%group_by(Name)%>%mutate(row=row_number())%>%tidyr::pivot_wider(names_from = Name, values_from = Value)
      dados_norm<-dados_norm[,-3]
      

      if(length(unique(table(dados_norm[,1])))>1){
        dados_norm<-datos$x%>%group_by(Name)%>%
          tidyr::pivot_wider(names_from = Name, values_from = Value)
      }
      dados_norm<-unnest(dados_norm,colnames(dados_norm[,3:ncol(dados_norm)]))
      
      dados_norm[which((rowSums(is.na(dados_norm))/ncol(dados_norm))<=0.1),]-> dados_norm
      dados_norm[,which((colSums(is.na(dados_norm))/nrow(dados_norm))<=0.1)]-> dados_norm
      na.omit(dados_norm)->dados_norm
      
      if(input$house_analysis=="Yes"){
        input$house_genes->house_genes
        strsplit(house_genes," ")->house_genes
        if(length(house_genes[[1]])>2){
          dados_norm[which(colnames(dados_norm) %in% house_genes[[1]])]->dados_norm2
          if(ncol(dados_norm2)>2){
            cbind(dados_norm[,1:2],dados_norm2)->dados_norm
          }
        }
      }
      
      
      as.data.frame(rbind(t(dados_norm[,-c(1:2)]),t(dados_norm[,2])))->dados_normfinder
      
      
      loadError = FALSE
      toto <- function(){
        
        tryCatch(resul_house<-Normfinder(dados_normfinder),error=function(e)loadError <<- TRUE)
      }
      
      toto()
      if(loadError){
        renderText({"There is no enough complete data for running NormFinder at this moment.
        If you know your housekeeping genes, you can include them in the space below."})->output$house_verif1
        output$out_normal2<-
          renderUI({
            renderText({"Now you need to make important decisions that will affect all further analysis."})
          })
        
        output$out_normal3<-
          renderUI({
            selectInput("which_norm","How do you want to normalize your data?",c("I want to define the housekeeping genes"))
          })
        
        # output$out_normal5<-
        #   renderUI({
        #     checkboxInput("which_outliers","2. Do you want to remove outliers per gene and group?")
        #   })
        
        output$out_normal6<-
          renderUI({
            actionButton("next_Analyses","Let's go to expression and correlation analyses",icon=icon("edge"))
          })
        
        renderText({paste0(" ")})->output$norm_text
        renderText({paste0(" ")})->output$norm_text1
        renderText({paste0(" ")})->output$norm_text2
        
      } else{
        resul_house<-Normfinder(dados_normfinder)
        
        renderText({"This code uses the NormFinder algorithm to determine the best housekeeping gene for your experiment. 
        Please note that the Normfinder algorithm uses only complete observations. So, this is only a suggestion for normalizing your genes.
        Then you can specify your preferred genes for normalization or leave the system to choose."})->output$house_verif1
        data.frame(Gene=rownames(head(resul_house$Ordered,10)),head(resul_house$Ordered,10))->resul_house_ord
        head(resul_house$PairOfGenes[order(resul_house$PairOfGenes[,ncol(resul_house$PairOfGenes)]),],10)->resul_house_pair
        renderTable({resul_house_ord})->output$house_table
        renderTable({resul_house_pair})->output$house_table2
        as.character(resul_house_pair[1,1:2])->datos$ref
        
        output$out_normal2<-
          renderUI({
            renderText({"Now you need to make important decisions that will affect all further analysis."})
          })
        
        output$out_normal3<-
          renderUI({
            selectInput("which_norm","How do you want to normalize your data?",c("Leave the system to choose","I want to define the housekeeping genes"))
          })
        
        # output$out_normal5<-
        #   renderUI({
        #     checkboxInput("which_outliers","2. Do you want to remove outliers per gene and group?")
        #   })
        
        output$out_normal6<-
          renderUI({
            actionButton("next_Analyses","Let's go to expression and correlation analyses",icon=icon("edge"))
          })
        
        renderText({paste0(" ")})->output$norm_text
        renderText({paste0(" ")})->output$norm_text1
        renderText({paste0(" ")})->output$norm_text2
      }
      }
      
      
  })

  output$out_normal4<-
    renderUI({
      req(input$which_norm=="I want to define the housekeeping genes")
      textInput("which_norm_genes","Please input the housekeeping genes for next analyses. Separate genes by space. i.e. ACTB GAPDH")
      })
  
  observeEvent(input$next_Analyses,{
    
    if(!is.null(datos$x)){
      dados_norm<-datos$x%>%
        group_by(Name)%>%
        mutate(row=row_number())%>%
        tidyr::pivot_wider(names_from = Name, values_from = Value)
      
      dados_norm<-dados_norm[,-3]
      
      if(length(unique(table(dados_norm[,1])))>1){
        dados_norm<-datos$x%>%group_by(Name)%>%
          tidyr::pivot_wider(names_from = Name, values_from = Value)
      }

      dados_norm<-unnest(dados_norm,colnames(dados_norm[,3:ncol(dados_norm)]))
      
      if(input$which_norm=="I want to define the housekeeping genes"){

        input$which_norm_genes->which_norm_genes
        strsplit(which_norm_genes," ")->which_norm_genes
        dados_norm[which(colnames(dados_norm) %in% which_norm_genes[[1]])]->dados_norm2
        
        if(ncol(dados_norm2)<1){
          renderText({"Please verify that the input housekeeping gene is present in the data table."})->output$norm_text
          renderText({"Please verify that the input housekeeping gene is present in the data table."})->output$norm_text1
          } else {
            renderText({paste0("Your data will be processed with the next housekeeping gene(s): ",
                               paste0(as.character(colnames(dados_norm2)),collapse=","))})->output$norm_text
            
            renderText({paste0("Your data will be processed with the next housekeeping gene(s): ",
                               paste0(as.character(colnames(dados_norm2)),collapse=","))})->output$norm_text1
            
            colnames(dados_norm2)->normalizers
            
            }
        } else {
          
          renderText({paste0("Your data will be processed with the next housekeeping gene(s): ",
                           paste0(datos$ref,collapse=","))})->output$norm_text
          
          renderText({paste0("Your data will be processed with the next housekeeping gene(s): ",
                           paste0(datos$ref,collapse=","))})->output$norm_text1
          
          datos$ref->normalizers
          }

      ##Normalizando os resultados
      dados_norm[which(colnames(dados_norm) %in% normalizers)]->c
      dados_norm$endo<-apply(c,1,mean)
      
      if(nrow(c[complete.cases(c), ]) !=nrow(c)){
        which(rowSums(is.na(dados_norm[which(colnames(dados_norm) %in% normalizers)]))>0)->c
        
        dados_norm[-c,]->dados_norm
      }
      
      dados_norm[which(!is.na(dados_norm$endo)),]->dados_norm
      
      ##Calculando valores -dCt
      
      partial<-data.frame(apply(dados_norm[,-c(1:2)],2,function(w){dados_norm[,ncol(dados_norm)] - w }))
      cbind(dados_norm[,1:2],partial)->partial
      colnames(partial)<-colnames(dados_norm)
      partial->dados_norm
      dados_norm[,-ncol(dados_norm)]->dados_norm
      unique(dados_norm[,2])->nome_grupo

      # if(input$which_outliers==TRUE){
      #   
      #   for (i in 3:ncol(dados_norm)){
      #     for(k in 1:length(nome_grupo)){
      #       quartiles <- quantile(dados_norm[which(dados_norm[,2]==as.character(nome_grupo[k])),i], probs=c(.25, .75), na.rm = TRUE)
      #       IQR <- IQR(as.matrix(dados_norm[which(dados_norm[,2]==as.character(nome_grupo[k])),i]), na.rm = TRUE)
      #       
      #       Lower <- quartiles[1] - 1.5*IQR
      #       Upper <- quartiles[2] + 1.5*IQR 
      #       
      #       dados_norm[which(dados_norm[,2]==as.character(nome_grupo[k]) & dados_norm[,i]<Lower),i]<-NA
      #       dados_norm[which(dados_norm[,2]==as.character(nome_grupo[k]) & dados_norm[,i]> Upper),i]<-NA
      #     }
      #   }
      #   
      #   renderText({paste0("Your data will be processed with the next housekeeping gene(s): ",
      #                      paste0(datos$ref,collapse=","))})->output$norm_text
      #   renderText({paste0("Your data will be processed with the next housekeeping gene(s): ",
      #                      paste0(datos$ref,collapse=","))})->output$norm_text1
      # }
      
      # dados_norm[which((rowSums(is.na(dados_norm))/ncol(dados_norm))<=0.3),]-> dados_norm
      # dados_norm[,which((colSums(is.na(dados_norm))/nrow(dados_norm))<=0.3)]-> dados_norm


      
      if(input$house_analysis=="Yes" && input$which_norm!="I want to define the housekeeping genes"){
        input$house_genes->house_genes
        strsplit(house_genes," ")->house_genes
        if(length(house_genes[[1]])>2){
          dados_norm[which(!(colnames(dados_norm) %in% house_genes[[1]]))]->dados_norm
        }
      }
      
      if(input$which_norm=="I want to define the housekeeping genes"){
        input$which_norm_genes->which_norm_genes
        strsplit(which_norm_genes," ")->which_norm_genes
        dados_norm[which(!(colnames(dados_norm) %in% which_norm_genes[[1]]))]->dados_norm
      }
      
      aggregate(x=dados_norm,
                by=list(dados_norm$Sample,
                        dados_norm$Group),
                FUN=mean)->dados_norm
      
      dados_norm[,-c(3:4)]->dados_norm
      
      colnames(dados_norm)<-c("Sample","Group",colnames(dados_norm[,-c(1:2)]))
      
      dados_norm->datos$final
      
      tail(unlist(strsplit(input$xlsx_input$datapath,"[.]")),1)->formato
      if(formato=="xlsx"){
      if(length(excel_sheets(input$xlsx_input$datapath))>1){
        renderText({paste0("We are processing secondary variables...")})->output$norm_text2
        other_variables<-read_excel(input$xlsx_input$datapath,sheet = 2)
        
        if(ncol(other_variables)!=0){
          colnames(other_variables)<-c("Samples",colnames(other_variables[,-1]))
          other_variables<-other_variables%>%filter(Samples %in% unique(dados_norm$Sample))
          
          if(nrow(other_variables)==length(unique(dados_norm$Sample))){
            as.data.frame(other_variables)->datos$o_var
            
            
          } else {
            renderText({paste0("Please revise secondary variables for all analyzed samples")})->output$norm_text2
          }
          
          }
        
      }
      }
      
      updateSelectInput(session,"ref.group","Select the reference group (for statistical comparisons)",
                                          paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))
       
      updateSelectInput(session,"gene","Select one gene from the list",
                        sort(colnames(datos$final[,-c(1:2)])))

      
      updateSelectInput(session,"gloss_specie","Select the specie you are working with",
                        c("Human"="hsapiens_gene_ensembl", "D. melanogaster"="dmelanogaster_gene_ensembl",
                          "Mouse"="mmusculus_gene_ensembl", "Rat"="rnorvegicus_gene_ensembl"))
      
      updateSelectInput(session,"gloss_type","Choose the format of your gene list",
                        c("Gene Name"="external_gene_name",
                          "Ensembl Gene ID"="ensembl_gene_id",
                          "Ensembl Transcript ID"="ensembl_transcript_id",
                          "NCBI Entrez ID"="entrezgene_id"))
      
      if(is.null(datos$o_var)){
        updateSelectInput(session,"group0","Select the variable to be compared",
                          c("Primary comparison"))
        
        updateSelectInput(session,"ref.group0","Select the variable to be compared",
                          c("Primary comparison"))
        
        updateSelectInput(session,"Treatment0","Select the variable to be compared",
                          c("Primary comparison"))
        
        updateSelectInput(session,"sc_group0","Select the variable to be compared",
                          c("Primary comparison"))
        

      } else {
        updateSelectInput(session,"group0","Select the variable to be compared",
                          c("Primary comparison",colnames(datos$o_var[,-1])))
        
        updateSelectInput(session,"ref.group0","Select the variable to be compared",
                          c("Primary comparison",colnames(datos$o_var[,-1])))
        
        updateSelectInput(session,"Treatment0","Select the variable to be compared",
                          c("Primary comparison",colnames(datos$o_var[,-1])))
        
        updateSelectInput(session,"sc_group0","Select the variable to be compared",
                          c("Primary comparison",colnames(datos$o_var[,-1])))
      }
      
      
      updateSelectInput(session,"group1","Select the first group for comparison",
                        paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))
      
      updateSelectInput(session,"group2","Select the second group for comparison",
                        paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))
      
      updateSelectInput(session,"Treatment2","Choose an specific group",
                        paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))
      
      updateSelectInput(session,"sc_group1","Choose an specific group",
                        paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))
      
      updateSelectInput(session,"sc_group2","Choose an specific group",
                        paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))
      
      updateSelectInput(session,"sc_gene1","Select the first gene for correlation analysis",
                        sort(colnames(datos$final[,-c(1:2)])))
      
      updateSelectInput(session,"sc_gene2","Select the second gene for correlation analysis",
                        sort(colnames(datos$final[,-c(1:2)])))
      
      
      renderText({paste0("Your data was correctly processed!")})->output$norm_text2
      
      output$out_resul1 <- renderUI({
        
        downloadButton("downloadConvertedData","Download -\u25B3CT values per sample")         
      })
      
      
      output$downloadConvertedData <- downloadHandler(
        filename = function() {
          gsub(":","_",Sys.time())->timed
          paste("normalized_Data-", timed, ".csv", sep="")
        },
        content = function(file) {
          write.csv(datos$final, file,row.names = F)
        }
      )
      
      
    }
  })
  

  observe({
    
    input$ref.group0->prim_com
    if(prim_com!="Primary comparison"){
      
      
      updateSelectInput(session,"ref.group",
                        "Select the reference group (for statistical comparisons)",
                        sort(paste0(t(data.frame(unique(datos$o_var[,colnames(datos$o_var) %in% prim_com]))))))
      
      updateCheckboxGroupInput(session,"ref.group0_5",
                               "Select groups to be compared",
                               sort(paste0(t(data.frame(unique(datos$o_var[,colnames(datos$o_var) %in% prim_com]))))))
      
    } else {
      updateSelectInput(session,"ref.group","Select the reference group (for statistical comparisons)",
                        paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))
      
      updateCheckboxGroupInput(session,"ref.group0_5",
                               "Select groups to be compared",
                               paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))
      
      
    }

    
    
    })
  
  # output$UIinputexpr1<-
  #   renderUI({
  #     req(input$specific_expr)
  #     selectInput("which_focus_genes","Please choose the gene",sort(colnames(datos$final[,-c(1:2)])))
  #   })
  # 
  # output$UIinputexpr2<-
  #   renderUI({
  #     req(input$specific_expr)
  #     selectInput("which_focus_level","Now choose the expression level",c("High expression (Q4)","Low expression (Q1)"))
  #   })

  observeEvent(input$run_button_Exp,{
    
    
    req(datos$final)
    
    if(!is.null(input$gene) & !is.null(input$ref.group)){
      
      input$gene->gene
      input$ref.group0_5 -> variables
      input$ref.group->ref.group
      dados_proc2<-datos$final
      input$ref.group0->prim_com
      
      
      if(length(variables)>0){
        if(prim_com!="Primary comparison"){
          
          datos$o_var[,colnames(datos$o_var) %in% c("Samples",prim_com)]->dados_proc3
          merge(dados_proc2,dados_proc3,by.x="Sample","Samples")->dados_proc2
          dados_proc2[,ncol(dados_proc2)]->dados_proc2[,2]
          dados_proc2[,-ncol(dados_proc2)]->dados_proc2
          colnames(dados_proc2)<-c("Sample","Group",colnames(dados_proc2[,-c(1:2)]))
        }
        
        if(ref.group %in% variables){
          
          dados_proc2[which(dados_proc2$Group %in% variables),]->dados_proc2
        }
        
        which(colnames(dados_proc2)==gene)->coluna
        
        filt<-dados_proc2[which(!is.na(dados_proc2[,coluna])),c(2,coluna)]
        
        names(table(as.matrix(filt[,1])))->rot
        NULL->resfin
        "parametric"->norm_test
        
        for(i in 1:length(rot)){
          ks.test(filt[which(filt$Group==rot[i]),2],"pnorm")[2]->res
          if(res < 0.01){
            "no_parametric"->norm_test
          }
          c(resfin,res)->resfin
        }
        
        data.frame(Group=rot,KS_p.value=as.numeric(resfin))->KS_result
        
        formatC(as.numeric(KS_result[,2]),format="e")->KS_result[,2]
        
        KS_result[which(KS_result[,2]< 1e-16),2]<-"< 1e-16"
        
        renderTable({KS_result})->output$NormTable
        
        if(is.null(datos$order_color)){
          niveis<-paste0(t(unique(dados_proc2[order(dados_proc2[,2]),][,2])))
          hue_pal()(length(unique(dados_proc2[,2])))->col.palette
          as.matrix(unique(dados_proc2[order(dados_proc2[,2]),][,2]))->lista_grupos
          # }else{
          #   niveis<-paste0(as.matrix(datos$order_color[order(datos$order_color[,2]),][,1]))
          #   as.character(as.matrix(datos$order_color[order(datos$order_color[,2]),][,3]))->col.palette
          #   as.matrix(datos$order_color[order(datos$order_color[,2]),][,1])->lista_grupos
        }
        
        if(input$method_exp=="Non-parametric test (Mann-Whitney)"){
          method_exp<-"wilcox.test"
          method_exp_t<-"Mann-Whitney test"
        } else if (input$method_exp=="Parametric test (T-test)"){
          method_exp<-"t.test"
          method_exp_t<-"T-test"
        } else {
          
          if(norm_test=="parametric"){
            method_exp<-"t.test"
            method_exp_t<-"T-test"
          }else{
            method_exp<-"wilcox.test"
            method_exp_t<-"Mann-Whitney test"
          }
          
        }
        
        
        dados_proc2$Group <- factor(dados_proc2$Group , levels=niveis)
        
        p2<-dados_proc2[which(!is.na(dados_proc2[,coluna])),c(2,coluna)]%>%
          ggplot(aes(Group,get(gene),color=Group))+geom_boxplot()+geom_jitter(width=0.1,size=1,alpha=0.7)+
          theme_classic(base_size=18)+
          scale_color_manual( values=col.palette,breaks=paste0(lista_grupos))+ 
          theme(axis.text.x = element_text(angle = 45,hjust=1))+
          labs(title=paste(gene,"levels across groups"),x="Groups",y=bquote(italic(.(gene)) ~ levels (-~Delta~Ct)))+
          stat_compare_means(aes(label = ..p.signif..),ref.group=ref.group,method=method_exp)
        
        renderPlot({p2}, height = 500, width = ifelse(length(unique(as.data.frame(dados_proc2)[,2]))>2,
                                                      150*length(unique(as.data.frame(dados_proc2)[,2])),
                                                      400) )->output$express2Plot
        
        
        plotdCt<-reactive({p2})
        
        output$UIexpr1 <- renderUI({
          
          downloadButton("downloadPlotdCt","Download this plot")         
        })
        
        output$downloadPlotdCt <- downloadHandler(
          filename = paste0(gene," levels (-dCt) across groups",".png"),
          content = function(file) {
            png(file,height = 5000, width = ifelse(length(unique(as.data.frame(dados_proc2)[,2]))>2,
                                                   1500*length(unique(as.data.frame(dados_proc2)[,2])),
                                                   4000),res=600)
            print(plotdCt())
            dev.off()
          }) 
        
        
        
        
        dados_proc2->dados_calc
        dados_proc2[,c(3:ncol(dados_proc2))]-median(as.matrix(dados_proc2[which(!is.na(dados_proc2[,coluna]) & dados_proc2$Group==ref.group),coluna]))->dados_calc[,3:ncol(dados_proc2)]
        2^dados_calc[,3:ncol(dados_proc2)]->dados_calc[,3:ncol(dados_proc2)]
        
        dados_calc$Group <- factor(dados_calc$Group , levels=niveis)
        
        
        p1<-dados_calc[which(!is.na(dados_calc[,coluna])),]%>% 
          ggplot(aes(Group,get(gene),color=Group))+geom_boxplot()+geom_jitter(width=0.1,size=1,alpha=0.7)+
          theme_classic(base_size=18)+
          scale_color_manual( values=col.palette,breaks=paste0(lista_grupos))+ 
          theme(axis.text.x = element_text(angle = 45,hjust=1))+
          labs(title=paste(gene,"levels across groups"),x="Groups",y=bquote(italic(.(gene)) ~ levels (2^(-~Delta~Delta~Ct))))+
          stat_compare_means(aes(label = ..p.signif..),ref.group=ref.group,method=method_exp)
        renderPlot({p1}, height = 500, width = ifelse(length(unique(as.data.frame(dados_proc2)[,2]))>2,
                                                      150*length(unique(as.data.frame(dados_proc2)[,2])),
                                                      400) )->output$expressPlot
        
        
        plotddCt<-reactive({p1})
        
        output$UIexpr2 <- renderUI({
          
          downloadButton("downloadPlotddCt","Download this plot")         
        })
        
        output$downloadPlotddCt <- downloadHandler(
          filename = paste0(gene," levels (2^-ddCt) across groups",".png"),
          content = function(file) {
            png(file,height = 5000, width = ifelse(length(unique(as.data.frame(dados_proc2)[,2]))>2,
                                                   1500*length(unique(as.data.frame(dados_proc2)[,2])),
                                                   4000),res=600)
            print(plotddCt())
            dev.off()
          }) 
        
        
        unique(dados_calc$Group)->grupos
        NULL->matriz
        
        for(i in 1:length(grupos)){
          if(ref.group!=grupos[i]){
            
            dados_calc[which(dados_calc$Group==ref.group | dados_calc$Group==grupos[i]),]->brief
            if(method_exp=="wilcox.test"){
              wilcox.test(get(gene)~Group,data=brief)->resul
            } else {
              t.test(get(gene)~Group,data=brief)->resul
            }
            
            median( as.matrix(na.omit(dados_calc[which(dados_calc$Group==ref.group),which(colnames(dados_calc)==gene)])))->a
            median( as.matrix(na.omit(dados_calc[which(dados_calc$Group==grupos[i]),which(colnames(dados_calc)==gene)])))->b
            round(IQR( as.matrix(na.omit(dados_calc[which(dados_calc$Group==ref.group),which(colnames(dados_calc)==gene)]))),2)->c
            round(IQR( as.matrix(na.omit(dados_calc[which(dados_calc$Group==grupos[i]),which(colnames(dados_calc)==gene)]))),2)->d
            
            c(ref.group,as.character(grupos[i]),
              round(length(as.matrix(na.omit(brief[which(brief$Group==ref.group),which(colnames(dados_calc)==gene)]))),digits=0),
              round(length(as.matrix(na.omit(brief[which(brief$Group==grupos[i]),which(colnames(dados_calc)==gene)]))),digits=0),
              paste0(round(a,2)," (",c,")"),
              paste0(round(b,2)," (",d,")"),
              a/b,
              resul[3],
              method_exp_t)->linha
            
            rbind(matriz,linha)->matriz
            
          }
          
        }
        
        colnames(matriz)<-c("Group1","Group2","n (G1)", "n (G2)","median (IQR) for G1 (2^-ddCt)","median (IQR) for G2 (2^-ddCt)","Fold Change","p.value", "test")  
        formatC(as.numeric(matriz[,8]),format="e")->matriz[,8]
        renderTable(matriz)->output$expressTable
        
        output$idExp1 <- renderUI({
          
          downloadButton("downloadData", "Download")
          
        })
        
        output$downloadData <- downloadHandler(
          filename = function() {
            gsub(":","_",Sys.time())->timed
            paste(gene," expressionData-", timed, ".csv", sep="")
          },
          content = function(file) {
            write.csv(matriz, file,row.names = F)
          }
        )
      }
      }
      
      
      
     
    
    
    
  })
  
  output$UIinputvolc1<-
    renderUI({
      req(input$specific_volc)
      selectInput("which_focus_genes_volc","Please choose the gene",sort(colnames(datos$final[,-c(1:2)])))
    })
  
  output$UIinputvolc2<-
    renderUI({
      req(input$specific_volc)
      selectInput("which_focus_level_volc","Now choose the expression level",c("High expression (Q4)","Low expression (Q1)"))
    })
  
  observe({

    input$group0->prim_com_volc
    if(prim_com_volc!="Primary comparison"){


      updateSelectInput(session,"group1",
                        "Select the first group for comparison",
                        sort(paste0(t(data.frame(unique(datos$o_var[,colnames(datos$o_var) %in% prim_com_volc]))))))

      updateSelectInput(session,"group2",
                        "Select the second group for comparison",
                        sort(paste0(t(data.frame(unique(datos$o_var[,colnames(datos$o_var) %in% prim_com_volc]))))))

    } else {
      updateSelectInput(session,"group1","Select the first group for comparison",
                        paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))

      updateSelectInput(session,"group2","Select the second group for comparison",
                        paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))

    }



  })
  
  observeEvent(input$run_button_volc,{
    
    req(datos$final)
   

    if(!is.null(input$group1) & !is.null(input$group2)){
      
      input$group1->group1
      input$group2->group2
      ifelse(is.na(as.numeric(input$p_val_volc)),1e-4,ifelse(between(as.numeric(input$p_val_volc),0,1),as.numeric(input$p_val_volc),1e-4))->p_volc
      ifelse(is.na(as.numeric(input$FC_volc)),2,as.numeric(input$FC_volc))->FC_volc

      datos$final->dados3
      input$group0->prim_com_volc

      # if(input$specific_volc==TRUE){
      #   
      #   if(length(which(is.na(dados3[colnames(dados3) %in% input$which_focus_genes_volc])))>0){
      #     dados3<-dados3[-which(is.na(dados3[colnames(dados3) %in% input$which_focus_genes_volc])),]
      #     
      #   }
      #   
      #   dados3<-dados3%>%group_by(Group)%>%mutate(expr_group=ifelse(get(input$which_focus_genes_volc)>quantile(get(input$which_focus_genes_volc),0.75),"Q4",
      #                                                                         ifelse(get(input$which_focus_genes_volc)>quantile(get(input$which_focus_genes_volc),0.5),"Q3",
      #                                                                                ifelse(get(input$which_focus_genes_volc)>quantile(get(input$which_focus_genes_volc),0.25),"Q2","Q1"))))
      #   if(input$which_focus_level_volc=="High expression (Q4)"){
      #     dados3<-dados3[which(dados3$expr_group=="Q4"),-ncol(dados3)]
      #   }else{
      #     dados3<-dados3[which(dados3$expr_group=="Q1"),-ncol(dados3)]
      #   }
      # }
      
      if(prim_com_volc!="Primary comparison"){
        
        datos$o_var[,colnames(datos$o_var) %in% c("Samples",prim_com_volc)]->dados4
        merge(dados3,dados4,by.x="Sample","Samples")->dados3
        dados3[,ncol(dados3)]->dados3[,2]
        dados3[,-ncol(dados3)]->dados3
        colnames(dados3)<-c("Sample","Group",colnames(dados3[,-c(1:2)]))
      }
      
      
      
      if(group1!=group2){
        dados3[which(dados3$Group==group1 | dados3$Group==group2),]->dados3
        dados3[,-1]->dados3
        
        colnames(dados3)->nomes
        
        NULL->matriz
        
        as.data.frame(dados3)->dados3
        apply(dados3[,-1],2,as.numeric)->dados3[,-1]
        
        if(input$method_volc=="Non-parametric test (Mann-Whitney)"){
          for(i in 2:ncol(dados3)){
            
            if(length(as.matrix(unique(dados3[which(!is.na(dados3[,i])),1])))==2){
              as.numeric(wilcox.test(get(nomes[i])~Group,dados3)[3])->result_vol
              
              2^median(as.matrix(na.omit(dados3[which(dados3$Group==group1),i])))->a
              2^median(as.matrix(na.omit(dados3[which(dados3$Group==group2),i])))->b
              
              c(paste(group1,"vs.",group2),nomes[i],length(as.matrix(na.omit(dados3[which(dados3$Group==group1),i]))),
                length(as.matrix(na.omit(dados3[which(dados3$Group==group2),i]))),a,b,a/b,result_vol)->linha
              
              rbind(matriz,linha)->matriz
            }
          }
          
          
        } else {
          
          for(i in 2:ncol(dados3)){
            
            if(length(as.matrix(unique(dados3[which(!is.na(dados3[,i])),1])))==2){
              
              if(table(!(is.na(dados3[which(dados3$Group==group1),][i])))["TRUE"]>1 &
                 table(!(is.na(dados3[which(dados3$Group==group2),][i])))["TRUE"]>1){
                
              as.numeric(t.test(get(nomes[i])~Group,dados3)[3])->result_vol
              
              2^median(as.matrix(na.omit(dados3[which(dados3$Group==group1),i])))->a
              2^median(as.matrix(na.omit(dados3[which(dados3$Group==group2),i])))->b
              
              c(paste(group1,"vs.",group2),nomes[i],length(as.matrix(na.omit(dados3[which(dados3$Group==group1),i]))),
                length(as.matrix(na.omit(dados3[which(dados3$Group==group2),i]))),a,b,a/b,result_vol)->linha
              
              rbind(matriz,linha)->matriz
              
              }
            }
          }
          
        }
        
        
        data.frame(matriz,adjusted_p_BH=p.adjust(matriz[,8],method="BH"))->matriz
        as.numeric(matriz[,7])->matriz[,7]
        
        data.frame(matriz[,1:7],log2FC=log2(matriz[,7]),matriz[,8:9])->matriz
        as.numeric(matriz[,9])->matriz[,9]
        as.numeric(matriz[,10])->matriz[,10]
        colnames(matriz)<-c("Test","gene","n_G1","n_G2","G1","G2","FC","log2FC","pval","adjP")
        
        
        if(input$p_plotted=="adjusted p-value (BH)"){
          "adjP"->test_p
          "p-values were adjusted using the BH method"->text_p
        }else{
          "pval"->test_p
          "p-values were not adjusted"->text_p
        }
        
        pvolc<-EnhancedVolcano(matriz,lab=matriz$gene,x="log2FC",y=test_p,title=paste(group1,"vs.",group2),
                               drawConnectors = TRUE,subtitle=paste0(text_p," || p-value<",as.character(p_volc)," || FC>",FC_volc),
                               pCutoff=p_volc,FCcutoff=log2(FC_volc))
        
        renderPlot({
          pvolc
        },height=600,width=600)->output$volcanoPlot
        
        
        plot_volc<-reactive({pvolc})
        
        output$UIvolcano1 <- renderUI({
          
          downloadButton("downloadPlotvolc","Download this plot")         
        })
        
        output$downloadPlotvolc <- downloadHandler(
          filename = paste0(paste(group1,"vs.",group2)," volcano_plot",".png"),
          content = function(file) {
            png(file,height = 6000, width = 6000,res=600)
            print(plot_volc())
            dev.off()
          }) 
        
        
        matriz[order(matriz$adjP),]->matriz2
        formatC(as.numeric(matriz2[,9]),format="e")->matriz2[,9]
        formatC(as.numeric(matriz2[,10]),format="e")->matriz2[,10]
        
        apply(matriz2[,1:2],2,as.character)->matriz2[,1:2]
        apply(matriz2[,3:8],2,as.numeric)->matriz2[,3:8]
        colnames(matriz2)<-c("Comparison","gene","n (G1)","n (G2)","expression levels in G1 (2^-dCt)","expression levels in G2 (2^-dCt)","Fold Change","log2FC","p-value","adjusted p-value")
        
        renderTable({matriz2})->output$volcanoTable
        
        output$UIvolcano2 <- renderUI({
          
          downloadButton("download_Data", "Download")
          
        })
        
        output$download_Data <- downloadHandler(filename = function() {
          gsub(":","_",Sys.time())->timed
          paste(paste(group1,"vs.",group2),timed, ".csv", sep="")
        },
        content = function(file1) {
          write.csv(matriz2, file1,row.names = F)
        }
        )
      }
    }
    
  })
  
  observe({
    
    input$Treatment0->prim_com_treat
    if(prim_com_treat!="Primary comparison"){
      
      updateSelectInput(session,"Treatment2","Choose an specific group",
                        sort(paste0(t(data.frame(unique(datos$o_var[,colnames(datos$o_var) %in% prim_com_treat]))))))
      
    } else {
      updateSelectInput(session,"Treatment2","Choose an specific group",
                        paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))
      
    }
    
    
    
  })
  
  
  observeEvent(input$run_button_correlat,{
    
    req(datos$final)
    
    
    if(!is.null(input$Treatment2)){

    datos$final->mat_correlat
    input$Treatment0->prim_com_treat
    mat_correlat[,-1]->mat_correlat
    
    
    if(prim_com_treat!="Primary comparison"){
      datos$final->mat_correlat
      datos$o_var[,colnames(datos$o_var) %in% c("Samples",prim_com_treat)]->mat_correlat2
      merge(mat_correlat,mat_correlat2,by.x="Sample","Samples")->mat_correlat
      mat_correlat[,ncol(mat_correlat)]->mat_correlat[,2]
      mat_correlat[,-ncol(mat_correlat)]->mat_correlat
      colnames(mat_correlat)<-c("Sample","Group",colnames(mat_correlat[,-c(1:2)]))
      mat_correlat[,-1]->mat_correlat
    }
    
    
    names(table(mat_correlat$Group))->grupos_corr
    
    NULL->b
    
    withProgress(message = 'Estimating correlation values...', value = 0, {
    for(i in 1:length(grupos_corr)){
      rcorr(as.matrix(mat_correlat[which(mat_correlat$Group==grupos_corr[i]),-1]))->a
      
      if(is.null(b)){
        flattenCorrMatrix(a$r, a$P)->b
      }else{
        flattenCorrMatrix(a$r, a$P)->c
        cbind(b,c[,3:4])->b
      }
      
      incProgress(1/length(grupos_corr), detail = paste("Data from group", i))
    }
      
    })
    
    b->corr_matrix
    
    corr_matrix->corr_matrix2

    limit<-input$valor_p/100
    
    
    
    for(i in c((2*(1:length(grupos_corr)))+2)){
      corr_matrix2[which(corr_matrix2[,i]>limit),(i-1)]<-0
    }
    
    grupos_corr->grupos_corr2
    
    
    which(grupos_corr2==input$Treatment2)->base_corr
    1:length(grupos_corr)->compar_corr
    compar_corr[-base_corr]->compar_corr
    
    matrix()->result
    NULL->labels
    
    if(length(compar_corr)==1){
      abs(corr_matrix2[(base_corr*2)+1]-corr_matrix2[(compar_corr*2)+1])->x
      cbind(result,x)->result
      
      paste0(grupos_corr2[base_corr],"\nvs.\n",grupos_corr2[compar_corr])->xx
      c(labels,xx)->labels
      result[,-1]->result
      
    }else{
      
      for(i in 1:length(compar_corr)){
        abs(corr_matrix2[(base_corr*2)+1]-corr_matrix2[(compar_corr[i]*2)+1])->x
        cbind(result,x)->result
        
        paste0(grupos_corr2[base_corr],"\nvs.\n",grupos_corr2[compar_corr[i]])->xx
        c(labels,xx)->labels
      }
      result[,-1]->result
      
      colnames(result)<-labels
    }

    

    data.frame(Par=paste0(corr_matrix2$row,"_",corr_matrix2$column),result)->result
    result2<-pivot_longer(result,!Par,names_to="comparison",values_to="values")
    
    ggplot(result2,aes(comparison,values,text=paste0("Pair: ",Par,"\n","dR-val:", round(values,digits=2))))+
      geom_jitter(width=0.25,size=1.5,color="deepskyblue4",alpha=0.7)+theme_classic()+
      ylim(c(-.1,max(result2$values)+.2))+
      geom_hline(yintercept=0.7,linetype="dashed",color="red")+
      labs(y="|Rval1-Rval2|")+
      scale_x_discrete(labels=labels)->p1
    

      ggplotly(p1, tooltip = "text", height = 700, width = 1100)->p1
      p1<-p1%>% layout(title = list(text = paste0('\nScatter plot of differentially correlated genes',
                                                  '<br>',
                                                  '<sup>',
                                                  "All correlation values with p-values higher than ",limit," were replaced with 0",'</sup>')))
    
    
    renderPlotly({p1})->output$correlatPlot
    
    }
  })
  
  observe({
    
    input$sc_group0->prim_com_sc
    if(prim_com_sc!="Primary comparison"){

      updateSelectInput(session,"sc_group1","Choose the first group to be plotted",
                        sort(paste0(t(data.frame(unique(datos$o_var[,colnames(datos$o_var) %in% prim_com_sc]))))))
      
      updateSelectInput(session,"sc_group2","Choose the second group to be plotted",
                        sort(paste0(t(data.frame(unique(datos$o_var[,colnames(datos$o_var) %in% prim_com_sc]))))))
    } else {
      
      updateSelectInput(session,"sc_group1","Choose the first group to be plotted",
                        paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))
      
      updateSelectInput(session,"sc_group2","Choose the second group to be plotted",
                        paste0(t(unique(datos$final[order(datos$final[,2]),][,2]))))

    }
    
    
    
  })
  
  observeEvent(input$run_button_sc,{
    
    req(datos$final)

    if(!is.null(input$sc_group1) & !is.null(input$sc_group2)){
      
      datos$final->dados_proc
      
      input$sc_group0->prim_com_sc
      
      if(prim_com_sc!="Primary comparison"){
        datos$final->dados_proc
        datos$o_var[,colnames(datos$o_var) %in% c("Samples",prim_com_sc)]->dados_proc21
        merge(dados_proc,dados_proc21,by.x="Sample","Samples")->dados_proc
        dados_proc[,ncol(dados_proc)]->dados_proc[,2]
        dados_proc[,-ncol(dados_proc)]->dados_proc
        colnames(dados_proc)<-c("Sample","Group",colnames(dados_proc[,-c(1:2)]))
      }
      
      
      
    dados_proc[which(dados_proc$Group==input$sc_group1 |dados_proc$Group==input$sc_group2 ),]->dados_sc
    
    dados_sc[which(colnames(dados_sc)==input$sc_gene1 | colnames(dados_sc)==input$sc_gene2 | colnames(dados_sc)=="Group")]->dados_sc
    
    
    if(is.null(datos$order_color)){
      hue_pal()(length(unique(datos$x[,2])))->col.palette
      as.matrix(unique(datos$x[order(datos$x[,2]),][,2]))->lista_grupos
    }else{
      as.character(as.matrix(datos$order_color[order(datos$order_color[,2]),][,3]))->col.palette
      as.matrix(datos$order_color[order(datos$order_color[,2]),][,1])->names(col.palette)
    }
    
    as.character(dados_sc$Group)->dados_sc$Group
    
    sc_p<-ggplot(dados_sc,aes(get(input$sc_gene1),get(input$sc_gene2),color=Group))+
      geom_point()+
      facet_wrap(Group~.,scales = "free")+
      theme_bw(base_size=18)+stat_cor(output.type = "text",size=6)+
      labs(x=bquote(italic(.(input$sc_gene1))),y=bquote(italic(.(input$sc_gene2))))+
      scale_color_manual(values=col.palette,breaks=c(input$sc_group1,input$sc_group2))+
      geom_smooth(method=lm, se=T, size=.6)
    
    renderPlot({sc_p},width=1000)->output$scatterPlot
    
    plot_scatter<-reactive({sc_p})
    
    output$UIscatter1 <- renderUI({
      
      downloadButton("downloadScatter","Download this plot")         
    })
    
    output$downloadScatter <- downloadHandler(
      filename = paste0(input$sc_gene1," and ", input$sc_gene2,"scatter plot.png"),
      content = function(file) {
        png(file,height = 4000, width = 10000,res=600)
        print(plot_scatter())
        dev.off()
      }) 
    
    }
    
  })
  
  
  
  observeEvent(input$run_gloss,{
    
    req(datos$final)
    
    # output$gloss_id0 <-
    #   renderUI({
    #     shinycssloaders::withSpinner(
    #       dataTableOutput("gloss")
    #     )
    #   })
    
    spec<-NA
    
    
    if(input$gloss_specie=="Human"){
      spec<-"hsapiens_gene_ensembl"
    }else if(input$gloss_specie=="dmelanogaster_gene_ensembl"){
      spec<-"Drosophila_melanogaster"
    }else if(input$gloss_specie=="mmusculus_gene_ensembl"){
      spec<-"Mus_musculus"
    }else if(input$gloss_specie=="rnorvegicus_gene_ensembl"){
      spec<-"Rattus_norvegicus"
    }
    

    if(!is.null(input$gloss_specie) & !is.null(input$gloss_type)){
    #Buscando os nomes de genes
      
    withProgress(message = 'Looking for your genes...', value = 0, {
      
     

    mart <- useDataset(input$gloss_specie, useMart("ensembl"))
    
    incProgress(1/4)

    xy<-getBM(filters= input$gloss_type, attributes= c("ensembl_gene_id",
                                                        "external_gene_name","description",input$gloss_type,"entrezgene_id",
                                                       "entrezgene_description","ensembl_gene_id_version"),
              values=as.matrix(colnames(datos$final[,-c(1:2)])),mart= mart)
    
    version_spec<-as.character(searchDatasets(mart = mart, pattern = input$gloss_specie)[3])
    
    incProgress(1/4)
    
   

    
    if(nrow(xy)>0){
      data.frame(Query=xy[,4],Gene=xy[,2],Description=xy[,3],
               Ensembl_ID=xy[,1],Entrez_ID=paste0(xy[,5],": ",xy[,6]), entrez=xy[,5],version=xy[,7],
               `Genome Version`=version_spec)->glosa

    glosa<-glosa %>% mutate(GeneCards = paste0("<a href=' https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                                               Gene,"' target='_blank'>",Gene,"</a>"))
    } else{
      data.frame(Query=NA,Gene=NA,Description=NA,
                 Ensembl_ID=NA,Entrez_ID=NA, entrez=NA,version=NA,
                 `Genome Version`=NA,GeneCards=NA)->glosa
    }
    
    incProgress(1/4)

   
    
    glosa<-glosa %>% mutate(Ensembl_ID = paste0("<a href='  https://www.ensembl.org/",spec,"/Gene/Summary?db=core;g=",
                                                Ensembl_ID,"' target='_blank'>",version,"</a>"))
    
    glosa<-glosa %>% mutate(NCBI = paste0("<a href='  https://www.ncbi.nlm.nih.gov/gene/",
                                                entrez,"' target='_blank'>",Gene,"</a>"))
    
   
    glosa[,c(1,8,5,2,3,4,9,10)]->glosa
    incProgress(1/4)
    
    
    })

    output$gloss<- renderDataTable({glosa[order(glosa$Gene),]},escape = FALSE)
      
      
    }
    
  })
  
}




shinyApp(ui = ui, server = server)
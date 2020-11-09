# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://www.rstudio.com/shiny/
#
rm(list=ls())
library(shiny)
library(GenomicAlignments)
library(GenomicFeatures)
library(ggplot2)
library(xtable)
require(rlist)


args <- commandArgs(trailing=TRUE)
 
#Args list:
#args[[1]]: first map file
#args[[2]]: second map file
#args[[3]]: chr

data_set <- args[[1]]
mode_set <- args[[2]]
# mode_set <- "ALL"
# source("ChimersBaseFrequences.R")
if (mode_set=="HOM"){
    source("~/scripts/r_scripts/chimeval_blocks_1_functions_hom.r")
} else {
    source("~/scripts/r_scripts/chimeval_blocks_1_functions.r")
}

# source("~/scripts/r_scripts/chimeval_blocks_1_functions.r")

# path_folder <- "/home/cocca/analyses/michelangelo/chimerismo/20AMP/100bp_length"
# path_folder <- "/home/cocca/analyses/michelangelo/chimerismo/20AMP/trimmed_reads/sorted_bam/100bp_length"
# path_folder <- "/home/cocca/analyses/michelangelo/chimerismo/20AMP/15072018_phased_bam_NO_A_85"
path_folder <- "/home/cocca/analyses/michelangelo/chimerismo/06092018"
# path_folder <- "/home/cocca/analyses/michelangelo/chimerismo/TEST_17042018"
outfolder <- paste(path_folder,"/TEST_06092018_chimeval_blocks_",mode_set,sep="")
# genfile <- "/home/cocca/analyses/michelangelo/chimerismo/5HotSpot_IAD132502_181_40amp.txt"
genfile <- "/home/cocca/analyses/michelangelo/chimerismo/06092018/chimeval_HLAandNOHLA_752snps.txt"

bamindexfileDonor=paste(path_folder,"/RIC_0%.bai",sep="")
bamfileDonor=paste(path_folder,"/RIC_0%.bam",sep="")

bamindexfileAcceptor=paste(path_folder,"/RIC_100.bai",sep="")
bamfileAcceptor=paste(path_folder,"/RIC_100.bam",sep="")

bamfileChimer_01=paste(path_folder,"/RIC_0,1%.bam",sep="")
bamindexfileChimer_01=paste(path_folder,"/RIC_0,1%.bai",sep="")

bamfileChimer_1=paste(path_folder,"/RIC_1%.bam",sep="")
bamindexfileChimer_1=paste(path_folder,"/RIC_1%.bai",sep="")

# bamfileChimer="/home/cocca/analyses/michelangelo/chimerismo/20AMP/chimera_1.bam"
# bamindexfileChimer="/home/cocca/analyses/michelangelo/chimerismo/20AMP/chimera_1.bai"

###################################################################################
# Set donor == chimera
# bamfileDonor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/Don.bam"
# bamindexfileDonor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/Don.bai"
# bamfileAcceptor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/RIC_40amp.bam"
# bamindexfileAcceptor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/RIC_40amp.bai"
# bamfileChimer="/home/cocca/analyses/michelangelo/chimerismo/20AMP/Don.bam"
# bamindexfileChimer="/home/cocca/analyses/michelangelo/chimerismo/20AMP/Don.bai"

###################################################################################
#We enned to test different combinations
# Set 1
R_D_R <- c(bamfileAcceptor,bamfileDonor,bamfileAcceptor,bamindexfileAcceptor,bamindexfileDonor,bamindexfileAcceptor)
R_D_C1 <- c(bamfileAcceptor,bamfileDonor,bamfileChimer_1,bamindexfileAcceptor,bamindexfileDonor,bamindexfileChimer_1)
R_D_C01<- c(bamfileAcceptor,bamfileDonor,bamfileChimer_01,bamindexfileAcceptor,bamindexfileDonor,bamindexfileChimer_01)
R_D_D <- c(bamfileAcceptor,bamfileDonor,bamfileDonor,bamindexfileAcceptor,bamindexfileDonor,bamindexfileDonor)
# Set 2
D_R_D <- c(bamfileDonor,bamfileAcceptor,bamfileDonor,bamindexfileDonor,bamindexfileAcceptor,bamindexfileDonor)
D_R_C1 <- c(bamfileDonor,bamfileAcceptor,bamfileChimer_1,bamindexfileDonor,bamindexfileAcceptor,bamindexfileChimer_1)
D_R_C01 <- c(bamfileDonor,bamfileAcceptor,bamfileChimer_01,bamindexfileDonor,bamindexfileAcceptor,bamindexfileChimer_01)
D_R_R <- c(bamfileDonor,bamfileAcceptor,bamfileAcceptor,bamindexfileDonor,bamindexfileAcceptor,bamindexfileAcceptor)

#all_set <- c("R_D_R","R_D_C1","R_D_C01","R_D_D","D_R_D","D_R_C1","D_R_C01","D_R_R")
# all_set <- c("R_D_D","D_R_R","R_D_R","R_D_C1","R_D_C01","D_R_D","D_R_C1","D_R_C01")
# data_set <- c("R_D_D")
all_set <- c(data_set)
    
CutoffHom <- 0.9
CutoffHet <- c(0.3,0.6)

#read aplicon files
# ampl_list <- read.table("/home/cocca/analyses/michelangelo/chimerismo/4HotSpot_IAD132502_181_40amp.list",he=F)
ampl_list <- read.table("/home/cocca/analyses/michelangelo/chimerismo/06092018/chimeval_HLAandNOHLA_752snps.list",he=F)
ampl_list[,1] <- as.character(ampl_list[,1])

for (samples_set in all_set){
    # samples_set <- all_set[1]
    current_set <- get(samples_set)
    print(current_set)

    bamfileAcceptor <- current_set[1]
    bamfileDonor <- current_set[2]
    bamfileChimer <- current_set[3]
    bamindexfileAcceptor <- current_set[4]
    bamindexfileDonor <- current_set[5]
    bamindexfileChimer <- current_set[6]

FromAccetor_blocks<-matrix(ncol=1,nrow=dim(ampl_list)[1],data=NA)
rownames(FromAccetor_blocks)<-(ampl_list[,1])

for (j in 1:dim(ampl_list)[1]) {
    # j <- 22
    ampl <- ampl_list[j,1]
    # bedfile <- paste("/home/cocca/analyses/michelangelo/chimerismo/BED_FILES/SNPS/4HotSpot_IAD132502_181_40amp_mod1_",ampl,"_sorted.bed",sep="")
    bedfile <- paste("/home/cocca/analyses/michelangelo/chimerismo/06092018/BED_FILES/SNPS/chimeval_HLAandNOHLA_752snps_",ampl,".bed",sep="")

    FromAccetor_block <- CountChimericDNA_block(bedfile,genfile,bamfileDonor,bamindexfileDonor,bamfileAcceptor,bamindexfileAcceptor,bamfileChimer,bamindexfileChimer,CutoffHom,CutoffHet)
    FromAccetor_blocks[rownames(FromAccetor_blocks)[j],1] <- FromAccetor_block
    
}

Media<-round(mean(FromAccetor_blocks[!is.na(FromAccetor_blocks)]),5)
DeviazioneStandard<-round(sd(FromAccetor_blocks[!is.na(FromAccetor_blocks)]),5)
Error<-round(qnorm(0.975)*DeviazioneStandard/sqrt(length(FromAccetor_blocks[!is.na(FromAccetor_blocks)])),5)
Design<-read.table(file=bedfile,sep="\t")
rownames(Design)<-Design[,4]
Design<-Design[rownames(Design)%in%rownames(FromAccetor_blocks)[!is.na(FromAccetor_blocks)],7:8]
Result<-paste("Predicted Recipient mean",Media*100,"%; Number of Blocks=",length(FromAccetor_blocks[!is.na(FromAccetor_blocks)]),";CI",(Media+Error)*100,"-",(Media-Error)*100,sep=" ")
ResultSummary<-list(Mean=Media,SD=DeviazioneStandard,Err=Error,Summary=Result,AcceptorTable=FromAccetor_blocks[!is.na(FromAccetor_blocks)],SNPS=rownames(FromAccetor_blocks)[!is.na(FromAccetor_blocks)],Design=Design)

print(Result)
dir.create(outfolder, recursive=T)

list.save(ResultSummary, paste(outfolder,"/",samples_set,'_list.yaml',sep=""))
        
}



# shinyServer(function(input, output) {
#     InputPar<-reactive(Parameters(input))
#     output$par<-renderTable(InputPar())
#     Result<-reactive({CountChimericDNA(input$fileBed,input$fileGen,input$fileD,input$fileDi,input$fileA,input$fileAi,input$fileC,input$fileCi,input$Hom,input$Het)})
#     PlotResult<-reactive({DataToPlot(Result())})
#     output$res<-renderPrint(print(unlist(Result())[4]))
#     output$resplot<-renderPlot(ggplot(PlotResult(), aes(x = treatment, y = mean,ymax=0.04,ymin=0)) +
#     geom_errorbar(aes(x = treatment, y = mean,ymin=(as.numeric(mean)-(as.numeric(me))), ymax=(as.numeric(mean)-(as.numeric(me)))), colour="black") +
#     geom_point(position = position_dodge(), stat="identity", fill="blue") +
#     ggtitle("Bar plot with 99% confidence intervals"))
#     Table<-reactive({SNPsTable(Result())})
#     output$tab<-renderTable(Table())
# })

# # This is the user-interface definition of a Shiny web application.
# # You can find out more about building applications with Shiny here:
# #
# # http://www.rstudio.com/shiny/
# #

# library(shiny)

# shinyUI(fluidPage(
# titlePanel("Custom Pipeline Evaluating Chimerism Using NGS"),
# fluidRow(
# sidebarPanel(
# helpText("Please select the .bam file path or ftp for Donor Sample"),
# textInput("fileD", label = h3("Donor bam file"), value = "file.bam"),
# helpText("Please select the .bam file path or ftp for Accettor Sample"),
# textInput("fileA", label = h3("Recipient bam file"), value = "file.bam"),
# helpText("Please select the .bam file path or ftp for Chimeric Sample"),
# textInput("fileC", label = h3("Chimeric bam file"), "value = file.bam "),
# helpText("Please select the Design .bed file path or ftp"),
# textInput("fileBed", label = h3("Bed file"), value = "file.bed"),
# helpText("Please select the Design Genotype file"),
# textInput("fileGen", label = h3("Gen file"), value = "file.txt"),
# helpText("Please select the Thresholds for Homozygous and Heterozygous call"),
# sliderInput("Hom","Threshold Hom:",min = 0.00,max = 1.00,value = 0.90),
# sliderInput("Het","Threshold Het:",min = 0.00,max = 1.00,value = c(0.30,0.60)),
# submitButton("Submit")
# ),
# sidebarPanel(
# helpText("Please select the .bai file obtained from Donor Sample"),
# textInput("fileDi", label = h3("Donor bai file"), value = "file.bai"),
# helpText("Please select the .bai file obtained from Accettor Sample"),
# textInput("fileAi", label = h3("Recipient bai file"), value = " file.bai "),
# helpText("Please select the .bai file obtained from Chimeric Sample"),
# textInput("fileCi", label = h3("Chimeric bai file"), value = " file.bai ")
# ),
# mainPanel(tabsetPanel(
# tabPanel("Info Summary",helpText("…"),
# tableOutput("par")),
# tabPanel("Result Status",textOutput("…")),
# tabPanel("Result Plot",plotOutput("…")),
# tabPanel("SNPs Value",tableOutput("…"))
# )
# )
# )
# ))


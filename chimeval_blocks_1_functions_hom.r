#Starting from donor, acceptor and chimera calculate the donor chimare frequncies.
CountChimericDNA_block<-function(bedfile,genfile,bamfileDonor,bamindexfileDonor,bamfileAcceptor,bamindexfileAcceptor,bamfileChimer,bamindexfileChimer,CutoffHom,CutoffHet){

    # bedfile="/home/cocca/analyses/michelangelo/chimerismo/4HotSpot_IAD132502_181_40amp_mod1.bed"
    # bedfile="/home/cocca/analyses/michelangelo/chimerismo/4HotSpot_IAD132502_181_40amp_mod1_AMPL7165264376_sorted.bed"
    # genfile="/home/cocca/analyses/michelangelo/chimerismo/5HotSpot_IAD132502_181_40amp.txt"
    # bamfileDonor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/Don.bam"
    # bamindexfileDonor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/Don.bai"
    # bamfileAcceptor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/RIC_40amp.bam"
    # bamindexfileAcceptor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/RIC_40amp.bai"
    # bamfileChimer="/home/cocca/analyses/michelangelo/chimerismo/20AMP/chimera_01.bam"
    # bamindexfileChimer="/home/cocca/analyses/michelangelo/chimerismo/20AMP/chimera_01.bai"

    ###################################################################################
    # Set donor == chimera
    # bamfileDonor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/Don.bam"
    # bamindexfileDonor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/Don.bai"
    # bamfileAcceptor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/RIC_40amp.bam"
    # bamindexfileAcceptor="/home/cocca/analyses/michelangelo/chimerismo/20AMP/RIC_40amp.bai"
    # bamfileChimer="/home/cocca/analyses/michelangelo/chimerismo/20AMP/Don.bam"
    # bamindexfileChimer="/home/cocca/analyses/michelangelo/chimerismo/20AMP/Don.bai"
    ###################################################################################
    
    # CutoffHom=0.9
    # CutoffHet=c(0.3,0.6)


    CountsDonor<-CountsBases(bedfile,bamfileDonor,bamindexfileDonor)
    CountsDonor<-CountsDonor[,1:4]
    CountsAccetor<-CountsBases(bedfile,bamfileAcceptor,bamindexfileAcceptor)
    CountsAccetor<-CountsAccetor[,1:4]
    CountsChimers<-CountsBases(bedfile,bamfileChimer,bamindexfileChimer)
    CountsChimers<-CountsChimers[,1:4]
    SNPs <- read.table(genfile,header=T,row.names=1)

    #Remove snps without information in Chimers/Donor/Accetor
    CountsDonor<-CountsDonor[(apply(CountsDonor,1,sum)>0)&(apply(CountsChimers,1,sum)>0),]
    CountsAccetor<-CountsAccetor[(apply(CountsAccetor,1,sum)>0)&(apply(CountsChimers,1,sum)>0),]
    CountsDonor<-CountsDonor[rownames(CountsDonor)%in%rownames(CountsAccetor),]
    CountsAccetor<-CountsAccetor[rownames(CountsAccetor)%in%rownames(CountsDonor),]
    CountsChimers<-CountsChimers[rownames(CountsChimers)%in%rownames(CountsDonor),]
    SNPs<-SNPs[rownames(SNPs)%in%rownames(CountsChimers),]

   if(dim(SNPs)[1] > 0) {
    print("We have some SNPs in common!!")
    # we need to have some snps in common between all the samples!
   #########################################################################
    #Select expected SNPs and Genotype
    DonorStep1<-CountsDonor[,1:4]
    colnames(DonorStep1)<-c("A","B","Min","Max")
    for(i in 1:dim(CountsDonor)[1]){
        DonorStep1[i,1:2]<-names(CountsDonor[i,order(CountsDonor[i,])[4:3]])
        DonorStep1[i,3:4]<-CountsDonor[i,order(CountsDonor[i,])[3:4]]
        if(DonorStep1[i,"Max"]>(sum(DonorStep1[i,3:4])*CutoffHom)){DonorStep1[i,"B"]<-"other";DonorStep1[i,"Min"]<-0}
        if(DonorStep1[i,"Max"]>(sum(DonorStep1[i,3:4])*CutoffHet[2])&DonorStep1[i,"Max"]<(sum(DonorStep1[i,3:4])*(CutoffHom-0.01))){DonorStep1[i,1:2]<-c("other","other");DonorStep1[i,1:2]<-c(0,0)}
    }

    AccettorStep1<-CountsAccetor[,1:4]
    colnames(AccettorStep1)<-c("A","B","Min","Max")
    for(i in 1:dim(CountsAccetor)[1]){
        AccettorStep1[i,1:2]<-names(CountsAccetor[i,order(CountsAccetor[i,])[4:3]])
        AccettorStep1[i,3:4]<-CountsAccetor[i,order(CountsAccetor[i,])[3:4]]
        if(AccettorStep1[i,"Max"]>(sum(AccettorStep1[i,3:4])*CutoffHom)){AccettorStep1[i,"B"]<-"other";AccettorStep1[i,"Min"]<-0}
        if(AccettorStep1[i,"Max"]>(sum(AccettorStep1[i,3:4])*CutoffHet[2])&AccettorStep1[i,"Max"]<(sum(AccettorStep1[i,3:4])*(CutoffHom-0.01))){AccettorStep1[i,1:2]<-c("other","other");AccettorStep1[i,1:2]<-c(0,0)}
    }

    ChimersStep1<-CountsChimers[,1:4]
    colnames(ChimersStep1)<-c("A","B","Min","Max")
    for(i in 1:dim(CountsChimers)[1]){
        ChimersStep1[i,1:2]<-names(CountsChimers[i,order(CountsChimers[i,])[4:3]])
        ChimersStep1[i,3:4]<-CountsChimers[i,order(CountsChimers[i,])[3:4]]
        if(ChimersStep1[i,"Max"]>(sum(ChimersStep1[i,3:4])*CutoffHom)){ChimersStep1[i,"B"]<-"other";ChimersStep1[i,"Min"]<-0}
        if(ChimersStep1[i,"Max"]>(sum(ChimersStep1[i,3:4])*CutoffHet[2])&ChimersStep1[i,"Max"]<(sum(ChimersStep1[i,3:4])*(CutoffHom-0.01))){ChimersStep1[i,1:2]<-c("other","other");ChimersStep1[i,1:2]<-c(0,0)}
    }

###################################################################################################################
#Select informative snps for the donor and the acceptor
    DonorInformative<-DonorStep1
    DonorInformative[,1:4]<-0
    select<-rownames(DonorInformative)
    for(i in 1:dim(DonorInformative)[1]){
        #select snp if its genotypes are hom/hom for different alleles in donor and acceptor -> AA vs BB 
        if((DonorStep1[select[i],2]=="other")&(AccettorStep1[select[i],2]=="other")&(DonorStep1[select[i],1]!=AccettorStep1[select[i],1])){DonorInformative[select[i],]<-DonorStep1[select[i],]}
        #select snp if its genotypes are het in donor and hom for acceptor -> AB vs BB|AA 
        if((DonorStep1[select[i],2]!="other")&(AccettorStep1[select[i],2]=="other")){DonorInformative[select[i],]<-DonorStep1[select[i],]}
    }
    DonorInformative<-DonorInformative[DonorInformative[,4]!="0",]
    rm(select)

    AccettorInformative<-AccettorStep1
    AccettorInformative[,1:4]<-0
    select<-rownames(AccettorInformative)
    for(i in 1:dim(AccettorInformative)[1]){
        #select snp if its genotypes are hom/hom for different alleles in donor and acceptor -> AA vs BB 
        if((AccettorStep1[select[i],2]=="other")&(DonorStep1[select[i],2]=="other")&(AccettorStep1[select[i],1]!=DonorStep1[select[i],1])){AccettorInformative[select[i],]<-AccettorStep1[select[i],]}
        #select snp if its genotypes are het in acceptor and hom for donor -> AB vs BB|AA 
        if((AccettorStep1[select[i],2]!="other")&(DonorStep1[select[i],2]=="other")){AccettorInformative[select[i],]<-AccettorStep1[select[i],]}
    }
    AccettorInformative<-AccettorInformative[AccettorInformative[,4]!="0",]
    rm(select)
###################################################################################################################

###################################################################################################################
#Calculate snps for the chimera: quantification for single SNP
# We need to work in block-mode: check if there are differences between chimera snps and acceptor in this block
# if there is only one difference, we nedd to perform a correction on the different genotype
    if(dim(AccettorInformative)[1] >= 3){
        print("We have more than 3 informative SNPs in this block!!")

        #we also want to discard those blocks where the accettor,the donor and the chimera are not always hom
        if (isTRUE(unique(unique(AccettorInformative$B)=="other"))&isTRUE(unique(unique(DonorInformative$B)=="other"))&isTRUE(unique(unique(ChimersStep1[rownames(AccettorInformative),"B"])=="other"))){


        # then calculate for each snp the acceptor reads count and divide by the total numbers of reads in this block
        #set a counter for to check block differences
        block_diff <- 0

        #count different snps in the current block and save the different ids to eventually correct
        for(i in 1:dim(AccettorInformative)[1]){
            #chek hom case first chimera and accettor are homo: we want to check if are the same or different
            if (ChimersStep1[rownames(AccettorInformative[i,]),"B"]=="other"&AccettorInformative[i,"B"]=="other") {
                if(ChimersStep1[rownames(AccettorInformative[i,]),"A"]!=AccettorInformative[i,"A"]){
                    # we have a mismatch, so we need to increase the block_diff counter
                    block_diff <- block_diff + 1
                    #and we need to save the mismatched snp to eventually correct it
                    to_correct <- CountsChimers[rownames(AccettorInformative[i,]),]
                }
            } else {

                if (ChimersStep1[rownames(AccettorInformative[i,]),"B"]!="other"&ChimersStep1[rownames(AccettorInformative[i,]),"B"]==AccettorInformative[i,"B"]&ChimersStep1[rownames(AccettorInformative[i,]),"A"]==AccettorInformative[i,"A"]) {
                #we have a match so no need to increase the block different
                } else {
                    if (ChimersStep1[rownames(AccettorInformative[i,]),"B"]!="other"&(ChimersStep1[rownames(AccettorInformative[i,]),"B"]==AccettorInformative[i,"A"]&ChimersStep1[rownames(AccettorInformative[i,]),"A"]==AccettorInformative[i,"B"])) {
                    #we still have a match so no need to increase the block different

                    } else {
                        #we have a final mismatch, so we need to increase the mismatch block count and save the variant
                        block_diff <- block_diff + 1
                        #and we need to save the mismatched snp to eventually correct it
                        to_correct <- CountsChimers[rownames(AccettorInformative[i,]),]
                    }
                }
            }

        }

        #We can do a correction only if we have 1 difference between Chimera and Acceptor
        if (block_diff == 1) {
            print("We can attemp a correction!!")
            CountsChimersCorr <- CountsChimers
            #correct the read, based on the difference we got:
            to_correct_snp <- rownames(to_correct)
            if (ChimersStep1[to_correct_snp,"B"]=="other"&AccettorInformative[to_correct_snp,"B"]!="other"){
                print("Correcting hom chimera to het!")
                #it means we have a hom site that has to be set to het
                #we need to know the other acceptor allele
                if (AccettorInformative[to_correct_snp,"A"]==ChimersStep1[to_correct_snp,"A"]) {
                    # in this case we need to calculathe the missing reads in the othe allele, and assign a new total
                    correction <- abs(CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"B"]] - (CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"A"]])/2)
                    CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"A"]] <- CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"A"]]/2 + 0.5
                    CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"B"]] <- correction + CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"B"]] - 0.5

                } else {
                 if (AccettorInformative[to_correct_snp,"B"]==ChimersStep1[to_correct_snp,"A"]) {
                    # in this case we need to calculathe the missing reads in the othe allele, and assign a new total
                    correction <- abs(CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"A"]] - (CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"B"]])/2)
                    CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"A"]] <- correction + CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"A"]] - 0.5
                    CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"B"]] <- CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"B"]]/2 + 0.5

                    }   
                }

            } else if (ChimersStep1[to_correct_snp,"B"]=="other"&AccettorInformative[to_correct_snp,"B"]=="other") {
                #it means we have a hom site that has to be set to the other hom site!!
                print("Correcting hom chimera to other hom!")
                #it means we have a hom site that has to be set to het
                #we need to know the other acceptor allele
                if (AccettorInformative[to_correct_snp,"A"]!=ChimersStep1[to_correct_snp,"A"]) {
                    # in this case we need to calculathe the missing reads in the othe allele, and assign a new total
                    correction <- CountsChimersCorr[to_correct_snp,ChimersStep1[to_correct_snp,"A"]]
                    CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"A"]] <- correction
                    CountsChimersCorr[to_correct_snp,ChimersStep1[to_correct_snp,"A"]] <- 0

                }   

            } else {
                # It means we have a het site to fix versus an hom site
                print("Correcting het chimera to hom!")
                if (AccettorInformative[to_correct_snp,"A"]==ChimersStep1[to_correct_snp,"A"]) {
                    # in this case we need to calculathe the missing reads in the othe allele, and assign a new total
                    # correction <- CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"A"]]
                    CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"A"]] <- 0
                    CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"B"]] <- CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"B"]]*2

                } else {
                 if (AccettorInformative[to_correct_snp,"B"]==ChimersStep1[to_correct_snp,"A"]) {
                    # in this case we need to calculathe the missing reads in the othe allele, and assign a new total
                    # correction <- CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"B"]]
                    CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"B"]] <- 0
                    CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"A"]] <- CountsChimersCorr[to_correct_snp,AccettorInformative[to_correct_snp,"A"]]*2
                    }   
                }

            }

        } else {
            #no can do any correction!!
            CountsChimersCorr <- CountsChimers

        }

        print("Let's quantify!!!")
        #we have to quantify, BY BLOCK the total Acceptor reads
        FromAccetor<-matrix(ncol=3,nrow=dim(AccettorInformative)[1],data=NA)
        rownames(FromAccetor)<-rownames(AccettorInformative)
        #quantification of acceptor in chimera in the block
        for(i in 1:dim(AccettorInformative)[1]){

            if(AccettorInformative[i,"Min"]==0&AccettorInformative[i,"B"]=="other"){
                FromAccetor[i,1]<-sum(CountsChimersCorr[rownames(AccettorInformative[i,]),AccettorInformative[i,"A"]])
                FromAccetor[i,2]<-sum(CountsChimersCorr[rownames(AccettorInformative[i,]),])
            }
            if(AccettorInformative[i,"B"]!="other"&AccettorInformative[i,"A"]%in%DonorStep1[rownames(AccettorInformative)[i],c("A","B")]){
                # FromAccetor[i,1]<-sum(CountsChimersCorr[rownames(AccettorInformative[i,]),AccettorInformative[i,"B"]]*2)
                FromAccetor[i,1]<-sum(CountsChimersCorr[rownames(AccettorInformative[i,]),AccettorInformative[i,"B"]]*2)
                FromAccetor[i,2]<-sum(CountsChimersCorr[rownames(AccettorInformative[i,]),])
            }
            if(AccettorInformative[i,"B"]!="other"&AccettorInformative[i,"B"]%in%DonorStep1[rownames(AccettorInformative)[i],c("A","B")]){
                # FromAccetor[i,1]<-sum(CountsChimersCorr[rownames(AccettorInformative[i,]),AccettorInformative[i,"A"]]*2)
                FromAccetor[i,1]<-sum(CountsChimersCorr[rownames(AccettorInformative[i,]),AccettorInformative[i,"A"]]*2)
                FromAccetor[i,2]<-sum(CountsChimersCorr[rownames(AccettorInformative[i,]),])
            }

        }

        FromAccetor[,3] <- FromAccetor[,1]/FromAccetor[,2]
        FromAccetor_block <- sum(FromAccetor[,1])/sum(FromAccetor[,2])

        # for(i in 1:dim(AccettorInformative)[1]){

        #     if(AccettorInformative[i,"Min"]==0&AccettorInformative[i,"B"]=="other"){
        #         FromAccetor[i,1]<-sum(CountsChimers[rownames(AccettorInformative[i,]),AccettorInformative[i,"A"]])/sum(CountsChimers[rownames(AccettorInformative[i,]),])
        #     }
        #     if(AccettorInformative[i,"B"]!="other"&AccettorInformative[i,"A"]%in%DonorStep1[rownames(AccettorInformative)[i],c("A","B")]){
        #         FromAccetor[i,1]<-sum(CountsChimers[rownames(AccettorInformative[i,]),AccettorInformative[i,"B"]]*2)/sum(CountsChimers[rownames(AccettorInformative[i,]),])
        #     }
        #     if(AccettorInformative[i,"B"]!="other"&AccettorInformative[i,"B"]%in%DonorStep1[rownames(AccettorInformative)[i],c("A","B")]){
        #         FromAccetor[i,1]<-sum(CountsChimers[rownames(AccettorInformative[i,]),AccettorInformative[i,"A"]]*2)/sum(CountsChimers[rownames(AccettorInformative[i,]),])
        #     }

        # }
        } else {
            # discard block!
            print("Block discarded!!!Not all HOM!!")
            FromAccetor_block <- NA
        }

        } else {
            # discard block!
            print("Block discarded!!!")
            FromAccetor_block <- NA
        }
   } else {
        print("No way!! Block discarded again!!!")
        FromAccetor_block <- NA
    
   }

        return(FromAccetor_block)
    # Media<-round(mean(FromAccetor[!is.na(FromAccetor)]),3)
    # DeviazioneStandard<-round(sd(FromAccetor[!is.na(FromAccetor)]),3)
    # Error<-round(qnorm(0.975)*DeviazioneStandard/sqrt(length(FromAccetor[!is.na(FromAccetor)])),3)
    # Design<-read.table(file=bedfile,sep="\t")
    # rownames(Design)<-Design[,4]
    # Design<-Design[rownames(Design)%in%rownames(FromAccetor)[!is.na(FromAccetor)],7:8]
    # Result<-paste("Predicted Recipient mean",Media*100,"%; Number of SNPs=",length(FromAccetor[!is.na(FromAccetor)]),";CI",(Media+Error)*100,"-",(Media-Error)*100,sep=" ")
    # ResultSummary<-List(Mean=Media,SD=DeviazioneStandard,Err=Error,Summary=Result,AcceptorTable=FromAccetor[!is.na(FromAccetor)],SNPS=rownames(FromAccetor)[!is.na(FromAccetor)],Design=Design)
    # return(ResultSummary)




}


CountsBases<-function(bedfile,bamfile,bamindexfile){
    # bedfile
    # bamfile=bamfileDonor
    # bamindexfile=bamindexfileDonor

    data <- read.table(bedfile,header=F)
    colnames(data) <- c('chr','start','end','id','score','strand')
    bed <- with(data, GRanges(chr, IRanges(end,width=1), strand="*", score, id=id))
    bamfile <- BamFile(bamfile,index=bamindexfile)
    seqinfo(bed) <- merge(seqinfo(bed), seqinfo(bamfile))
    seqlevels(bed) <- seqlevelsInUse(bed)
    which <- as(seqinfo(bed), "GRanges")
    flag <- scanBamFlag(isNotPassingQualityControls=NA,isDuplicate=NA)
    param <- ScanBamParam(flag=flag, what=c("seq", "qual"), which=which)
    # gal <- readGAlignmentsFromBam(bamfile, param=param)
    gal <- readGAlignments(bamfile, param=param)
    seqlevels(gal) <- seqlevels(bed)
    qseq <- mcols(gal)$seq
    qual <- mcols(gal)$qual
    nucl_piles <- pileLettersAt(qseq, seqnames(gal), start(gal), cigar(gal), bed)
    qual_piles <- pileLettersAt(qual, seqnames(gal), start(gal), cigar(gal), bed)
    mcols(bed)$nucl_piles <- nucl_piles
    mcols(bed)$qual_piles <- qual_piles
    #results<-data.frame(id=data$id,alphabetFrequency(nucl_piles, baseOnly=TRUE),alphabetFrequency(nucl_piles, baseOnly=TRUE)/apply(alphabetFrequency(nucl_piles, baseOnly=TRUE),1,sum))
    results<-data.frame(alphabetFrequency(nucl_piles, baseOnly=TRUE))
    rownames(results)<-data$id
    #write.table(results,file=paste(bamindexfile,".txt",sep=""),col.names=TRUE,row.names=TRUE,sep="\t")
    return(results)
}

DataToPlot<-function(Result){
    data.summary <-data.frame(treatment=as.factor("MRD"),mean=unlist(Result)[1],n=length(unlist(Result)[1]),sd=sd(unlist(Result)[5]))
    data.summary$sem <- data.summary$sd/sqrt(data.summary$n)
    data.summary$me  <- as.numeric(qt(1-0.01/2, df=data.summary$n)*data.summary$sem)
    return(data.summary)
}

Parameters<-function(input){
    InputData<-rbind(input$fileBed,input$fileGen,input$fileD,input$fileDi,input$fileA,input$fileAi,input$fileC,input$fileCi,input$Hom,input$Het[1],input$Het[2])
    rownames(InputData)<-c("Bed","Gen","Donor Bam","Donor Bai","Recipient Bam","Recipient Bai","Chimer BAM","Chimer Bai","Threshold_Hom","Threshold_Het_Mim","Threshold_Het_Max")
    colnames(InputData)<-"Parameters"
    return(InputData)
}

SNPsTable<-function(Result){
    data.summary <-data.frame(round(Result[[5]]*100,1),Result[[7]])
    rownames(data.summary)<-Result[[6]]
    colnames(data.summary)<-c("Values","Chr.Ampl","Type")
    return(data.summary)
}

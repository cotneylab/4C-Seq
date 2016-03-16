
#parameters
args<-commandArgs(TRUE)
parameterfile<-args[1] #
source(parameterfile)

#Code starts
source("/home/jy344/Programs/Scripts/4CSeq/fourcseq_revised_code.R")
library(FourCSeq)

#description of the fourC viewpoint
exptData <- SimpleList(projectPath = projectname,
	fragmentDir = "re_fragments",
	referenceGenomeFile = referenceGenomeFile,
	reSequence1 = reSequence1,
	reSequence2 = reSequence2,
	primerFile = primerFile,
	bamFilePath = bamFilePath)

#description of the exp
colData <- DataFrame(viewpoint = firstprimer,
condition = factor(conditions,levels = unique(conditions)),
replicate = replicates,
bamFile = bamfiles,
sequencingPrimer="first")

#create FourC object
fc <- FourC(colData, exptData)

#in silico digest of the genome
#default, remove fragments<20, removefragments that do not contain a cutting site of the second restriction enzyme
fc <- addFragments(fc)

#find viewpoint
findViewpointFragments(fc)
fc <- addViewpointFrags(fc)

#Count the reads for fragments #need to trim 4 bp for the restriction site
fc<-countFragmentOverlaps_edi(fc, trim=4, minMapq=-1) #need to trim 4bp for GATC or CATG

#No filter for direction of reads
fc <- combineFragEnds(fc,filter = FALSE)

#Smooth the couts
fc <- smoothCounts(fc)

#z-score
#default parameter, remove counts<40, remove counts till local minimum
fcf <- getZScores(fc) 
zScore <- assay(fcf, "zScore")

#use z-score threshold of 3, and fdr 0.01
fcf <- addPeaks(fcf, zScoreThresh=3, fdrThresh=0.1) #change here!!!

#write the results
results<-cbind(as.vector((rowData(fcf)@seqnames)),start(rowData(fcf)@ranges),end(rowData(fcf)@ranges),assay(fcf,"counts"),assay(fcf,"zScore"),assay(fcf,"pValue"),
assay(fcf,"pAdjusted"),assay(fcf,"peaks"))

repnames<-colnames(assay(fcf,"counts"))

colnames(results)<-c("chromsome","start","end",paste(c("Counts"),repnames),paste(c("zScore"),repnames),paste(c("pValue"),repnames),paste(c("BH pValue"),repnames),paste(c("Peaks"),repnames))

write.table(results,file=resultfile,col.names=NA,sep="\t",quote=F)


#scatter plot of the first two replicates
pdf(scatterplotfile)
plotScatter(fc[,c(colnames(fc)[1],colnames(fc)[2])],xlab=colnames(fc)[1], ylab=colnames(fc)[2], asp=1)
dev.off()

#plot fit figure
for(num in 1:ncol(fcf)) {
	fitplotfile.rep<-sub(".pdf",paste("_",num,".pdf",sep=""),fitplotfile)
	pdf(fitplotfile.rep,width=12, height=4)
	par(mfrow=c(1,3),cex=0.6)
	plotFits(fcf[,num])
	dev.off()
}


#plot zscores
pdf(peakplotfile)
plotZScores(fcf)
dev.off()

##FilterInteractions.R: Take the results.txt from the 4CSeq pipeline and 
# outputs the interactions as bed files that pass the input thresholds
#inputs = {ResultFile} {# Conditions} {# Replicates} {Z Cutoff} {P Cutoff} {Q Cutoff}

##Setup the important stuff
#Input Parameters
args<-commandArgs(TRUE)

numConds <- as.numeric(args[2])      #Number of different conditions per viewpoint
numReps  <- as.numeric(args[3])      #Number of replicates per conditions
cutoff.Z <- as.numeric(args[4])
cutoff.P <- as.numeric(args[5])
cutoff.Q <- as.numeric(args[6])

#Import necessary packages

#Input file
result <- read.table(args[1], sep="\t", header=TRUE, row.names=1)

#Create results folder and switch to it
dir.create(file.path(getwd(), "results"), showWarnings = FALSE)
setwd(file.path(getwd(), "results"))

##Calculate interactions that pass
#For the number of input conditions...
for(currCond in 0:(numConds-1)){                    #Starts with zero!!
   #Set columns with corresponding data
   condCols <- c(4+(0*numConds)+(currCond*numReps), 5+(0*numConds)+(currCond*numReps), #Counts
                 4+(2*numConds)+(currCond*numReps), 5+(2*numConds)+(currCond*numReps), #Z-Score
                 4+(4*numConds)+(currCond*numReps), 5+(4*numConds)+(currCond*numReps), #p-Value
                 4+(6*numConds)+(currCond*numReps), 5+(6*numConds)+(currCond*numReps), #Q-Score
                 4+(8*numConds)+(currCond*numReps), 5+(8*numConds)+(currCond*numReps)) #Peaks?
   
   #Apply cutoffs to p-value and q-score columns; save rows that pass
   ii <- which(
               #Z-Score Cutoff
               #result[,condCols[3]] <= cutoff.P & result[,condCols[4]] <= cutoff.P &
                  
               #p-Value Cutoff
               result[,condCols[5]] <= cutoff.P & result[,condCols[6]] <= cutoff.P &
               
               #Q-Score Cutoff
               result[,condCols[7]] <= cutoff.Q & result[,condCols[8]] <= cutoff.Q)

   #Create bed file of interactions that pass and add "chr#:start-end" column
   preBed <- result[ii,1:3]
   if(length(preBed[,1]) != 0){
      bed <- cbind(preBed, paste(preBed[,1], ":", preBed[,2], "-", preBed[,3], sep=''))}
   
    #Create positive results file with fragment info, p-val, q-val, and counts
   pass <- cbind(preBed, result[ii,condCols[c(1:2, 5:8)]]) 
   
    #Check to see if the columns are correct
   if(substr(colnames(result)[condCols[3]], 0, 6) != "zScore" |
      substr(colnames(result)[condCols[4]], 0, 6) != "zScore" |
      substr(colnames(result)[condCols[5]], 0, 6) != "pValue" |
      substr(colnames(result)[condCols[6]], 0, 6) != "pValue" |
      substr(colnames(result)[condCols[7]], 0, 9) != "BH.pValue" |
      substr(colnames(result)[condCols[8]], 0, 9) != "BH.pValue")
      {print("Columns don't match!"); ...}
   
   ##Output
   #Extract output file name from 1st column name (Assumes "Counts.viewpoint_cond" is first)
   currName <- substr(colnames(result)[condCols][1], 8, nchar(colnames(result)[condCols][1])-2)
   fileName <- paste(currName, "_Z", cutoff.Z,"_p", cutoff.P, "_Q", cutoff.Q, sep='')
  
   #Assign output file names
   output.bed <- paste(fileName, ".bed", sep='')
   output.txt <- paste(fileName, ".txt", sep=''); 
   
   #Output files (bed & positive results)
   #If there are results, output the bed; otherwise, output empty file with same name
   if(length(preBed[,1] != 0)){
      write.table(bed, output.bed, sep = "\t", col.names=FALSE, row.names=FALSE, quote=FALSE)}
   else{write.table(NULL, output.bed, quote=FALSE)}
   
   #Output expanded interaction region data (positive results)
   write.table(pass, output.txt, sep = "\t", row.names=FALSE, quote=FALSE)
}
1) you need to install the FourCSeq package from bioconductor:
source("http://bioconductor.org/biocLite.R")
biocLite("FourCSeq”)

2) you need to put these two scripts in the folder where you want to do the analysis.
run_fourcseq_workflow.R
fourcseq_revised_code.R

3)In this folder, you need to create two folders:
bam
primers

You need to put the bam files from 4C-seq in the bam folders, and primer information in the primers folder.
I also attached an example of the primer file (Wnt7aPromoter_3_primers.fas).

4) Edit the parameters in run_fourcseq_workflow.R

5) run the script in the shell
$Rscript run_fourcseq_workflow.R parameters.txt

You will find the peak calling results and figures in the same folder.

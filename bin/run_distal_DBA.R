args=commandArgs(T)
otu_file = args[1]
sample_file = args[2]
dlimits = args[3] 

source("../lib/distal_DBA.R")

distal_DBA(otu_file, sample_file,as.numeric(dlimits)) #0.0001
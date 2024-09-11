requireNamespace("Biostrings", "entropy", "BiocManager")
AMIP<-function(fasta_file,n1=1,n2=4){
  x=Biostrings::readDNAStringSet(fasta_file)
  #calculating frequency of occurence of nucleotides k bases apart
  AMI_fun<-function(x){
    y=Biostrings::oligonucleotideFrequency(x, width=1, step=1)
    z<-matrix(nrow=(n2-n1)+1,ncol=4)
    F<-list()
    length(F)<-(n2-n1)+1
    R<-numeric((n2-n1)+1)
    for (i in 1:((n2-n1)+1)){
      z[i,]=Biostrings::oligonucleotideFrequency(x, width=1, step=i+n1-1)
      F[[i]]=rbind(y,z[i,])
      R[i]=entropy::mi.plugin(F[[i]])
    }
    R=round(R,4)
    mean_AMI<-round(mean(R),4)
    AMI<-list( mean_AMI)
    return(AMI)
  }
  res<-lapply(x,AMI_fun)
  ress= data.frame(res)
  ress = t(ress)
  row.names(ress) <- names(x)
  colnames(ress) <- c("Mean Of AMIP")
  return(ress)
}


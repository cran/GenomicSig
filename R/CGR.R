requireNamespace("seqinr","kaos")
CGR <- function(data){
  
  A=function(data){
    sequence <-unlist(strsplit(data,""))
    sequence.cgr = kaos::cgr(sequence,res = 100)
    rplot_cgr <- kaos::cgr.plot(sequence.cgr, mode = "points")
    cgrvector=kaos::vectorize(sequence.cgr)
    result=list(sequencecgr=sequence.cgr,allplots=rplot_cgr, output_of_cgr_vector=cgrvector)
    return(result)
  }
  B=lapply(data,A)
  return(B)
}

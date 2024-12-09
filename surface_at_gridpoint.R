### ADAPTED FROM : Schmiedel and Lehner
### calculate fitness quantiles for nearest neighbours of gridpoint
#subfunction for call_epistasis class of functions
surface_at_gridpoint = function(i) {
  require(data.table)
  A = unlist(subDT[,.(D=sqrt((s_1-Fq[i,s_1])^2 +(s_2-Fq[i,s_2])^2),fitness_norm)][D < quantile(D,probs=span,na.rm=T),
                                                                                                      .(quantile(fitness_norm,p=c(0.5,Q,1-Q),na.rm=T))])
  return(A)
}

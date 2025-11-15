BiomassData=function(start_year,end_year){
  
  # load functions and data used in the following
  source("Z:/Rdata/buildcapelindata.R")
  capTab=buildcapelindata()

  source("Z:/Rdata/ADMButils.s")
  
  # -------------------------------------------------------------------------------------------
  # create a file called "Biomass.dat"
  meanlength=as.vector(data.matrix(capTab[[1]]["Meanlength(cm)"]))
  start_index=start_year-1971
  end_index=end_year-1971
  
  # Bl is the reported Biomass of capelin, in all the age classes, start_year to end_year
  # It has the following strucutre: Bl(year,l,a)=Bl((year-1)*length_l+l,a)
  # dim(Bl): length(meanlength)*(end_year-start_year) x 1
  Bl=data.matrix(capTab[[start_index]][,7])
  for(i in (start_index+1):end_index){
    # rbind takes to matrices A and B and outputs one matrix (A,B)'
    Bl=rbind(Bl,data.matrix(capTab[[i]][,7]))
  }
  
  # B is the total Biomass of reported capelin for one age class and one year:
  # B(year,a)=sum_l Bl(year,l,a)
  B=sum(Bl[1:length(meanlength)])
  for(i in 2:(1+(end_index-start_index))){        #(start_index+1):end_index
    tmp=sum(Bl[((i-1)*length(meanlength)+1):(i*length(meanlength))])
    B=rbind(B,tmp)
  }
  

  dat_write("Biomass.dat",list(start_year=start_year,end_year=end_year,Bl=Bl,B=B))
  return(list(start_year=start_year,end_year=end_year,Bl=Bl,B=B))
}
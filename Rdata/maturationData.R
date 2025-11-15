# maturationData is the main function to read data for the maturation process and save it in a '*.dat'-file
# INPUT: start_year, end_year, min_age, max_age

# !'ADMButils.s' has to be in the workspace

maturationData=function(start_year,end_year,min_age,max_age){

# load functions and data used in the following
source("Z:/Rdata/buildcapelindata.R")
capTab=buildcapelindata()
source("Z:/Rdata/catchMaturation.R")
catchM=catchMaturation()
source("Z:/Rdata/ADMButils.s")

# -------------------------------------------------------------------------------------------
# create a file called "maturation.dat", including the data necessary for parameter estimation

meanlength=as.vector(data.matrix(capTab[[1]]["Meanlength(cm)"]))
start_index=start_year-1971
end_index=end_year-1971

# Nl is the reported number of capelin, in age classes min_age to max_age, start_year to end_year
# It has the following strucutre: Nl(year,l,a)=Nl((year-1)*length_l+l,a)
# dim(Nl): length(meanlength)*(end_year-start_year) x (max_age-min_age)
Nl=data.matrix(capTab[[start_index]][,min_age:max_age])
for(i in (start_index+1):end_index){
  # rbind takes to matrices A and B and outputs one matrix (A,B)'
  Nl=rbind(Nl,data.matrix(capTab[[i]][,min_age:max_age]))
}

# N is the total number of reported capelin for one age class and one year:
# N(year,a)=sum_l Nl(year,l,a)
N=colSums(Nl[1:length(meanlength),])
for(i in 2:(1+(end_index-start_index))){        #(start_index+1):end_index
  tmp=colSums(Nl[((i-1)*length(meanlength)+1):(i*length(meanlength)),])
  N=rbind(N,tmp)
}

# catchM are the monthly catches in years start_year to end_year; for all age classes and months october-september
# It has the following strucutre: C(year,a,m)=C((year-1)*(min_age-max_age)+(a-min_age),m) 
catchM=catchM[((start_index-1)*5+1):(end_index*5),]
C=catchM[min_age:max_age,]
for(i in 1:((dim(catchM)[1]/5)-1)){
  C=rbind(C, catchM[(i*5+min_age):(i*5+max_age),])
}

# save dimensions, as well
length_l=length(meanlength)
index_Nl=start_year+dim(Nl)[1]-1
index_C=start_year+dim(C)[1]-1

dat_write("maturation.dat",list(min_age=min_age,max_age=max_age,start_year=start_year,end_year=end_year,length_l=length_l,index_Nl=index_Nl,index_C=index_C,meanlength=meanlength,Nl=Nl,N=N,C=C))
return(list(min_age=min_age,max_age=max_age,start_year=start_year,end_year=end_year,length_l=length_l,index_Nl=index_Nl,index_C=index_C,meanlength=meanlength,Nl=Nl,N=N,C=C))
}

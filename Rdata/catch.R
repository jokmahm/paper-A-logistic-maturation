# The function catch_xls_to_dat takes a *.xls file and outputs a list containing:
# a vector "scale" describing the scales of the columns in the data frames
# a data frame "loddedata" describing the catches
# an additional data frame of catches in 1976.

catch_xls_to_dat=function(xls_file){
# library including 'read.xls'
require(gdata)
# read from excel-file
loddedata=read.xls(xls_file,perl = "Z:/Rdata/strawberry/perl/bin/perl.exe",header=TRUE,blank.lines.skip=TRUE,as.is=TRUE,verbose=TRUE,skip=4)

# delete empty columns
for(i in dim(loddedata)[2]:1){
  if(all(is.na(loddedata[,i]))){
    loddedata=loddedata[,-i]}}

# initializations
scale=loddedata[1,2:23]
years=c(1972:2019)

# Rephrase column names
tmp=colnames(loddedata)[1:23]
tmp[2]="Total.catch.spring<14cm"
tmp[3]="Total.catch.spring<14cm.1"
tmp[4]="Total.catch.spring>14cm"
tmp[5]="Total.catch.spring>14cm.1"
tmp[10]="Total.catch.Oct.Dec<14cm"
tmp[11]="Total.catch.Oct.Dec<14cm.1"
tmp[12]="Total.catch.Oct.Dec>14cm"
tmp[13]="Total.catch.Oct.Dec>14cm.1"

# remove all lines which are not necessary
for(i in (dim(loddedata)[1]):1){
  if(loddedata[i,1]!="1" && loddedata[i,1]!="2" && loddedata[i,1]!="3" && loddedata[i,1]!="4" && loddedata[i,1]!="5"){
    loddedata=loddedata[-i,]
  }
}
loddedata=as.data.frame(sapply(loddedata,as.numeric))

loddedata=loddedata[,1:23]
colnames(loddedata)=tmp

# add years
yearss=c()
for(i in 1:dim(loddedata)[1]){
  yearss[i]=years[((i-1)%/%5)+1]}
loddedata[,"Year"]=yearss
#add_data[,"Year"]=1976

# build matrix and .dat file for each table 1-15
liste=list(scale,loddedata)
return(liste)
}


# Call: e.g. catch=catch_xls_to_dat("./Lodde/Data/Fangstdata/CapCatchUpdated.xls")

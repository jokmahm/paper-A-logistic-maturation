# Function to generate capelin data in form of a data.frame
# OUTPUT: data.frame containing stock information for several years, as specified in the first line.  For each year, estimates for stock sizes and properties of the stock are stored. The elements of the list of data.frames are the data.frames of one year. For each year,  the row represents the length and the column represents the age or measure. 
# PURPOSE:  To store information of all years considered in one data.frame. 
# USES:  captab

buildcapelindata=function(){

# specification of the years considered
years=c(1972:2019)

# load function captab_xls_to_dat
source("Z:/Rdata/captab.R")

# create a list containing all files to be read
files=list("Z:/Rdata/Bestandsdata/CapTab20111.csv","Z:/Rdata/Bestandsdata/CapTab20112.csv","Z:/Rdata/Bestandsdata/CapTab20113.csv","Z:/Rdata/Bestandsdata/CapTab20114.csv","Z:/Rdata/Bestandsdata/CapTab20115.csv","Z:/Rdata/Bestandsdata/CapTab20116.csv","Z:/Rdata/Bestandsdata/CapTab20117.csv","Z:/Rdata/Bestandsdata/CapTab20118.csv","Z:/Rdata/Bestandsdata/CapTab20119.csv","Z:/Rdata/Bestandsdata/CapTab201110.csv","Z:/Rdata/Bestandsdata/CapTab201111.csv","Z:/Rdata/Bestandsdata/CapTab201112.csv","Z:/Rdata/Bestandsdata/CapTab201113.csv","Z:/Rdata/Bestandsdata/CapTab201114.csv","Z:/Rdata/Bestandsdata/CapTab201115.csv","Z:/Rdata/Bestandsdata/CapTab201116.csv","Z:/Rdata/Bestandsdata/CapTab201117.csv","Z:/Rdata/Bestandsdata/CapTab201118.csv","Z:/Rdata/Bestandsdata/CapTab201119.csv","Z:/Rdata/Bestandsdata/CapTab201120.csv","Z:/Rdata/Bestandsdata/CapTab201121.csv","Z:/Rdata/Bestandsdata/CapTab201122.csv","Z:/Rdata/Bestandsdata/CapTab201123.csv","Z:/Rdata/Bestandsdata/CapTab201124.csv","Z:/Rdata/Bestandsdata/CapTab201125.csv","Z:/Rdata/Bestandsdata/CapTab201126.csv","Z:/Rdata/Bestandsdata/CapTab201127.csv","Z:/Rdata/Bestandsdata/CapTab201128.csv","Z:/Rdata/Bestandsdata/CapTab201129.csv","Z:/Rdata/Bestandsdata/CapTab201130.csv","Z:/Rdata/Bestandsdata/CapTab201131.csv","Z:/Rdata/Bestandsdata/CapTab201132.csv","Z:/Rdata/Bestandsdata/CapTab201133.csv","Z:/Rdata/Bestandsdata/CapTab201134.csv","Z:/Rdata/Bestandsdata/CapTab201135.csv","Z:/Rdata/Bestandsdata/CapTab201136.csv","Z:/Rdata/Bestandsdata/CapTab201137.csv","Z:/Rdata/Bestandsdata/CapTab201138.csv","Z:/Rdata/Bestandsdata/CapTab201139.csv","Z:/Rdata/Bestandsdata/CapTab201140.csv","Z:/Rdata/Bestandsdata/CapTab201141.csv","Z:/Rdata/Bestandsdata/CapTab201142.csv","Z:/Rdata/Bestandsdata/CapTab201143.csv","Z:/Rdata/Bestandsdata/CapTab201144.csv","Z:/Rdata/Bestandsdata/CapTab201145.csv","Z:/Rdata/Bestandsdata/CapTab201146.csv","Z:/Rdata/Bestandsdata/CapTab201147.csv","Z:/Rdata/Bestandsdata/CapTab201148.csv")

# read the files
capTab=list()
for(i in 1:length(years)){
  capTab[[i]]=captab_xls_to_dat(files[[i]],years[i])
}

# replace NA by 0
for(i in 1:length(capTab)){
for(j in 1:dim(capTab[[i]])[1]){
for(k in 1:dim(capTab[[i]])[2]){
if(is.na(capTab[[i]][j,k])){
capTab[[i]][j,k]=0}}}}

# if wished, save the data.frame:
# save(capTab,file="./R/Lodde/CapTab2011.Rda")
return(capTab)
}


#----------------------------------------------------------------------------------------
#Example for use of the function:
#buildcapelincatch=function(){
#source("./R/Lodde/catch.R")
#catch=catch_xls_to_dat("./Lodde/Data/Fangstdata/CapCatchUpdated.xls")
#save(capTab,file="./R/Lodde/CapTab2011.Rda")
#}

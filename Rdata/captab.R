# PURPOSE: The function reads a single table containing stock information of one year.
# INPUT: path for the data describing the stock, named "capTab", e.g. './Lodde/Data/Bestandsdata/CapTab20111.csv' and year of measurement
# OUTPUT:  data.frame consisting of 10 columns and 32 rows. The first 5 columns contain estimates for stock sizes for capelin of age 1 to 5+, for a certain length. The latter is given by the row. Column 6 to 10 represent total stock sizes, biomass, mean weight, mean length and year of measurement, respectively.
#
# The *.xls file contains: 
# Stockestimates for fish of age 1 to 5+ (columns 4 to 9, respectively)
# and lengths 5.00-21.00 (rows 4 to 35).
# Summaries of the stock as sum and biomass per length group / age group,
# mean lengths and weights. 

captab_xls_to_dat=function(csv_file,year){
# library including 'read.xls'
require(gdata)
loddedata=read.csv(csv_file,header=TRUE,blank.lines.skip=TRUE,as.is=TRUE,skip=1,nrow=33)

# Delete obsolent rows and columns
loddedata=loddedata[-1,]
if(dim(loddedata)[2]==13){loddedata=loddedata[,-13]}
if(year==1982){loddedata=loddedata[,18:26]
}else{
loddedata=loddedata[,-3]
loddedata=loddedata[,-2]
loddedata=loddedata[,-1]}

# Rephrase column names
tmp=colnames(loddedata)
tmp[1]=1
tmp[2]=2
tmp[3]=3
tmp[4]=4
tmp[5]=5
tmp[6]="Sum(10e9)"
tmp[7]="Biomass(10e3t)"
tmp[8]="Meanweight(g)"
tmp[9]="Meanlength(cm)"
colnames(loddedata)=tmp

loddedata=as.data.frame(sapply(loddedata[,1:9],as.numeric))

# add a column containing the year
loddedata[,"Year"]=year

return(loddedata)
}


# Call: e.g. capTab=captab_xls_to_dat("./Lodde/Data/Bestandsdata/CapTab20111.csv",1972)
# save(capTab,file="CapTab2007.Rda")
# load("CapTab2007.Rda")
# csv_file="./Lodde/Data/Bestandsdata/CapTab2011.csv"

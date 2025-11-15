# This function computes the monthly catch data.
# IN: a table containing all catch information (obtained by function "catch.R")
# OUT: total catches for each month: catchM, dim(catchM)=(number of years * number of ages) x 12

# Explanation: The catch is reported for seasons January-March, August-September, October-December.
# The catch in one month is a constant proportion of the total catch in the season.
# The distribution key for January-March is [0.1,0.5,0.4],
# for August-September [0.8,0.2] and [0.5,0.3,0.2] for winter.
# For recent keys, see the ICES report.
# The year starts in october, i.e. month 1 is october, month 12 is september.


catchMaturation=function(){

source("Z:/Rdata/catch.R")
catch=catch_xls_to_dat("Z:/Rdata/Fangstdata/CapCatchUpdated.xls")

# Total catch per year: catchM=catch[[2]][,23]
catchSpring=catch[[2]][,2]
# keySpring=c(0.1,0.5,0.4)
catchAutumn=catch[[2]][,8]
# keyAutumn=c(0.8,0.2)
catchWinter=catch[[2]][,10]
# keyWinter=c(0.5,0.3,0.2)
fill=mat.or.vec(205,1)

catchM=cbind(0.5*catchWinter,0.3*catchWinter,0.2*catchWinter,0.1*catchSpring,0.5*catchSpring,0.4*catchSpring,fill,fill,fill,fill,0.8*catchAutumn,0.2*catchAutumn)

return(catchM)
}

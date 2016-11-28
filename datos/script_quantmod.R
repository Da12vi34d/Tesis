#
#	Consulta de datos financieros usando "quantmod"
#

# install.packages("quantmod")

rm(list=ls())

require("quantmod")

path.data <- "C:/Users/jmartineov/Documents/GitHub/Tesis/datos/"

#	Apple
getSymbols("AAPL",src="google")

#	Microsoft
getSymbols("MSFT",src="google")

#	Yahoo
getSymbols("YHOO",src="google")

#	Data 4 analysis (levels)
data4analysis_l <- cbind(AAPL$AAPL.Close,MSFT$MSFT.Close,YHOO$YHOO.Close)
colnames(data4analysis_l) <- c("AAPL","MSFT","YHOO")
write.csv(data4analysis_l, file = paste(path.data,"data4analysis_l.csv",sep=""), row.names = TRUE)

#	Data 4 analysis (returns)
data4analysis_r <- cbind(Delt(AAPL$AAPL.Close),Delt(MSFT$MSFT.Close),Delt(YHOO$YHOO.Close))
colnames(data4analysis_r) <- c("AAPL","MSFT","YHOO")
write.csv(data4analysis_r, file = paste(path.data,"data4analysis_r.csv",sep=""), row.names = TRUE)

save.image(paste(path.data,"RealData.Rdata",sep=""))

#
#	End: script_quantmod.R
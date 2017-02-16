#install.packages("MVN")#Quitar el comentario si aun 
#no se ha instalado este paquete
library("MVN") #Paquete para realizar pruebas de 
#distribucion normal multivariada

#Cargamos los datos a usar
datos<-read.csv("datos.txt",header=TRUE,sep=",") 

#Quitamos la columna de indices
datos1<-datos[,2:4] 

#Vemos que hay datos faltantes
summary(datos1) 

#Quitamos los datos faltantes
datos1<-datos1[complete.cases(datos1),] 

#Graficos univariados

# graficos qq univariados
uniPlot(datos1,type="qqplot")

# histigramas univariados
uniPlot(datos1,type="histogram")

#Pruebas Shapiro-Wilk
#Prueba para Apple
shapiro.test(datos1[,1])
#Prueba para Microsoft
shapiro.test(datos1[,2])
#Prueba para Yahoo
shapiro.test(datos1[,3])

#Estadistica descriptiva
uniNorm(datos1,type="SW",desc=TRUE)

#Prueba Mardia
result<-mardiaTest(datos1,qqplot=TRUE)
result

#Prueba Henze-Zirklers
result1<-hzTest(datos1,qqplot=FALSE)
result1


x<-seq(10,25,by=.1)  #crea el rango a graficar, debe ser simétrico%
y<-seq(12,27,by=.1)


s<-matrix(c(1,.3,.3,1.5),ncol=2) # matriz de varianza covarianza%
m<-c(3,5) # vector de medias%
b<-c(2,2) # vector de tendencia%
u<-rpareto(1,5,2) # variable de mezcla%


z<-GN(x,y,m+u*b,u*s) # creando la malla%
plot_ly(x = x, y = y, z = z, type = "surface") # gráfica de la superficie con el paquete plotly%

plot_ly( x=x,y=y,z=z , type = "contour") # gráfica de los contornos con plotly%


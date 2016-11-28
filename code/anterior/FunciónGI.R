# devuelve el el logaritmo natural del valor de la función objetivo 
f<-function(m,a,b,l,n,f){
  x<-m^{a-1}*exp(-f*l*m^{-2}+n*l*m^{-1}-b*m)
  
  if (x >0)
  x<-log(x,base=exp(1))
else
  x<- -Inf
  
x
}
#Variable para actualizar los datos de entrada de f
a<-1
b<-1
c<-1
d<-1
e<-1

#Redefinimos la función para que solamente dependa de x
g<-function(x){
  x<-f(x,a,b,c,d,e)
  x
}

# devuelve el el logaritmo natural del valor de la función objetivo 
h<-function(u,r,m,b,s,l,f){
  x<-exp(-.5(r-m-u*b)*inv(u*s)*(r-m-u*b))*exp(-l*(u-f)^{2}/(2*f^{2}*u))*u^{-1.5-length(r)/2}
  
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
f<-1
#Redefinimos la función para que solamente dependa de x
g<-function(x){
  x<-f(x,a,b,c,d,e,f)
  x
}

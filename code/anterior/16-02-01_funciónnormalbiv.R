# evalua una función normal bivariada dada una matriz nx2 de datos, un vector de medias 
# y una matriz de varianza-covarianza %

Normal<-function(X,M,S){
  
  
  x<-(1/sqrt(2*pi*abs(det(S))))*exp(-(1/2)*(X-M)%*%solve(S)%*%(X-M))
  x
}


# crea una malla para graficar una superficioe normal bivariada%
GN<-function(t,r,M,S){
  j=1
  
  
  n<-length(t)
 
  Z<-matrix(seq(1,n*n),ncol=n)
  
  for(j in 1:n){
    i=1  
    for (i in 1:n) {
                          
           Z[i,j]<-Normal(c(t[i],r[j]),M,S)
                      }          
              }
  Z
}

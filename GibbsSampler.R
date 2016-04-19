#Implementación del muestreador de Gibbs 



#Repositorio para los parámetros de la GI
l<-array(NA,dim=c(M))
f<-array(NA,dim=c(M))
#repositorio de variables latentes
urep<-array(NA,dim=c(M,n))
#Repositorio para los parametros de la normal
m<-array(NA,dim=c(M,p))
b<-array(NA,dim=c(M,p))
S<-array(NA,dim=c(M,p,p))
#-------------------------------------------
#Hiperparámetros
#parámetros correspondietes a alpha
al
bl
#parámetros correspondientes a beta
af
bf
#Parámetros correspondientes a mu
mm<-array(NA,dim=c(p))
SS<-array(NA,dim=c(p,p))
#Parámetros correspondientes a beta
bm<-array(NA,dim=c(p))
bS<-array(NA,dim=c(p,p))
#Parámetros correspondientes a S
Z<-array(NA,dim=c(n,n))
#-------------------------------------------
#Valores iniciales
linicial<-array(NA,dim=c(1))
finicial<-array(NA,dim=c(1))
minicial<-array(NA,dim=c(p))
binicial<-array(NA,dim=c(p))
Sinicial<-array(NA,dim=c(p,p))
u<-array(NA,dim=c(n))
#-------------------------------------------
#Observaciones
r<-array(Na,dim=c(p,n))
#-------------------------------------------
#algoritmo recursivo

  for (j in 1:M){
                    
  #parámetros correspondientes a la distribución normal
                      i=0
                      for(i in 1:n){
                      P<-P+((u[i]-finicial)^2)/2*finicial^2*u[i]
                      X<-X+(r[i,]-minicial-binicial)/sqrt(u[i])
                      V<-V+u[i]
                                 U<-U+1/uinicial[i]
                                 R<-R+r[i,]/uinicial[i]           
                                 N<-N+(r[i,]-minicial)
                      }
                      V<-1/V
                      U<-U^{-1}
                      R<-R*U
                      minicial<-mvrnorm(1,(R*solve(U*Sinicial)+mm*solve(SS))*solve(U*solve(Sinicial)+solve(SS))),U*Sinicial+SS)
                      binicial<-mvrnorm(1,(solve(N*V*Sinicial)+bm*solve(bs))*solve(solve(V*Sinicial)+solve(bs)),V*Sinicial+bs)
                      Sinicial<-rWishart(1,2*n,X*t(X)+Z)
    
      #Parámetros correspondientes a la GI
                      linicial<-rgamma(1,n/2 +bl,P+al )
                      finicial<-uni.slice.alt(x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
                      
                
      #Variables latentes                
                                          i=0
              #Actualiza los parámetros de la función g                            
                                          a<-bf
                                          b<-af
                                          c<-linicial
                                          d<-n
                                          e<-
                                        
                                          
                                          
                                          for(i in 1:n){
                                            u[i]<-uni.slice.alt(x0, h, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
                                          }
                      
  #Llenando los repositorios
  m[j,]<-minicial
  b[j,]<-binicial
  l[j,]<-linicial
  f[j,]<-finicial
  S[j,,]<-Sinicial
  urep[j,]<-u
}














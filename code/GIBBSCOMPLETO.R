#Llamando las funciones adecuadas
#Kernel de la variable latente
setwd("C:/Users/David/Desktop")
getwd()
source("lat.r")
source("latente.r")
#Kernel del parámetro de la distribución Gaussiana Inversa
source("h.r")
source("g.r")
#Slice-Sampler
source("slice.r")

#Librerias
#El paquete MASS/mvrnorm contiene la distribuci???n normal multivariada
library(MASS)
#El paquete stats contiene la distribuci???n Wishart
library(stats)

#Implementaci???n del muestreador de Gibbs 
---------------------------------------------------------------------------
  p<<-2
  M<-50
  n<-50
  j<-0
  #Repositorio para los par???metros de la GI
  lrep<-array(NA,dim=c(M))
  frep<-array(NA,dim=c(M))
  #repositorio de variables latentes
  urep<-array(NA,dim=c(n,M))
  #Repositorio para los parametros de la normal
  mrep<-array(NA,dim=c(p,M))
  brep<-array(NA,dim=c(p,M))
  Srep<-array(NA,dim=c(p,p,M))
  #-------------------------------------------
  #Hiperpar???metros
  #par???metros correspondietes a alpha
  al<-1
  bl<-1
  #par???metros correspondientes a beta
  af<-1
  bf<-1
  #Par???metros correspondientes a mu
  mm<-array(1,dim=c(p))
  SS<-array(c(1,0,0,1),dim=c(p,p))
  #Par???metros correspondientes a beta
  bm<-array(1,dim=c(p))
  bs<-array(c(1,0,0,1),dim=c(p,p))
  #Par???metros correspondientes a S
  Z<-array(c(1,0,0,1),dim=c(p,p))
  #-------------------------------------------
  #Valores iniciales
  linicial<-array(1,dim=c(1))
  finicial<-array(1,dim=c(1))
  minicial<-array(.5,dim=c(p))
  binicial<-array(1,dim=c(p))
  Sinicial<-array(c(1,0,0,1),dim=c(p,p))
  u<-array(.1,dim=c(n))
  #-------------------------------------------
  #Observaciones
  r<-array(1,dim=c(p,n))
  #-------------------------------------------
  #algoritmo recursivo
  j<-1
  for (j in 1:M){
    
    #par???metros correspondientes a la distribuci???n normal
    i<-0
    P<-0
    X<-array(0,dim=c(p))
    V<-0
    U<-0
    R<-array(0,dim=c(p))
    N<-array(0,dim=c(p))
    
    for(i in 1:n){
      P<-P+((u[i]-finicial)^2)/(2*(finicial^2)*u[i])
      X<-X+(r[,i]-minicial-binicial)/sqrt(u[i])
      V<-V+u[i]
      U<-U+1/u[i]
      R<-r[,i]/u[i] +R           
      N<-(r[,i]-minicial)+N
    }
    V<-1/V
    U<-U^{-1}
    H<-array(c(U,0,0,U),dim=c(p,p)) ##Matriz auxiliar para multiplicar una matriz por una constante
    G<-array(c(V,0,0,V),dim=c(p,p))
    R<-R*U
    minicial<-mvrnorm(1,(R%*%solve(H*Sinicial)+mm%*%solve(SS))%*%solve(H*solve(Sinicial)+solve(SS)),H*Sinicial+SS)
    binicial<-mvrnorm(1,(N%*%(solve(G*Sinicial))+bm%*%solve(bs))%*%solve(solve(G*Sinicial)+solve(bs)),V*Sinicial+bs)
    Maux<-rWishart(1,2*n,X%*%t(X)+Z)
    Sinicial<-Maux[,,1]
    
    #Par???metros correspondientes a la GI
    
    
    linicial<-rgamma(1,n/2 +bl,P+al )
    #Actualiza los par???metros de la funci???n g                            
    i<-0
    e<-0
    for(i in 1:n) {
      e<-u[i]+e
      
    }
    a<<-bf
    b<<-af
    cd<<-linicial
    d<<-length(u)
    x0<-1
    
    #h(x0,a,b,cd,d,e,f)
    finicial<-uni.slice.alt(x0, g, w=1, m=Inf, lower=0, upper=+Inf, gx0=NULL)
    
    
    #Variables latentes                
    i<-0
    
    
    #Simulando las variables Latentes                                   
    for(i in 1:n){
      
      rr<<-r[,i]
      mmm<<-minicial
      bb<<-binicial
      ss<<-Sinicial
      ll<<-linicial
      ff<<-finicial
      x0<-1
      #lat(.0001,rr,mmm,bb,ss,ll,ff,p)
      #latente(.0001)
      u[i]<-uni.slice.alt(x0, latente, w=1, m=Inf, lower=0, upper=+Inf, gx0=NULL)
    }
    
    #Llenando los repositorios
    mrep[,j]<-minicial
    brep[,j]<-binicial
    lrep[j]<-linicial
    frep[j]<-finicial
    Srep[,,j]<-Sinicial
    urep[,j]<-u
  }
  #head(lrep,10)
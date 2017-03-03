# Llamando las funciones adecuadas
# Kernel de la variable latente
# setwd("C:/Users/David/Desktop")
rm(list=ls())

setwd("C:/JCMO.Trabajo/@Estudiantes/David Castillo/Tesis/Algoritmo_Final")
getwd()
ptm <- proc.time()
set.seed(1)
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
#library(MCMCpack) No funciona!
library(LaplacesDemon)
#Implementaci???n del muestreador de Gibbs 
#---------------------------------------------------------------------------
#Observaciones
#r<-array(1,dim=c(p,n))
datos<-read.csv("datos.txt",header=TRUE,sep=",") #Cargamos los datos a usar
head(datos)
datos1<-datos[1000:nrow(datos),2:4] #Quitamos los valores que no nos interean
r<-t(datos1) #Transponemos la base de datos
#nrow(datos2)
#ncol(datos2)

#-------------------------------------------
set.seed(76652)
p<<-nrow(r)
#p<<-2   #dimensión de los vectores
M<-10000 #Número de iteraciones
n<<-ncol(r)
#n<-1000 #número de observaciones
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
SS<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p))
#Par???metros correspondientes a beta
bm<-array(1,dim=c(p))
bs<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p))
#Par???metros correspondientes a S
Z<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p))
#-------------------------------------------
#Valores iniciales
linicial<-array(1,dim=c(1))
finicial<-array(1,dim=c(1))
minicial<-array(.5,dim=c(p))
binicial<-array(1,dim=c(p))
Sinicial<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p)) #Lo necesitamos para simular Mu
u<-rep(.1, times=n)
#u<-array(.1,dim=c(n))
#-------------------------------------------
#algoritmo recursivo
for (j in 1:M){
  print(paste("Iteracion: ",j,sep=" - "))
  #par???metros correspondientes a la distribuci???n normal
  i<-0
  P<-0
  X<-array(0,dim=c(p))
  V<-0
  U<-0
  R<-array(0,dim=c(p))
  N<-array(0,dim=c(p))
  
  U<-sum(1/u) #1.1 Parámetro correspondiente a la dist de Mu 
   U<-1/U
  #R<-sum(r/u) #1.2 Parámetro correspondiente a la dist de Mu
   #R<-R*U
  V<-sum(u)   #2.1 Parámetro correspondiente a la dist de B
   V<-1/V
  
   P <- sum((u-finicial)^{2}/(2*u*finicial^{2})) #3.1 Parámetro correspondiente a la distribución de lambda
   
   for(i in 1:n){
#3.1    #P<-P+((u[i]-finicial)^2)/(2*(finicial^2)*u[i])
    X<-X+(r[,i]-minicial-binicial)/sqrt(u[i])
#2.1    V<-V+u[i] Parámetro correspondiente a la dist de B
#1.1    U<-U+1/u[i] #Parámetro correspondiente a la dist de Mu 
    R<-r[,i]/u[i] +R #Parámetro correspondiente a la dist de Mu #1.2           
    N<-(r[,i]-minicial)+N #2.2 Parámetro correspondiente a la distribución de B
  }
   R<-R*U
   #V<-1/V
  #U<-U^{-1}
  H<-array(c(U,0,0,U),dim=c(p,p)) # 1.3 Matriz auxiliar para multiplicar una matriz por una constante
  G<-array(c(V,0,0,V),dim=c(p,p)) # 2.3 Matriz auxiliar para multiplicar una matriz por una constante
    minicial<-mvrnorm(1,(R%*%solve(H*Sinicial)+mm%*%solve(SS))%*%solve(solve(H*Sinicial)+solve(SS)),solve(solve(H*Sinicial)+solve(SS))) # 1.4)Se corrigió está parte del código
    binicial<-mvrnorm(1,(t(G%*%N)%*%(solve(G*Sinicial))+bm%*%solve(bs))%*%solve(solve(G*Sinicial)+solve(bs)),solve(solve(V*Sinicial)+solve(bs)))
  #binicial<-mvrnorm(1,(N%*%(solve(G*Sinicial))+bm%*%solve(bs))%*%solve(solve(G*Sinicial)+solve(bs)),V*Sinicial+bs)
  
  Sinicial<-rinvwishart(2*n,X%*%t(X) + Z)
  #det(X%*%t(X)+Z)
  #Maux<-rWishart(1,2*n,X%*%t(X)+Z)
  #det(Sinicial)
  #Sinicial<-Maux[,,1]
  
  #Par???metros correspondientes a la GI
  
  
  linicial<-rgamma(1,n/2 +bl,P+al)
  #Actualiza los par???metros de la funci???n g                            
  i<-0
  e<-0
  e<-sum(u)
  #for(i in 1:n) {
   # e<-u[i]+e
    
  #}
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
proc.time() - ptm

save.image(file="tesis_davidcastillo_170303_1.Rdata")

#-------------------------------------------
#-------------------------------------------
set.seed(762783)
p<<-nrow(r)
#p<<-2   #dimensión de los vectores
M<-10000 #Número de iteraciones
n<<-ncol(r)
#n<-1000 #número de observaciones
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
SS<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p))
#Par???metros correspondientes a beta
bm<-array(1,dim=c(p))
bs<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p))
#Par???metros correspondientes a S
Z<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p))
#-------------------------------------------
#Valores iniciales
linicial<-array(1,dim=c(1))
finicial<-array(1,dim=c(1))
minicial<-array(.5,dim=c(p))
binicial<-array(1,dim=c(p))
Sinicial<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p)) #Lo necesitamos para simular Mu
u<-rep(.1, times=n)
#u<-array(.1,dim=c(n))
#-------------------------------------------
#algoritmo recursivo
for (j in 1:M){
  print(paste("Iteracion: ",j,sep=" - "))
  #par???metros correspondientes a la distribuci???n normal
  i<-0
  P<-0
  X<-array(0,dim=c(p))
  V<-0
  U<-0
  R<-array(0,dim=c(p))
  N<-array(0,dim=c(p))
  
  U<-sum(1/u) #1.1 Parámetro correspondiente a la dist de Mu 
  U<-1/U
  #R<-sum(r/u) #1.2 Parámetro correspondiente a la dist de Mu
  #R<-R*U
  V<-sum(u)   #2.1 Parámetro correspondiente a la dist de B
  V<-1/V
  
  P <- sum((u-finicial)^{2}/(2*u*finicial^{2})) #3.1 Parámetro correspondiente a la distribución de lambda
  
  for(i in 1:n){
    #3.1    #P<-P+((u[i]-finicial)^2)/(2*(finicial^2)*u[i])
    X<-X+(r[,i]-minicial-binicial)/sqrt(u[i])
    #2.1    V<-V+u[i] Parámetro correspondiente a la dist de B
    #1.1    U<-U+1/u[i] #Parámetro correspondiente a la dist de Mu 
    R<-r[,i]/u[i] +R #Parámetro correspondiente a la dist de Mu #1.2           
    N<-(r[,i]-minicial)+N #2.2 Parámetro correspondiente a la distribución de B
  }
  R<-R*U
  #V<-1/V
  #U<-U^{-1}
  H<-array(c(U,0,0,U),dim=c(p,p)) # 1.3 Matriz auxiliar para multiplicar una matriz por una constante
  G<-array(c(V,0,0,V),dim=c(p,p)) # 2.3 Matriz auxiliar para multiplicar una matriz por una constante
  minicial<-mvrnorm(1,(R%*%solve(H*Sinicial)+mm%*%solve(SS))%*%solve(solve(H*Sinicial)+solve(SS)),solve(solve(H*Sinicial)+solve(SS))) # 1.4)Se corrigió está parte del código
  binicial<-mvrnorm(1,(t(G%*%N)%*%(solve(G*Sinicial))+bm%*%solve(bs))%*%solve(solve(G*Sinicial)+solve(bs)),solve(solve(V*Sinicial)+solve(bs)))
  #binicial<-mvrnorm(1,(N%*%(solve(G*Sinicial))+bm%*%solve(bs))%*%solve(solve(G*Sinicial)+solve(bs)),V*Sinicial+bs)
  
  Sinicial<-rinvwishart(2*n,X%*%t(X) + Z)
  #det(X%*%t(X)+Z)
  #Maux<-rWishart(1,2*n,X%*%t(X)+Z)
  #det(Sinicial)
  #Sinicial<-Maux[,,1]
  
  #Par???metros correspondientes a la GI
  
  
  linicial<-rgamma(1,n/2 +bl,P+al)
  #Actualiza los par???metros de la funci???n g                            
  i<-0
  e<-0
  e<-sum(u)
  #for(i in 1:n) {
  # e<-u[i]+e
  
  #}
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
proc.time() - ptm

save.image(file="tesis_davidcastillo_170303_2.Rdata")

#-------------------------------------------
set.seed(9287)
#p<<-2   #dimensión de los vectores
M<-10000 #Número de iteraciones
n<<-ncol(r)
#n<-1000 #número de observaciones
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
SS<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p))
#Par???metros correspondientes a beta
bm<-array(1,dim=c(p))
bs<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p))
#Par???metros correspondientes a S
Z<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p))
#-------------------------------------------
#Valores iniciales
linicial<-array(1,dim=c(1))
finicial<-array(1,dim=c(1))
minicial<-array(.5,dim=c(p))
binicial<-array(1,dim=c(p))
Sinicial<-array(c(1,0,0,0,1,0,0,0,1),dim=c(p,p)) #Lo necesitamos para simular Mu
u<-rep(.1, times=n)
#u<-array(.1,dim=c(n))
#-------------------------------------------
#algoritmo recursivo
for (j in 1:M){
  print(paste("Iteracion: ",j,sep=" - "))
  #par???metros correspondientes a la distribuci???n normal
  i<-0
  P<-0
  X<-array(0,dim=c(p))
  V<-0
  U<-0
  R<-array(0,dim=c(p))
  N<-array(0,dim=c(p))
  
  U<-sum(1/u) #1.1 Parámetro correspondiente a la dist de Mu 
  U<-1/U
  #R<-sum(r/u) #1.2 Parámetro correspondiente a la dist de Mu
  #R<-R*U
  V<-sum(u)   #2.1 Parámetro correspondiente a la dist de B
  V<-1/V
  
  P <- sum((u-finicial)^{2}/(2*u*finicial^{2})) #3.1 Parámetro correspondiente a la distribución de lambda
  
  for(i in 1:n){
    #3.1    #P<-P+((u[i]-finicial)^2)/(2*(finicial^2)*u[i])
    X<-X+(r[,i]-minicial-binicial)/sqrt(u[i])
    #2.1    V<-V+u[i] Parámetro correspondiente a la dist de B
    #1.1    U<-U+1/u[i] #Parámetro correspondiente a la dist de Mu 
    R<-r[,i]/u[i] +R #Parámetro correspondiente a la dist de Mu #1.2           
    N<-(r[,i]-minicial)+N #2.2 Parámetro correspondiente a la distribución de B
  }
  R<-R*U
  #V<-1/V
  #U<-U^{-1}
  H<-array(c(U,0,0,U),dim=c(p,p)) # 1.3 Matriz auxiliar para multiplicar una matriz por una constante
  G<-array(c(V,0,0,V),dim=c(p,p)) # 2.3 Matriz auxiliar para multiplicar una matriz por una constante
  minicial<-mvrnorm(1,(R%*%solve(H*Sinicial)+mm%*%solve(SS))%*%solve(solve(H*Sinicial)+solve(SS)),solve(solve(H*Sinicial)+solve(SS))) # 1.4)Se corrigió está parte del código
  binicial<-mvrnorm(1,(t(G%*%N)%*%(solve(G*Sinicial))+bm%*%solve(bs))%*%solve(solve(G*Sinicial)+solve(bs)),solve(solve(V*Sinicial)+solve(bs)))
  #binicial<-mvrnorm(1,(N%*%(solve(G*Sinicial))+bm%*%solve(bs))%*%solve(solve(G*Sinicial)+solve(bs)),V*Sinicial+bs)
  
  Sinicial<-rinvwishart(2*n,X%*%t(X) + Z)
  #det(X%*%t(X)+Z)
  #Maux<-rWishart(1,2*n,X%*%t(X)+Z)
  #det(Sinicial)
  #Sinicial<-Maux[,,1]
  
  #Par???metros correspondientes a la GI
  
  
  linicial<-rgamma(1,n/2 +bl,P+al)
  #Actualiza los par???metros de la funci???n g                            
  i<-0
  e<-0
  e<-sum(u)
  #for(i in 1:n) {
  # e<-u[i]+e
  
  #}
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
proc.time() - ptm

save.image(file="tesis_davidcastillo_170303_3.Rdata")

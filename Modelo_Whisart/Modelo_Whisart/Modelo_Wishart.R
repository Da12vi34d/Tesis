getwd()
ptm <- proc.time()
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
library(MVN)
#Implementaci???n del muestreador de Gibbs 
#---------------------------------------------------------------------------
#Observaciones
prueba1<-read.csv("Datos_Modelo_Wishart.csv",stringsAsFactors = FALSE)
r<-t(prueba1[,-1])
#-------------------------------------------
p<<-nrow(r)
M<-15000 #Número de iteraciones
n<<-ncol(r)
matriz<-array(0,dim=c(M,15))
super_rep<-list(matriz,matriz,matriz,matriz,matriz,matriz,matriz,matriz,matriz,matriz)
vsuper_rep<-list(matriz,matriz,matriz,matriz,matriz,matriz,matriz,matriz,matriz,matriz)
#n<-1000 #número de observaciones
for(k in 1:10){
  j<-0
  set.seed(k)
  #Repositorio para los par???metros de la GI
  lrep<-array(NA,dim=c(M))
  frep<-array(NA,dim=c(M))
  #repositorio de variables latentes
  urep<-array(NA,dim=c(n,M))
  #Repositorio para los parametros de la normal
  mrep<-array(NA,dim=c(p,M))
  brep<-array(NA,dim=c(p,M))
  Srep<-array(NA,dim=c(p,p,M)) 
  ###############################################################################################
  ###############################################################################################
  #Repositorios para guardar los resultados de cada iteración
  Super_rep_lrep<-list(lrep,lrep,lrep,lrep,lrep,lrep,lrep,lrep,lrep,lrep)
  Super_rep_frep<-list(frep,frep,frep,frep,frep,frep,frep,frep,frep,frep)
  Super_rep_urep<-list(urep,urep,urep,urep,urep,urep,urep,urep,urep,urep)
  Super_rep_mrep<-list(mrep,mrep,mrep,mrep,mrep,mrep,mrep,mrep,mrep,mrep)
  Super_rep_brep<-list(brep,brep,brep,brep,brep,brep,brep,brep,brep,brep)
  Super_rep_Srep<-list(Srep,Srep,Srep,Srep,Srep,Srep,Srep,Srep,Srep,Srep)
  #-------------------------------------------
  #Hiperpar???metros
  #par???metros correspondietes a alpha
  al<-1
  bl<-0
  #par???metros correspondientes a beta
  af<-1
  bf<-0
  #Par???metros correspondientes a mu
  mm<-array(.00001,dim=c(p))
  SS<-array(c(.00001,0,0,0,.00001,0,0,0,.00001),dim=c(p,p))
  #Par???metros correspondientes a beta
  bm<-array(.00001,dim=c(p))
  bs<-array(c(.00001,0,0,0,.00001,0,0,0,.00001),dim=c(p,p))
  #Par???metros correspondientes a S
  Z<-array(c(.00001,0,0,0,.00001,0,0,0,.00001),dim=c(p,p))
  #-------------------------------------------
  #minicial<-array(.5,dim=c(p))
  binicial<-rep(1,times=p)
  #binicial<-array(1,dim=c(p))
  Sinicial<-array(c(.00001,0,0,0,.00001,0,0,0,.00001),dim=c(p,p)) #Lo necesitamos para simular Mu
  u<-rep(.0001, times=n)
  #u<-array(.1,dim=c(n))
  #Valores iniciales
  linicial<-array(n/sum(1/u - 1/mean(u)),dim=c(1))
  finicial<-array(mean(u),dim=c(1))
  minicial<-rep(1,times=p)
  #-------------------------------------------
  j<-1
  #algoritmo recursivo
  for (j in 1:M){
    
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
    V<-sum(u)   #2.1 Parámetro correspondiente a la dist de B
    V<-1/V
    
    P <- sum((u-finicial)^{2}/(2*u*finicial^{2})) #3.1 Parámetro correspondiente a la distribución de lambda
    
    R<-apply(r/u,1,sum)
    N<-apply(r-minicial,1,sum)
    R<-R*U
    H<-array(c(U,0,0,0,U,0,0,0,U),dim=c(p,p)) # 1.3 Matriz auxiliar para multiplicar una matriz por una constante
    G<-array(c(V,0,0,0,V,0,0,0,V),dim=c(p,p)) # 2.3 Matriz auxiliar para multiplicar una matriz por una constante
    minicial<-mvrnorm(1,(R%*%solve(H*Sinicial)+mm%*%solve(SS))%*%solve(solve(H*Sinicial)+solve(SS)),solve(solve(H*Sinicial)+solve(SS))) # 1.4)Se corrigió está parte del código
    binicial<-mvrnorm(1,(t(G%*%N)%*%(solve(G*Sinicial))+bm%*%solve(bs))%*%solve(solve(G*Sinicial)+solve(bs)),solve(solve(V*Sinicial)+solve(bs)))
    i<-0
    for(i in 1:ncol(r)){
      X<-X+(r[,i]-minicial-u[i]*binicial)/sqrt(u[i]) 
    }
    
    Sinicial<-rinvwishart(2*n,X%*%t(X) + Z)
    linicial<-rgamma(1,n/2 + bl,(P+al)) 
    #Actualiza los par???metros de la funci???n g                            
    e<-0
    e<-sum(u)
    a<<-bf
    b<<-af
    cd<<-linicial
    d<<-length(u)
    x0<-1
    
    finicial<-uni.slice.alt(x0, g, w=1, m=1000, lower=0, upper=+Inf, gx0=NULL)
    
    #Variables latentes                
    i<-1
    
    #Simulando las variables Latentes                                   
    for(i in 1:n){
      
      rr<<-r[,i]
      mmm<<-minicial
      bb<<-binicial
      ss<<-Sinicial
      ll<<-linicial
      ff<<-finicial
      x0<-1
      u[i]<-uni.slice.alt(x0, latente, w=1, m=1000, lower=0, upper=+Inf, gx0=NULL)
    }
    
    #Llenando los repositorios
    mrep[,j]<-minicial
    brep[,j]<-binicial
    lrep[j]<-linicial
    frep[j]<-finicial
    Srep[,,j]<-Sinicial
    urep[,j]<-u
    
    print(c(k,j))
  }
  #################################################################################################
  #Estimadores
  #Estimadores de la matriz de varianza y covarianza
  maz<-array(c(summary(Srep[1,1,2500:M])[4],
               summary(Srep[1,2,2500:M])[4],
               summary(Srep[1,3,2500:M])[4],
               summary(Srep[2,1,2500:M])[4],
               summary(Srep[2,2,2500:M])[4],
               summary(Srep[2,3,2500:M])[4],
               summary(Srep[3,1,2500:M])[4],
               summary(Srep[3,2,2500:M])[4],
               summary(Srep[3,3,2500:M])[4]
  ),dim=c(3,3))
  
  #Estimador de lambda
  lambda<-summary(lrep[2500:M])[4]
  #Estimador de fi
  fi<-summary(frep[2500:M])[4]
  
  #Estimador de beta
  beta<-c(
    summary(brep[1,2500:M])[4],
    summary(brep[2,2500:M])[4],
    summary(brep[3,1000:M])[4])
  
  #Estimador de mu
  mu<-c(summary(mrep[1,2500:M])[4],
        summary(mrep[2,2500:M])[4],
        summary(mrep[3,2500:M])[4])
  #################################################################################################
  #Gráficos convergencia en media (Promedios ergódicos)
  mu1<-cumsum(mrep[1,])/seq_along(mrep[1,])
  mu2<-cumsum(mrep[2,])/seq_along(mrep[1,])
  mu3<-cumsum(mrep[3,])/seq_along(mrep[1,])
  
  beta1<-cumsum(brep[1,])/seq_along(brep[1,])  
  beta2<-cumsum(brep[2,])/seq_along(brep[1,])
  beta3<-cumsum(brep[3,])/seq_along(brep[1,])
  
  lambda_p<-cumsum(lrep)/seq_along(lrep)
  fi_p<-cumsum(frep)/seq_along(frep)
  
  s1<-cumsum(Srep[1,1,])/seq_along(Srep[1,1,])
  s12<-cumsum(Srep[1,2,])/seq_along(Srep[1,1,])
  s13<-cumsum(Srep[1,3,])/seq_along(Srep[1,1,])
  s2<-cumsum(Srep[2,2,])/seq_along(Srep[1,1,])
  s3<-cumsum(Srep[3,3,])/seq_along(Srep[1,1,])
  s23<-cumsum(Srep[2,3,])/seq_along(Srep[1,1,])
  
  #Guardamos las series de los promedios ergódicos 
  parametros<-cbind(iteracion<-c(1:length(mu1)),mu1,mu2,mu3,beta1,beta2,beta3,lambda_p,fi_p,s1,s2,s3,s12,s13,s23)
  ##############################################################################################
  #Gráficos estabilización de la varianza
  vmu1<-cumsum((mrep[1,]-mu[1])^2)/seq_along(mrep[1,])^2
  vmu2<-cumsum((mrep[2,]-mu[2]))^2/seq_along(mrep[1,])^2
  vmu3<-cumsum((mrep[3,]-mu[3])^2)/seq_along(mrep[1,])^2
  
  vbeta1<-cumsum((brep[1,]-beta[1])^2)/seq_along(brep[1,])^2  
  vbeta2<-cumsum((brep[2,]-beta[2])^2)/seq_along(brep[1,])^2
  vbeta3<-cumsum((brep[3,]-beta[3])^2)/seq_along(brep[1,])^2
  
  vlambda_p<-cumsum((lrep-lambda)^2)/seq_along(lrep)^2
  vfi_p<-cumsum((frep-fi)^2)/seq_along(frep)^2
  
  vs1<-cumsum((Srep[1,1,]-maz[1,1])^2)/seq_along(Srep[1,1,])^2
  vs12<-cumsum((Srep[1,2,]-maz[1,2])^2)/seq_along(Srep[1,1,])^2
  vs13<-cumsum((Srep[1,3,]-maz[1,3])^2)/seq_along(Srep[1,1,])^2
  vs2<-cumsum((Srep[2,2,]-maz[2,2])^2)/seq_along(Srep[1,1,])^2
  vs3<-cumsum((Srep[3,3,]-maz[3,3])^2)/seq_along(Srep[1,1,])^2
  vs23<-cumsum((Srep[2,3,]-maz[2,3])^2)/seq_along(Srep[1,1,])^2
  
  #Guardamos las series de estabilización de varianza
  vparametros<-cbind(iteracion<-c(1:length(mu1)),vmu1,vmu2,vmu3,vbeta1,vbeta2,vbeta3,vlambda_p,vfi_p,vs1,vs2,vs3,vs12,vs13,vs23)
  ######Guardamos las series de los promedios ergódicos en super_rep, y las series
  #de estabilización de varianza en vsuper_rep
  
  super_rep[[k]]<-parametros
  vsuper_rep[[k]]<-vparametros
  #Guardamos la serie de los valores simulados
  Super_rep_lrep[[k]]<-lrep
  Super_rep_frep[[k]]<-frep
  Super_rep_urep[[k]]<-Srep
  Super_rep_mrep[[k]]<-mrep
  Super_rep_brep[[k]]<-brep
  Super_rep_Srep[[k]]<-urep
}
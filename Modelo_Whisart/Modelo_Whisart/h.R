h<-function(m,a,b,l,n,f){
  x<-m^{a-1}*exp(-.5*f*l*m^{-2}+n*l*m^{-1}-b*m)

  if (x >0) 
    x<-log(x,base=exp(1))
  else
    x<- -Inf

  x
}

#h<-function(m,a,b,l,n,f){
#  x<-m^{a-1}*exp(-.5*f*l*m^{-2}+n*l*m^{-1}-b*m)
#  x
#}





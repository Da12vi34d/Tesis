lat<-function(x,rr,mmm,bb,ss,ll,ff,p){
  bb<-x*bb
  ss<-solve(x*ss)
  A<-(rr-mmm-bb)%*%ss%*%(rr-mmm-bb)
  x<-exp(-.5*A)*x^{-.5*p-1.5}*exp(-(x-ff)^{2}/(2*x*ff^{2}))
  
  if (x >0) 
    x<-log(x,base=exp(1))
  else
    x<- -Inf
  
  x
}

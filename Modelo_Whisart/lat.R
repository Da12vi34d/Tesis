lat<-function(x,rr,mmm,bb,ss,ll,ff,p){
  bb<-x*bb
  ss<-solve(x*ss)
  A<-(rr-mmm-bb)%*%ss%*%(rr-mmm-bb)
  x<-exp(-.5*A)*x^{-.5*p-1.5}*exp(-ll*(x-ff)^{2}/(2*x*ff^{2}))
  
  if (x >0) 
    x<-log(x,base=exp(1))
  else
    x<- -Inf
  
  x
}


#lat<-function(x0,rr,mmm,bb,ss,ll,ff,p){
#  x0<-100
#  bb<-x0*bb
#  HS<-array(c(x0,0,0,0,x0,0,0,0,x0),dim=c(p,p))
#  ss<-solve(HS%*%ss)
#  A<-(rr-mmm-bb)%*%ss%*%(rr-mmm-bb)
#  ss%*%(rr-mmm-bb)
#  c(crossprod(ss[1,],rr-mmm-bb),
#    crossprod(ss[2,],rr-mmm-bb),
#    crossprod(ss[2,],rr-mmm-bb)
#  )
#  aux<-matrix(rr-mmm-bb,nrow=3,ncol=3)
#  aux_1<-c(crossprod(ss[1,],aux[,1]),crossprod(ss[2,],aux[,2]),crossprod(ss[2,],aux[,3]))
#  aux_2<-c(rr-mmm-bb)
  
#  A<-aux_1[1]*aux_2[1]+aux_1[2]*aux_2[2]+aux_1[3]*aux_2[3]
#  x<-exp(-.5*A)*x0^{-.5*p-1.5}*exp(-ll*(x0-ff)^{2}/(2*x0*ff^{2}))
#  x
}#
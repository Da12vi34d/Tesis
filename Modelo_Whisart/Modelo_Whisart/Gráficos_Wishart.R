#En Total_1 ordenamos de manera conveniente las series de las diez simulaciones de convergencia
#en media.
Total<-rbind(super_rep[[1]],super_rep[[2]],super_rep[[3]],
             super_rep[[4]],super_rep[[5]],super_rep[[6]],super_rep[[7]],super_rep[[8]],
             super_rep[[9]],super_rep[[10]])

Total_1<-as.data.frame(Total)
Total_1$Simulación<-rep(1:10, each = 15000)
Total_1$indice<-rep(1:15000, each = 1, len = 150000)

###############################################################################################
#En VTotal_1 ordenamos de manera conveniente las series de las diez simulaciones de convergencia
#en Varianza.

#Varianza
vTotal<-rbind(vsuper_rep[[1]],vsuper_rep[[2]],vsuper_rep[[3]],vsuper_rep[[4]],vsuper_rep[[5]],
              vsuper_rep[[6]],vsuper_rep[[7]],vsuper_rep[[8]],vsuper_rep[[9]],vsuper_rep[[10]])
vTotal_1<-as.data.frame(vTotal)
vTotal_1$Simulación<-rep(1:10, each = 15000)
vTotal_1$indice<-rep(1:15000, each = 1, len = 150000)
################################################################################################

library("ggplot2")
library("gridExtra")
#install.packages("gridExtra")
#Gráficos de convergencia en media ggplot
p1<-ggplot(Total_1,aes(x=indice,y=mu1))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("m_1")#+scale_y_continuous(limits = c(-.0004, -.0002))+scale_x_continuous(limits = c(13000, 15000))
p2<-ggplot(Total_1,aes(x=indice,y=mu2))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("m_2")#+scale_y_continuous(limits = c(-.0004, -.0002))+scale_x_continuous(limits = c(13000, 15000))
p3<-ggplot(Total_1,aes(x=indice,y=mu3))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("m_3")#+scale_y_continuous(limits = c(-.0004, -.0002))+scale_x_continuous(limits = c(13000, 15000))
grid.arrange(p1, p2, p3, nrow=3)

p4<-ggplot(Total_1,aes(x=indice,y=beta1))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("b_1")+scale_x_continuous(limits = c(14985, 15000))+scale_y_continuous(limits = c( -.0607,-.0605))
p5<-ggplot(Total_1,aes(x=indice,y=beta2))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("b_2")+scale_x_continuous(limits = c(14985, 15000))+scale_y_continuous(limits = c( -.0607,-.0605))
p6<-ggplot(Total_1,aes(x=indice,y=beta3))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("b_3")+scale_x_continuous(limits = c(14985, 15000))+scale_y_continuous(limits = c( -.0607,-.0605))
grid.arrange(p4, p5, p6, nrow=3)

p7<-ggplot(Total_1,aes(x=indice,y=lambda_p))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("lambda")+scale_y_continuous(limits = c(0, 2.5))
p8<-ggplot(Total_1,aes(x=indice,y=fi_p))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("phi")
grid.arrange(p7, p8, nrow=2)


p9<-ggplot(Total_1,aes(x=indice,y=s1))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_1")
p10<-ggplot(Total_1,aes(x=indice,y=s2))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_2")
p11<-ggplot(Total_1,aes(x=indice,y=s3))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_3")
grid.arrange(p9, p10,p11, nrow=3)

p12<-ggplot(Total_1,aes(x=indice,y=s12))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_12")
p13<-ggplot(Total_1,aes(x=indice,y=s13))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_13")
p14<-ggplot(Total_1,aes(x=indice,y=s23))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_23")
grid.arrange(p12,p13,p14,nrow=3)

#################################################################################################
##################################################################################################
#Gráficos de convergencia en varianza
vp1<-ggplot(vTotal_1,aes(x=indice,y=vmu1))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("m_1")
vp2<-ggplot(vTotal_1,aes(x=indice,y=vmu2))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("m_2")
vp3<-ggplot(vTotal_1,aes(x=indice,y=vmu3))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("m_3")
grid.arrange(vp1, vp2, vp3, nrow=3)

vp4<-ggplot(vTotal_1,aes(x=indice,y=vbeta1))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("b_1")
vp5<-ggplot(vTotal_1,aes(x=indice,y=vbeta2))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("b_1")
vp6<-ggplot(vTotal_1,aes(x=indice,y=vbeta3))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("b_1")
grid.arrange(vp4, vp5, vp6, nrow=3)

vp7<-ggplot(vTotal_1,aes(x=indice,y=vlambda_p))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("lambda")
vp8<-ggplot(vTotal_1,aes(x=indice,y=vfi_p))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("phi")
grid.arrange(vp7, vp8, nrow=2)

vp9<-ggplot(vTotal_1,aes(x=indice,y=vs1))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_1")
vp10<-ggplot(vTotal_1,aes(x=indice,y=vs2))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_2")
vp11<-ggplot(vTotal_1,aes(x=indice,y=vs3))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_3")
grid.arrange(vp9, vp10, vp11, nrow=3)

vp12<-ggplot(vTotal_1,aes(x=indice,y=vs12))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_12")
vp13<-ggplot(vTotal_1,aes(x=indice,y=vs13))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_13")
vp14<-ggplot(vTotal_1,aes(x=indice,y=vs23))+geom_line(aes(group=Simulación, colour=Simulación))+xlab(" ")+ylab("")+ggtitle("s_23")
grid.arrange(vp12, vp13, vp14, nrow=3)


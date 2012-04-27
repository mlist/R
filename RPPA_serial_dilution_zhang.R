##############################
#Sips model          ##########
##############################

Sips = function(a=100,M=200, gamma=1,x) {
#a background level
#M saturation level
#gamma non-linear adjustment power
(x^gamma)/((x^gamma)/(M-a)+1)+a
}
 

########################################################################
#fitting the non-linear regression model and making the serial dilution plot
#############################################################################
plot.dilution.series <- function (D0=2, data.dilutes,sensible.min=5, sensible.max=1.e9,minimal.err=5) {
K = ncol(data.dilutes) # number of dilution steps in a dilution series
nsample = nrow(data.dilutes) # number of samples
D0 # this is a preset value in diluteion experiments, typical D0=10, 3, or 2.
#estimate starting parameters
a= max(5, min(data.dilutes, na.rm=T))
M= max(data.dilutes, na.rm =T) 
D=D0

#compile data for model fitting
S.y = as.numeric(data.dilutes[,-K]) # get all columns except the last 
S.x = as.numeric(data.dilutes[,-1])            # get all columns except the first


# making the dilution series plot 
plot(S.x,S.y,pch='+',col='black',xlab='Signal at next dilution step', ylab='Signal')
abline(0,1) # identity line

#testing only setting, sensible value range
filter = ((S.x < sensible.min) | (S.y <sensible.min))
points(S.x[filter],S.y[filter], pch='+',col='red')
filter = ((S.x > sensible.max) | (S.y > sensible.max))
points(S.x[filter],S.y[filter], pch='o',col='red')

S.x1=S.x
S.x1[S.x1<max(0,sensible.min)] =max(0,sensible.min)
S.y1=S.y
S.y1[S.y1<max(0,sensible.min)] =max(0,sensible.min)
ratio= log(abs(S.y1)/abs(S.x1)) #log ratio
ratio.median= median(ratio, na.rm=T)
ratio.mad =mad(ratio,na.rm=T)
filter  =  abs(ratio - ratio.median)/3 > ratio.mad
points(S.x[filter],S.y[filter], pch='+',col='red')

#hist(ratio,nclass=50)


#remove nonsensible values
S.x[filter] = NA
S.x[filter] = NA
S.x[S.x <sensible.min] =NA
S.y[S.y <sensible.min] =NA
S.x[S.x >sensible.max] =NA
S.y[S.y >sensible.max] =NA


#now fitting the mighty nls, model
fit = nls(S.y ~ a +1/((1/(S.x -a) -c)/D+c),start=list(a=a,D=D,c=1/M),alg="port", lower=list(minimal.err,1,0),weights=1/(minimal.err+abs(S.x)))

summary(fit)
a=summary(fit)$parameter[1]
D=summary(fit)$parameter[2]
c=summary(fit)$parameter[3]
d.a= summary(fit)$parameter[4]
d.D= summary(fit)$parameter[5]
d.c= summary(fit)$parameter[6]
M=a+1/summary(fit)$parameter[3]

# making the dilution series plot, adding dilution serial curve 
S.pred = a + 1/((1/(sort(S.x) -a) -1/(M-a))/D+1/(M-a))
points(sort(S.x), S.pred, pch='+', col='blue', type='l',lwd=2)
 

##write parameters on the figure
mtext(paste('a=', format(a,digits=3),' M =', format(M,digits=3),' Dilution factor',format(D,digits=3),'\n\n')) 


#return fitted parameters
list(D=D,c=c, a=a, d.D=d.D, d.c=d.c,d.a=d.a)
} #end of function 

######################################################
#estimate protein concentrations
#####################################################

protein.con <- function (D0,D,c,a,d.D,d.c, d.a, data.dilutes,r=1.2,minimal.err=5) {
#D0 = dilution.factor # this is a preset value in diluteion experiments, typical D0=10, 3, or 2.
#D fitted dilution factor. Ideally, D = D0 ^ gamma, where gamma is a parameter in Sips model
# k 1:ncol(data.dilutes)    # index of dilute dilution steps in each dilution sereies
# Np 1:nrow(data.dilutes) 	# index of samples
x.weighted.mean= rep(NA,nrow(data.dilutes))
x.err = x.weighted.mean
xflag = x.weighted.mean		# takes values of 0,1,2, which means under detection, OK, saturated
K = ncol(data.dilutes) 		# number of total dilution steps for a sample
igamma = log(D0)/log(D)	#where gamma is 1/gamma, a parameter in Sips model
M =min(1e9,1/c+a)			#when M is too large, take 1e9.

x.saturation.level=   D0^(K-1)/((1/( M/r - a)- 1/(M-a)))^igamma 
x.nodetection.level = D0^(1-1)/((1/( r*a - a)- 1/(M-a)))^igamma

for (Np in 1:nrow(data.dilutes)){ # for each sample
   x=rep(NA,K); w=x; xL=x; xH=x;  #initialization
   y = data.dilutes[Np,]
   if((y[K]> M/r) && length(y[y<M/r])<2) {#condition to call saturation
		xflag[Np] = 2 
		x.weighted.mean[Np] = x.saturation.level # Use M/r value
		x.err[Np] = NA
   } else {
   if((y[1]<r*a) & length(y[y>r*a])<2) {#condition to call undetected

		xflag[Np] = 0
		x.weighted.mean[Np] = x.nodetection.level # Use r*a value
		x.err[Np] = NA
   } else {

   y[y>M/r] = NA # for removing signals near saturation 
   y[y<a*r] = NA # for removing signals near bg noise
   
   for (k in 1:K){# for each signal in a dilution series
      y[k] =max(min(M/1.01,y[k]), a+minimal.err) # limit y[k] to be within a+minimal.err and M/1.01
     	x[k] =   D0^(k-1) /(1/(y[k]-a)- c)^igamma #estimated protein concentration prior dilution
	#estimate the derivitives
	de.x.over.de.a = igamma * D0^(k-1)*(1/(y[k]-a)- c)^(-igamma-1)/(y[k]-a)^2

      de.x.over.de.c = igamma * D0^(k-1)*(1/(y[k]-a)- c)^(-igamma-1)
      de.x.over.de.D = x[k] *log(1/(y[k]-a)- c) * igamma/D/log(D)/D0^(k-1)
	w[k] = (de.x.over.de.a * d.a)^2 + ( de.x.over.de.c * d.c)^2 + (de.x.over.de.D * d.D)^2
	}
   w = w[!is.na(x)] # removing signals near saturation or bg noise
   x = x[!is.na(x)] # removing signals near saturation or bg noise

   if(length(x) > 0 ) {
	x.range = 3* max(1, median(x)*0.01,mad(x)) 
	x.f = (abs(x-median(x)) < x.range) # removing outliers
	#c(x,w)
	x=x[x.f]
	w=w[x.f] # removing outliers
	w= 1/w
	x.weighted.mean[Np] = sum (x*w) /sum(w)
	x.err[Np]=1/sqrt(sum(w))
	}
  }#end of else saturation
  }# end of else below detection 

}#end of for each sample Np
#return value:
cbind(x.weighted.mean, x.err,xflag)
}#end of function

##################################################################
#Making plot of Sips Model, dilution series, for the paper:
#Sipps model
x=exp(seq(-5, 5, 0.1))
S =  1 + x/(1+x)
plot(log(x),S, typ='l')
S1=1+x/2/(1+x/2)
plot(S1,S,type='l', lwd=3,xlim=c(.9,2.1),ylim=c(0.9,2.1), xlab='', ylab='')
abline(0,1)
###############################################################

########################################################3##
############
############  examples
####################################################################
r=1.2
#############################################################
#Create simulated data:
#############################################################

D0=2 # preset dilution factor, known from experiments
a=100;M=50000;gamma=1.0
x.t= 2^(seq(1,20,0.02)) # true concentrations
error.rate =0.15		# for introduing mutiplicative error

S = Sips(a=a,M=M,gamma=gamma, x.t) # error free signals
S = S* exp(rnorm(length(S)) *error.rate) #Multiplicate error added
S1 = Sips(a=a,M=M,gamma=gamma, x.t/D0) # signals of sample diluted by a half, error free
S1 = S1 * exp(rnorm(length(S1)) *error.rate) #Multiplicative error added
S2 = Sips(a=a,M=M,gamma=gamma, x.t/D0^2)# signals of sample diluted by 1/4, error free
S2 = S2 * exp(rnorm(length(S2)) *error.rate) #Multiplicative error added
S3 = Sips(a=a,M=M,gamma=gamma, x.t/D0^3)# signals of sample diluted by 1/8, error free
S3 = S3 * exp(rnorm(length(S3)) *error.rate) #Multiplicative error added

data.dilutes= cbind(S,S1,S2,S3)

#op <- par(mfrow = c(3, 2),  pty = "s") 
#plot simulated data
plot(log2(x.t), log2(S), xlab='log2(True Concentration)', ylab='log2(Signal)',pch='+')
#data.dilutes[1,]=c(106.22090,  96.73194,  2000 , 97.44940 ) #make an outlier

param= plot.dilution.series(D0=2, data.dilutes)
a=param$a; c=param$c; D=param$D; 
d.a=param$d.a; d.c=param$d.c; d.D=param$d.D; 

X.result = protein.con(D0=D0,D=D,c=c,a=a,d.a=d.a, d.D=d.D, d.c=d.c,data.dilutes=data.dilutes,r=r) 
x.saturation.level=   D0^(3)/((1/( M/r - a)- 1/(M-a)))^(log(D0)/log(D)) 
x.nodetection.level = D0^(0)/((1/( r*a - a)- 1/(M-a)))^(log(D0)/log(D))

#plot(log2(X.result[,1]),log2(X.result[,2]))
plot( log2(X.result[,1]) , X.result[,2]/X.result[,1],xlab='log2(Estimated Con)', ylab='CV',ylim=c(0,0.2),xlim=c(4,20))
abline(v=log2(x.saturation.level),lty=4,col='green')
abline(v=log2(x.nodetection.level),lty=4,col='blue')


#plot(log2(X.result[,1]), (X.result[,2]/X.result[,1]))
plot(log2(x.t),log2(X.result[,1]),pch='+', xlab='log2(True Concentration)', ylab='log2(Estimated Con)')
abline(0,1,col='red',lwd=3)

abline(log2(x.saturation.level),0,lty=4,col='green')
abline(log2(x.nodetection.level),0,lty=4,col='blue')

#making plot SIgnal vs estimated concentration
#X.result[,1] =x.t
plot(log2(X.result[,1]), log2(S),col='red',pch='+', xlab='log2(Estimated Con)', ylab='log2(Signal)')
points(log2(X.result[,1]/2), log2(S1),col='yellow',pch='+')
points(log2(X.result[,1]/4), log2(S2),col='green',pch='+')
points(log2(X.result[,1]/8), log2(S3),col='blue',pch='+')

############## processing Wei Qing YI data
data.x=read.table("//Bcbwffw7y61/X/Consulting/Qingyi Wei/repair_NER.txt",sep="\t",header=T,fill=T)
k=19
x1=data.x[data.x[,11] ==1,k] 
x2=data.x[data.x[,11] ==2,k] 
x3=data.x[data.x[,11] ==3,k] 
x4=data.x[data.x[,11] ==4,k] 
x5=data.x[data.x[,11] ==5,k] 
x6=data.x[data.x[,11] ==6,k] 

x1[x1<100] =NA
x2[x2<100] =NA
x3[x3<100] =NA
x4[x4<100] =NA
x5[x5<100] =NA
x6[x6<100] =NA

x1[x1>110000] =NA
x2[x2>110000] =NA
x3[x3>11000] =NA
x4[x4>110000] =NA
x5[x5>110000] =NA
x6[x6>110000] =NA

 
plot(x2,x1, xlab='Dilution series k+1', ylab='Dilution series k',xlim=c(0,xM),ylim=c(0,xM))
points(x3,x2, pch='+', col='blue')

points(x5,x4, pch='+', col='blue')
points(x6,x5, pch='+', col='blue')

data.dilutes=cbind(x1+x4,x2+x5,x3+x6)
#data.dilutes=rbind(cbind(x1,x2,x3),cbind(x4,x5,x6))
plot(x2+x5,x1+x4)
D0=2
param= plot.dilution.series(D0 =D0, data.dilutes,sensible.max=110000)
a=param$a; M=param$M; D=param$D; r=1.2
d.a=param$d.a; d.c=param$d.c; d.D=param$d.D; 

X.result = protein.con(D0=D0,D=D,c=c,a=a,r=r,d.a=d.a, d.D=d.D, d.c=d.c,data.dilutes=data.dilutes)
 
plot(log2(X.result[,1]), log2(data.dilutes[,1]),col='red',pch='+', xlab='log2(Estimated Con)', ylab='log2(Signal)',xlim=c(4,11),ylim=c(6,17))
points(log2(X.result[,1]/D0), log2(data.dilutes[,2]),col='green',pch='+')
points(log2(X.result[,1]/D0^2), log2(data.dilutes[,3]),col='blue',pch='+')

#######################cell cycle data ########################################### 
data.x=read.table("//Bcbwffw7y61/X/Consulting/Qingyi Wei/cell cycle array-protein-only.txt",sep="\t",header=T,fill=T)
x1=matrix(data.x[,3], ncol=6,byrow=T)
data.dilutes=rbind(x1[,c(1,3,5)],x1[,c(2,4,6)])

D0=2
param= plot.dilution.series(D0 =D0, sensible.min=1000,sensible.max=63000,data.dilutes)
a=param$a; c=param$c; D=param$D; r=1.2
M=1/c+a
d.a=param$d.a; d.c=param$d.c; d.D=param$d.D; 

X.result = protein.con(D0=D0,D=D,c=c,a=a,r=r,d.a=d.a, d.D=d.D, d.c=d.c,data.dilutes=data.dilutes)
 
plot(log2(X.result[,1]), log2(data.dilutes[,1]),col='red',pch='+', xlab='log2(Estimated Con)', ylab='log2(Signal)',xlim=c(5,13),ylim=c(10,16))
points(log2(X.result[,1]/D0), log2(data.dilutes[,2]),col='green',pch='+')
points(log2(X.result[,1]/D0^2), log2(data.dilutes[,3]),col='blue',pch='+')

#plot error estimate
plot( log2(X.result[,1]) , X.result[,2]/X.result[,1],xlab='log2(Estimated Con)', ylab='CV',ylim=c(0.03,0.1),xlim=c(7.5,12.5))
x.saturation.level =   D0^(3)/((1/( M/r - a)- 1/(M-a)))^(log(D0)/log(D)) 
x.nodetection.level = D0^(0)/((1/( r*a - a)- 1/(M-a)))^(log(D0)/log(D))
abline(v=log2(x.saturation.level),lty=4,col='green')
abline(v=log2(x.nodetection.level),lty=4,col='blue')

plot(matrix(log2(X.result[,1]),ncol=2),xlab='log2(Est.Con1)', ylab='log2(Est.Con2)')
abline(0,1)




 
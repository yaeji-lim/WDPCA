source('~/functions.R', chdir = TRUE)

################################################
############## 			Example 1			 ############## 
###############################################

### data generation
cc =c(.59, .36,  .87, .21)
p=2
n=512

X<-list()

for(model in 1:2){
	if(model==1){
		Y<-matrix(nrow=n, ncol=p)
		for(j in 1:p){
		  if(j <(p/2)){  Y[,j]=  HaarConcat_new(n, c(1,3)) }
		   if(j >=(p/2)){  Y[,j]= HaarConcat_new(n, c(1,5)) } 
		}
	}
		if(model==2){
		Y<-matrix(nrow=n, ncol=p)
		for(j in 1:p){
		  if(j <(p/2)){  Y[,j]=  HaarConcat_new(n, c(1,2,3,4)) }
		   if(j >=(p/2)){  Y[,j]= HaarConcat_new(n, c(2,3)) } 
		}
	}
	
	YY=Y	
	for(i in 1:ncol(Y)){
	  YY[,i]= Y[,i]-apply(Y, 2 , mean)[i]
	}
	
	Y=YY
	X[[model]]=Y
}

### methods

par(mfrow=c(2,2))
for(model in 1:2){
	WPCA_fit=WDPCA(X[[model]])
	}
for(model in 1:2){
	DPCA_fit= FDPCA(X[[model]])
	}

########################################################
############## 	      		Example 2        			 ############## 
#########################################################

### data generation
X=NA
for(sim in 1:10){
  data1=arima.sim(n = 256, list(ar = c(0.8, -0.81) ),sd = 1)
  data2=arima.sim(n = 256, list(ar = c(-0.9, -0.81) ),sd = 1)
  data3=arima.sim(n = 256, list(ar = c(0.8, -0.81) ),sd = 1)
  data11=c( data1[1:63], data2[64:128] ,data3[129:256])
  
  
  data1=arima.sim(n = 256, list(ar = c(0.8, -0.81) ),sd = 1)
  data2=arima.sim(n = 256, list(ar = c(0.6, -0.81) ),sd = 1)
  data3=arima.sim(n = 256, list(ar = c(0.8, -0.81) ),sd = 1)
  data22=c( data1[1:63], data2[64:128] ,data3[129:256])
  
  X=cbind(X, cbind(data11, data22) )
}
X=X[,-1]
dim(X)


### methods
par(mfrow=c(2,2))
WPCA_fit=WDPCA(X)
DPCA_fit= FDPCA(X)


library(IRISSeismic)
library(TSA)
library(wavethresh)
library(mvLSW)
library(eegkit)
library(astsa)
library(MASS)
library(imputeTS)
library(freqdom)
library(condMVNorm)
library(MASS)
library(nowcasting)





 HaarConcat_new<-function(n, order_index=c(1,3)) {
 	
 	y<-vector()
 	for(k in 1:length(order_index)){
    x1 <- HaarMA(n = n/length(order_index), order = order_index[k])
  	 y=c(y,x1)
	}
return(y)
}


close.eigen = function(M,Prev){

  if (!is.matrix(M) || dim(M)[1] != dim(M)[2])

    stop ("M must be a square matrix")

  if (!is.matrix(Prev) || dim(Prev)[1] != dim(Prev)[2])

    stop ("Prev must be a square matrix")

  if (dim(Prev)[1] != dim(M)[2])

    stop ("Dimensions of M and P must be equal")

  

  Eg = eigen(M)

	V = Re(Eg$vectors)

	W = Re(Prev)

	nbasis = dim(M)[1]



	for (col in 1:nbasis){		

		if (sum(V[,col] * W[,col]) < 0)

			Eg$vectors[,col] = -Eg$vectors[,col]

	}



	Eg

}



############## ############## ############## ############## 
############## FDPCA############## ############## 
############## ############## ############## ############## 

FDPCA_global<-function(X, Ndpc = 2){

  res.dpca = dpca(X, Ndpc =2,freq=(-10:10/10) * pi)
  filter_dpca=res.dpca$spec.density$operators  ## f(w) : p *p at each w 
 
 eigenvalue_f<-matrix(nrow= dim(filter_dpca)[3], ncol=Ndpc)
 eigenvector_f<-array(NA,  c(dim(filter_dpca)[3],ncol(X), Ndpc) )
 	for(w in 1:dim(filter_dpca)[3]){
 		eigenvector_f[w,  , ]= prcomp(filter_dpca[,,w])$rotation[, 1: Ndpc]  ## 
 		eigenvalue_f[w,]= (prcomp(filter_dpca[,,w])$sdev[ 1: Ndpc])^2   ## w*q
 	}
  yhat_dpca=dpca.KLexpansion(X, res.dpca$filter)
  
  return(list(freq=res.dpca$spec.density$freq, eigenvector_f =eigenvector_f , eigenvalue_f =eigenvalue_f , yhat_dpca =yhat_dpca))

}



FDPCA<-function(X,q= ncol(X)){

  T=nrow(X)
  p=ncol(X)
  
   N=round(T/35)  ## num of block
  ep=round(T/(N-1))-10
  M=T-(N-1)*ep
  M ## block length
  print(M)
  for(b in 1 ){
    
    print(b)
    t= seq( 1+(b-1)*ep, (b-1)*ep+M)
    block_X= X[t,]
    fit= FDPCA_global(block_X, Ndpc=q)
  } 
 num_f= nrow(fit$eigenvalue_f)
  
  eigenvalues<-array(NA,  c(num_f,  q,N) )
  eigenvectors<-array(NA,  c(num_f,p,q, N) )
  
  for(b in 1:N){
    
    print(b)
    t= seq( 1+(b-1)*ep, (b-1)*ep+M)
    block_X= X[t,]
    fit= FDPCA_global(block_X,Ndpc=q)
    ## frequency band (eg. delta freq, ...)
    eigenvectors[,,,b]=fit$eigenvector_f  ## q th  eigenvector at time u_b , (omega * p * q)
    eigenvalues[,,b] =fit$eigenvalue_f  ##  eigenvalue at time u_b (omega*q)
  } 
  


  for(j in 1:2){
	  image(y=fit$freq, z=t(Re(eigenvalues[,j,])), xlab='time', ylab='frequency', main=paste('Time varying spectrum: ', j,'-th PC', sep=''), col=gray.colors(20), axes=F)
	   axis(1  ,  seq(0, nrow(X), length.out=5)/nrow(X)  , round( seq(0, nrow(X), length.out=5) ) )
    box()
  
  
  }


   return(list( eigenvectors= eigenvectors, eigenvalues= eigenvalues) )
  

}




############## ############## ############## ############## 
############## Wavelet DPCA ############## ############## 
############## ############## ############## ############## 
WDPCA<-function(X){

  T=nrow(X)
  p=ncol(X)
  EWS_smooth <- mvEWS(X, filter.number = 1, family = "DaubExPhase", kernel.name = "modified.daniell", kernel.param = round(sqrt(T)),optimize = FALSE)
  summary(EWS_smooth)
  #plot(EWS_smooth, style = 2, info = 2)                   
  dim(EWS_smooth$spectrum)   ## p times p matrix tilde I(j,omega) for each J(level)=6, k(time)=T            
  
  ## eigen-analysis
  
  wavelet_eigenvalue<-array(NA,   dim(EWS_smooth$spectrum )[c(1,3,4)] )
  wavelet_eigenvectors<-array(NA,   dim(EWS_smooth$spectrum ) )
  for(j in 1:dim(EWS_smooth$spectrum)[3]){ ## for each level
    wav_spectrum=EWS_smooth$spectrum[ , ,j, ]
    E.vectors = array(0, c(p, p, T))
    E.values = array(0, c(p, T))
    Prev = diag(p)
    for (time_id in 1:T) {
      Eg = close.eigen(wav_spectrum[, , time_id], Prev)
      Prev = Eg$vectors
      E.vectors[, , time_id] = Re(Eg$vectors)
      E.values[, time_id] = Re(Eg$values)
    }
    wavelet_eigenvectors[,,j,]=(E.vectors)  # p*q*time, for each scale
    wavelet_eigenvalue[ ,j,]=E.values # q* time, for each scale
    
  }
  


  for(j in 1:2){
    image(t(Re(wavelet_eigenvalue[j,,])), xlab='time', ylab='scale j', main=paste('Time scale dependent spectrum: ', j,'-th PC', sep=''), 
          axes = FALSE,col=gray.colors(10))  
    axis(2  ,seq(1:dim(wavelet_eigenvalue)[2])/(dim(wavelet_eigenvalue)[2]-1)-1/(dim(wavelet_eigenvalue)[2]-1), c(1:dim(wavelet_eigenvalue)[2]))
    axis(1  ,seq(0, dim(wavelet_eigenvalue)[3], length.out=5)/dim(wavelet_eigenvalue)[3], seq(0, dim(wavelet_eigenvalue)[3], length.out=5))
    box()
  }


  return(list(wavelet_eigenvectors=wavelet_eigenvectors,wavelet_eigenvalue=wavelet_eigenvalue ))
  
}




#install.packages("MASS")
#install.packages("sfsmisc")
library(MASS)
library(sfsmisc)

set.seed(123456)

#Q1
N=50
mean=-3
sd=2
LN_amostra=rlnorm(N,meanlog = mean,sdlog = sd)

parametros=c(mean,sd)

#Q2
funcao=function(theta,dados){
  MV=(-length(dados)/2*log(2*pi*theta[2]^2)-sum(log(dados))-sum(log(dados)^2)/(2*theta[2]^2)+sum(log(dados)*theta[1])/theta[2]^2-length(dados)*theta[1]^2/(2*theta[2]^2))
  return(MV)
}

funcao(parametros,LN_amostra)

#Q3
grid_mean=seq(-10,10,1)
grid_sd=seq(1,5,1)

vetor_MV=vector()
vetor_media=vector()
vetor_sd=vector()

for (k in 1:length(grid_sd)){
  for (i in 1:length(grid_mean)){
    vetor_MV[length(grid_mean)*(k-1)+i]=funcao(c(grid_mean[i],grid_sd[k]),LN_amostra)
    vetor_media[length(grid_mean)*(k-1)+i]=grid_mean[i]
    vetor_sd[length(grid_mean)*(k-1)+i]=grid_sd[k]
  }
}
plot(vetor_MV)

#Q4
match(max(vetor_MV),vetor_MV)
vetor_media[match(max(vetor_MV),vetor_MV)]
vetor_sd[match(max(vetor_MV),vetor_MV)]


#Q5
sd0=sqrt((N-1)/N)*sd(log(rlnorm(N,vetor_media[match(max(vetor_MV),vetor_MV)],vetor_sd[match(max(vetor_MV),vetor_MV)])))
sds=c(sd0/sqrt(N),sd0/sqrt(2*N))
VAR_COV=matrix(c(sds[1]^2,0,0,sds[2]^2),ncol=2)
VAR_COV

#teste
fit_LN=fitdistr(LN_amostra,"lognormal")
sd0=sqrt((N-1)/N)*sd(log(LN_amostra))
sds=c(sd0/sqrt(N),sd0/sqrt(2*N))
VAR_COV=matrix(c(sds[1]^2,0,0,sds[2]^2),ncol=2)
VAR_COV
vcov(fit_LN)

#intervalo Normal 95% (supondo que a estimativa para desvio é valor vdd do parametro)
IC_inf_Esp=vetor_media[match(max(vetor_MV),vetor_MV)]-1.96*(vetor_sd[match(max(vetor_MV),vetor_MV)]/sqrt(N))
IC_sup_Esp=vetor_media[match(max(vetor_MV),vetor_MV)]+1.96*(vetor_sd[match(max(vetor_MV),vetor_MV)]/sqrt(N))
c(IC_inf_Esp,IC_sup_Esp)

IC_inf_Sd=(N-1)*vetor_sd[match(max(vetor_MV),vetor_MV)]/71.42
IC_sup_Sd=(N-1)*vetor_sd[match(max(vetor_MV),vetor_MV)]/32.357
c(IC_inf_Sd,IC_sup_Sd)

#teste
IC_inf_Esp=fit_LN[[1]][1]-1.96*(fit_LN[[1]][2]/sqrt(N))
IC_sup_Esp=fit_LN[[1]][1]+1.96*(fit_LN[[1]][2]/sqrt(N))
IC_inf_Sd=fit_LN[[1]][2]*sqrt(N-1)/(sqrt(31.555))
IC_sup_Sd=fit_LN[[1]][2]*sqrt(N-1)/(sqrt(70.222))
confint(fit_LN)

#Q6
MV_est=vector()
Esp_est=vector()
Sd_est=vector()

for (j in 1:100){
  N=50
  mean=-3
  sd=2
  LN_amostra=rlnorm(N,meanlog = mean,sdlog = sd)
  
  parametros=c(mean,sd)

  funcao=function(dados,theta){
    MV=(-length(dados)/2*log(2*pi*theta[2]^2)-sum(log(dados))-sum(log(dados)^2)/(2*theta[2]^2)+sum(log(dados)*theta[1])/theta[2]^2-length(dados)*theta[1]^2/(2*theta[2]^2))
    return(MV)
  }
  
  funcao(parametros,LN_amostra)
  
  grid_mean=seq(-10,10,.1)
  grid_sd=seq(1,10,.1)
  
  vetor_MV=vector()
  vetor_media=vector()
  vetor_sd=vector()
  
  for (k in 1:length(grid_sd)){
    for (i in 1:length(grid_mean)){
      vetor_MV[length(grid_mean)*(k-1)+i]=funcao(LN_amostra,c(grid_mean[i],grid_sd[k]))
      vetor_media[length(grid_mean)*(k-1)+i]=grid_mean[i]
      vetor_sd[length(grid_mean)*(k-1)+i]=grid_sd[k]
    }
  }
  
  MV_est[j]=match(max(vetor_MV),vetor_MV)
  Esp_est[j]=vetor_media[match(max(vetor_MV),vetor_MV)]
  Sd_est[j]=vetor_sd[match(max(vetor_MV),vetor_MV)]
}

hist(Esp_est)
hist(Sd_est)


#Q7
MV_est=vector()
Esp_est=vector()
Sd_est=vector()

for (j in 1:1000){
  N=50
  mean=-3
  sd=2
  LN_amostra=rlnorm(N,meanlog = mean,sdlog = sd)
  
  parametros=c(mean,sd)
  
  funcao=function(dados,theta){
    MV=(-length(dados)/2*log(2*pi*theta[2]^2)-sum(log(dados))-sum(log(dados)^2)/(2*theta[2]^2)+sum(log(dados)*theta[1])/theta[2]^2-length(dados)*theta[1]^2/(2*theta[2]^2))
    return(MV)
  }
  
  funcao(parametros,LN_amostra)
  
  grid_mean=seq(-10,10,.1)
  grid_sd=seq(1,10,.1)
  
  vetor_MV=vector()
  vetor_media=vector()
  vetor_sd=vector()
  
  for (k in 1:length(grid_sd)){
    for (i in 1:length(grid_mean)){
      vetor_MV[length(grid_mean)*(k-1)+i]=funcao(LN_amostra,c(grid_mean[i],grid_sd[k]))
      vetor_media[length(grid_mean)*(k-1)+i]=grid_mean[i]
      vetor_sd[length(grid_mean)*(k-1)+i]=grid_sd[k]
    }
  }
  
  MV_est[j]=match(max(vetor_MV),vetor_MV)
  Esp_est[j]=vetor_media[match(max(vetor_MV),vetor_MV)]
  Sd_est[j]=vetor_sd[match(max(vetor_MV),vetor_MV)]
}

hist(Esp_est)
hist(Sd_est)

#Kolgomorov test
Esp_est_ordered=Esp_est[order(Esp_est,decreasing = F)]
z=(Esp_est_ordered-mean(Esp_est_ordered))/sd(Esp_est_ordered)
F_z=pnorm(z,mean = mean(z),sd=sd(z))
F_s=seq(1/1000,1,1/1000)
abs_Fz_Fs=abs(F_z-F_s)
D=max(abs_Fz_Fs)
D>KSd(length(Esp_est))

#teste
ks.test(Esp_est_ordered, "pnorm", mean(Esp_est), sd(Esp_est))

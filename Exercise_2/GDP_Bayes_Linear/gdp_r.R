df <- read.csv("gdpgrowth.csv")
library(mvtnorm)
y <- df$GR6096
n <- length(y)
X <- matrix(NA,ncol=2,nrow=n)
X[,1] <- rep(1,n)
X[,2] <-df$DEF60

Lambda <- diag(1,n)
K <- diag(0.001,2)

d <- 1
eta <- var(y)
m <- matrix(c(0,0),nrow=2)

d_n <- n + d
K_n <- t(X)%*%Lambda%*%X +K
Sigma_n <- solve(K_n)
m_n <- Sigma_n%*%(t(X)%*%Lambda%*%y + K%*%m)
eta_n <- t(y)%*%Lambda%*%y+ t(m)%*%K%*%m - t(m_n)%*%K_n%*%m_n + eta


################
S <- 10000
omega <-rgamma(S,d_n/2,eta_n/2)
Beta_1 <- matrix(NA,ncol=2,nrow=S)
for(i in 1:S){
  Beta_1[i,] <- rmvnorm(1,m_n,sigma= Sigma_n/omega[i])
}

# plot(X[,2],y=y)
# for(i in 1:S){
#   abline(a=Beta_1[i,1],b=Beta_1[i,2])
# }
############################


##Initial Values
B <- 10000
cur_beta <- c(0,0)
cur_omega <- 1
cur_lambda <- rep(1,n)
num_par <- length(cur_beta)+length(cur_lambda) +1 
par_mat <- matrix(NA,nrow=B,ncol=num_par)



##Prior Param
Lambda <- diag(cur_lambda,n)
K <- diag(0.001,2)
h <- 0.01
d <- 1
eta <- var(y)
m <- matrix(c(0,0),nrow=2)

##Posterior Param
d_n <- n + d
K_n <- t(X)%*%Lambda%*%X +K
eta_n <- t(y)%*%Lambda%*%y+ t(m)%*%K%*%m - t(m_n)%*%K_n%*%m_n + eta
m_n <- Sigma_n%*%(t(X)%*%Lambda%*%y + K%*%m)
Sigma_n <- solve(K_n)



for(i in 1:B){
  ##Update all the parameters for full conditional of beta
  Lambda <- diag(cur_lambda,n)
  K_n <- t(X)%*%Lambda%*%X +K
  Sigma_n <- solve(K_n)
  m_n <- Sigma_n%*%(t(X)%*%Lambda%*%y + K%*%m)

  cur_beta <-rmvnorm(1,mean=m_n,sigma= Sigma_n/cur_omega)
  
  ##Update all the parameters for full conditional of omega
  
  eta_n <- t(y)%*%Lambda%*%y+ t(m)%*%K%*%m - t(m_n)%*%K_n%*%m_n + eta
  cur_omega <- rgamma(1,d_n/2,eta_n/2)
  
  ##Update all the parameters for full conditional of lambda
  mu <- X%*%t(cur_beta)
  cur_lambda <- rgamma(n,rep((h+1)/2,n),(cur_omega*(y-mu)^2+h)/2)
  
  par_mat[i,] <- c(cur_beta,cur_omega,cur_lambda)

}
Beta_2 <- par_mat[,1:2]

final_idx <- seq(from=100,to=B,by=5)
plot(X[,2],y=y,col="darkorange",pch=19,main="Regression plot",xlab="X")

# plot(par_mat[final_idx,1])
abline(a=mean(Beta_2[final_idx,1]),b=mean(Beta_2[final_idx,2]),col="salmon",lwd=3)
abline(a=mean(Beta_1[final_idx,1]),b=mean(Beta_1[final_idx,2]),col="cyan4",lwd=3)
legend(0.12,-0.01, legend=c("Same Lambda","Different Lambda"),
       col=c("cyan4","salmon"), lty=1, cex=0.8)


get_beta_no_thin <- function(X,y){
  B <- 10000
  cur_beta <- c(0,0)
  cur_omega <- 1
  cur_lambda <- rep(1,n)
  num_par <- length(cur_beta)+length(cur_lambda) +1 
  par_mat <- matrix(NA,nrow=B,ncol=num_par)
  
  
  
  ##Prior Param
  Lambda <- diag(cur_lambda,n)
  K <- diag(0.001,2)
  h <- 0.01
  d <- 1
  eta <- var(y)
  m <- matrix(c(0,0),nrow=2)
  
  ##Posterior Param
  d_n <- n + d
  K_n <- t(X)%*%Lambda%*%X +K
  eta_n <- t(y)%*%Lambda%*%y+ t(m)%*%K%*%m - t(m_n)%*%K_n%*%m_n + eta
  m_n <- Sigma_n%*%(t(X)%*%Lambda%*%y + K%*%m)
  Sigma_n <- solve(K_n)
  
  
  
  for(i in 1:B){
    ##Update all the parameters for full conditional of beta
    Lambda <- diag(cur_lambda,n)
    K_n <- t(X)%*%Lambda%*%X +K
    Sigma_n <- solve(K_n)
    m_n <- Sigma_n%*%(t(X)%*%Lambda%*%y + K%*%m)
    
    cur_beta <-rmvnorm(1,mean=m_n,sigma= Sigma_n/cur_omega)
    
    ##Update all the parameters for full conditional of omega
    
    eta_n <- t(y)%*%Lambda%*%y+ t(m)%*%K%*%m - t(m_n)%*%K_n%*%m_n + eta
    cur_omega <- rgamma(1,d_n/2,eta_n/2)
    
    ##Update all the parameters for full conditional of lambda
    mu <- X%*%t(cur_beta)
    cur_lambda <- rgamma(n,rep((h+1)/2,n),(cur_omega*(y-mu)^2+h)/2)
    
    par_mat[i,] <- c(cur_beta,cur_omega,cur_lambda)
    
  }
  Beta_2 <- par_mat[,1:2]
  return(c(mean(Beta_2[,1]),mean(Beta_2[,2])))
}

get_beta_thin <- function(X,y){
  B <- 10000
  cur_beta <- c(0,0)
  cur_omega <- 1
  cur_lambda <- rep(1,n)
  num_par <- length(cur_beta)+length(cur_lambda) +1 
  par_mat <- matrix(NA,nrow=B,ncol=num_par)
  
  
  
  ##Prior Param
  Lambda <- diag(cur_lambda,n)
  K <- diag(0.001,2)
  h <- 0.01
  d <- 1
  eta <- var(y)
  m <- matrix(c(0,0),nrow=2)
  
  ##Posterior Param
  d_n <- n + d
  K_n <- t(X)%*%Lambda%*%X +K
  eta_n <- t(y)%*%Lambda%*%y+ t(m)%*%K%*%m - t(m_n)%*%K_n%*%m_n + eta
  m_n <- Sigma_n%*%(t(X)%*%Lambda%*%y + K%*%m)
  Sigma_n <- solve(K_n)
  
  
  
  for(i in 1:B){
    ##Update all the parameters for full conditional of beta
    Lambda <- diag(cur_lambda,n)
    K_n <- t(X)%*%Lambda%*%X +K
    Sigma_n <- solve(K_n)
    m_n <- Sigma_n%*%(t(X)%*%Lambda%*%y + K%*%m)
    
    cur_beta <-rmvnorm(1,mean=m_n,sigma= Sigma_n/cur_omega)
    
    ##Update all the parameters for full conditional of omega
    
    eta_n <- t(y)%*%Lambda%*%y+ t(m)%*%K%*%m - t(m_n)%*%K_n%*%m_n + eta
    cur_omega <- rgamma(1,d_n/2,eta_n/2)
    
    ##Update all the parameters for full conditional of lambda
    mu <- X%*%t(cur_beta)
    cur_lambda <- rgamma(n,rep((h+1)/2,n),(cur_omega*(y-mu)^2+h)/2)
    
    par_mat[i,] <- c(cur_beta,cur_omega,cur_lambda)
    
  }
  final_idx <- seq(from=100,to=B,by=5)
  Beta_2 <- par_mat[,1:2]
  return(c(mean(Beta_2[final_idx,1]),mean(Beta_2[final_idx,2])))
  
}
num_rep<- 100
no_thin<-t(replicate(num_rep,get_beta_no_thin(X,y)))
thin<-t(replicate(num_rep,get_beta_thin(X,y)))
ans1<- apply(no_thin,2,var)
ans2<- apply(thin,2,var)

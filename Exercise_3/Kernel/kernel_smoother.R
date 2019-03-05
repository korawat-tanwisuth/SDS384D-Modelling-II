
compute_kernel <- function(d,h){
  ##Input: d- distance (scalar), h-bandwith (scalar)
  ##Output: kernel value (scalar)
  return(1/h*dnorm(d/h,mean=0,sd=1))
}

get_weight <- function(x_vec,new_x,ker_func,h){
  ##Input: x_vec- observation (vector), new_x - new x value (scalar), ker_func - function to compute kernel(function)
  dist_x <- x_vec-new_x
  return(ker_func(dist_x,h))
}

n <- 100
x <- runif(n,-2,2)
y <- x^3 + rnorm(n,0,1) 

x_cen <- x - mean(x)
y_cen <- y-mean(y)

h_vec <- c(0.01,0.1,0.2,1,2)
B <- n/2
x_star <- runif(B,-2,2)
get_y_pred<- function(x,x_star,y,h){
  x_cen <- x - mean(x)
  y_cen <- y-mean(y)
  W <- matrix(NA,nrow=length(x),ncol=length(x_star))
  for(i in 1:length(x_star)){
    W[,i]<-get_weight(x_cen,x_star[i],compute_kernel,h)
  }
  
  normalize <- function(x){
    return(x/sum(x))
  }
  W_n <- apply(W,2,normalize)
  y_pred <- rep(NA,length(x_star))
  for(i in 1:length(x_star)){
    y_pred[i] <-sum(W_n[,i]*y_cen)
  }
  return(y_pred)
}
Y_P <- matrix(NA,nrow=B,ncol=length(h_vec))
for(j in 1:length(h_vec)){
  Y_P[,j]<-get_y_pred(h_vec[j])
}
plot(x_cen,y_cen)
col_vec <- 1:length(h_vec)+1
for(k in 1:length(h_vec)){
  x_idx <- order(x_star)
  lines(x_star[x_idx],Y_P[,k][x_idx],col=col_vec[k])
}

#################################

df <- cbind(x,y)

train_idx <- sample(1:n,0.8*n)
train <-  df[train_idx,]
test <- df[-train_idx,]
h <- 0.01

get_y_pred<- function(x_train,x_star,y_train,h){
  ##Input: take vectors x and y from training data and (x_star) from testing data
  ##Return: Predicted values of new y
  W <- matrix(NA,nrow=length(x_train),ncol=length(x_star))
  for(i in 1:length(x_star)){
    W[,i]<-get_weight(x_train,x_star[i],compute_kernel,h)
  }
  
  normalize <- function(x){
    return(x/sum(x))
  }
  W_n <- apply(W,2,normalize)
  y_pred <- rep(NA,length(x_star))
  for(i in 1:length(x_star)){
    y_pred[i] <-sum(W_n[,i]*y_train)
  }
  return(y_pred)
}
get_sse <- function(y_pred,y_star){
  return(sum((y_pred-y_star)^2))
}

get_pred_err <- function(x_train,x_star,y_train,y_star,h){
  
  y_pred <-get_y_pred(x_train,x_star,y_train,h)
  sse <- get_sse(y_pred,y_star)
  err_pred <- cbind(y_star-y_pred,y_pred)
  return(err_pred)
}
err_pred <-get_pred_err(train[,1],test[,1],train[,2],test[,2],0.05)
err <- err_pred[,1]
y_pred <- err_pred[,2]
plot(train[,1],train[,2])
points(test[,1],test[,2],col="red")
points(test[,1],y_pred,col="blue",cex=2)

       
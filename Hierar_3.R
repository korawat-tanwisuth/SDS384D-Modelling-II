library(truncnorm)
library(MASS)
library(MCMCpack)
df <- read.csv("polls.csv")
df <- df[order(df$state),]
df <- df[, -c(1,2,3)]
df <- na.omit(df)
n_vec <- as.numeric(table(df$state))
K <- length(n_vec)
df[, "idx"] <- rep(1:length(n_vec),n_vec)
df[df[,"idx"] == 1,]
X_lst <- vector("list", K) 
group_idx <- df[,"idx"]
for(i in 1:K){
    idx <- group_idx == i
    X_lst[[i]] <- as.matrix(cbind(model.matrix(~edu+age,df[idx,]),df[idx,c("female", "black")]))
}

p <- ncol(X_lst[[1]])
y <- df[,c("bush")]

mu_0 <- rep(0,p)
L_0 <- diag(1, p)
n_0 <- p
s_0 <- diag(1,p)

z <- rep(NA, length(y))

B <- 10000

cur_beta <- matrix(0, K,p)
cur_theta <- rep(0, p)
cur_sigma <- diag(1,p)

beta_arr <- array(dim = c(B,K,p))
theta_mat <- matrix(NA, B, p)
sigma_arr <- array(dim = c(B, p, p))
z_mat <- matrix(NA, B, length(z))
for(iter in 1:B){
    # count <- 0
    # for(i in 1:K){
    #     for(j in 1:n_vec[i]){
    #         count <- count + 1
    #         if(y[count]== 1){
    #             z[count] <- rtruncnorm(1, a = 0, b = Inf, mean = sum(X_lst[[i]][j, ]*cur_beta[i,]), sd = 1)
    #         }
    #         else{
    #             z[count] <- rtruncnorm(1, a = -Inf, b = 0, mean = sum(X_lst[[i]][j, ]*cur_beta[i,]), sd = 1)
    #         }
    #     }   
    # }
    for(i in 1:K){
        idx <- group_idx == i
        y_i <- y[idx]
        z_i <- z[idx]
        idx_1 <- y_i == 1
        idx_0 <- y_i == 0
        mean_z_i <- X_lst[[i]]%*%cur_beta[i, ]
        if(sum(idx_1) > 0){
            z_i[idx_1] <- rtruncnorm(sum(idx_1), a = 0, b = Inf, mean = mean_z_i[idx_1], sd = 1)
        }
        if(sum(idx_0) > 0){
            z_i[idx_0] <- rtruncnorm(sum(idx_0), a = -Inf, b = 0, mean = mean_z_i[idx_0], sd = 1)
        }
        z[idx] <- z_i
    }
    
    ##Update beta
    for(i in 1:K){
        X_i <- X_lst[[i]]
        z_i <- z[group_idx == i]
        temp_sigma <- solve(crossprod(X_i, X_i) + solve(cur_sigma))
        temp_mu <- temp_sigma%*%(t(X_i)%*%z_i + solve(cur_sigma)%*%cur_theta)
        cur_beta[i, ] <- mvrnorm(1, mu = temp_mu, Sigma = temp_sigma)
    }
    
    ##Update theta
    temp_sigma <- solve(K*solve(cur_sigma) + solve(L_0))
    temp_mu <- temp_sigma%*%(solve(cur_sigma)%*%apply(cur_beta,2,sum) + solve(L_0)%*%mu_0)
    cur_theta <- mvrnorm(1, mu = temp_mu, Sigma = temp_sigma)

    ##Update sigma
    s <- matrix(0, p, p)
    for(i in 1:K){
        s <- s + (cur_beta[i,] - cur_theta)%*%t(cur_beta[i,] - cur_theta)
    }
    cur_sigma <- riwish(K + n_0, s + s_0)
    
    ##Store variables
    beta_arr[iter, ,] <- cur_beta
    theta_mat[iter, ] <- cur_theta
    sigma_arr[iter, ,] <- cur_sigma
    z_mat[iter, ] <- z
}
burn_idx <- 0.1*B
z_mat <- z_mat[burn_idx:B, ]
theta_mat <- theta_mat[burn_idx:B, ]
beta_arr <- beta_arr[burn_idx:B,,]
prob_vec <- rep(NA, length(z))
prob_mat <- matrix(NA, nrow = length(burn_idx:B),ncol = length(z))
    
for(iter in 1:length(burn_idx:B)){
    for(i in 1:K){
        temp_beta_i <- beta_arr[iter, i,  ]
        idx <- group_idx == i
        X_i <- X_lst[[i]]
        mean_vec <- (X_i)%*%temp_beta_i
        prob_vec[idx] <- pnorm(mean_vec, sd = 1)
    }
    prob_mat[iter, ] <- prob_vec
}

hist(prob_vec, probability = T)

prob_mean <- apply(prob_mat, 2, mean)
accuracy <- 1- (sum(abs((apply(prob_mat,2,mean)>0.5)-y))/length(y))



pos_theta <- apply(theta_mat,2,mean)
u_theta <- apply(theta_mat, 2, quantile, probs = 0.975)
l_theta <- apply(theta_mat, 2, quantile, probs = 0.025)
significant_var <- !((0<u_theta)&(l_theta<0))
names(significant_var) <- colnames(X_lst[[1]])
##Main effect
significant_var
sig_var_mat <- matrix(NA, nrow = K, ncol = p)
pos_beta_mat <- matrix(NA, nrow = K, ncol = p)

for(i in 1:K){
    u_beta <- apply(beta_arr[,i,], 2, quantile, probs = 0.975)
    l_beta <- apply(beta_arr[,i,], 2, quantile, probs = 0.025)
    temp_var  <- !((0<u_theta)&(l_theta<0))
    sig_var_mat[i, ] <- temp_var  
    pos_beta_mat[i, ] <- apply(beta_arr[,i,], 2, mean)
}


var_names <- colnames(X_lst[[1]])
par(mfrow = c(3,3))
for(i in 1:p){
    hist(theta_mat[,i], probability = T, main = paste("coefficient for", var_names[i]))
}
state_names <- unique(df[,"state"])
par(mar=c(1,1,1,1))
par(mfrow = c(7,7))
prob_by_state <- rep(NA, K)
for(i in 1:K){
    idx <- group_idx ==i
    prob_by_state[i] <- mean(prob_mean[idx])
    hist(prob_mean[idx], main = state_names[i])
}
names(prob_by_state) <- unique(df[,"state"])


library(usmap)
library(ggplot2)

prob_state <- cbind(data.frame(prob_by_state), unique(df[,"state"]))
plot_usmap(
    data = prob_state, values = "prob_by_state", include = unique(df[,"state"]), lines = "red"
) + 
    scale_fill_continuous(
        low = "white", high = "red", name = "Population (2015)", label = scales::comma
    ) + 
    labs(title = "Western US States", subtitle = "These are the states in the Pacific Timezone.") +
    theme(legend.position = "right")




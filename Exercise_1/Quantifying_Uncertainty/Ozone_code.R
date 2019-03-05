# Load the library
# you might have to install this the first time
library(mlbench)

# Load the data
ozone = data(Ozone, package='mlbench')


# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
x = as.matrix(ozone[,2:10])

# add an intercept
x = cbind(1,x)

# compute the estimator
betahat = solve(t(x) %*% x) %*% t(x) %*% y

# Fill in the blank
# betacov = ?

# Now compare to lm
# the 'minus 1' notation says not to fit an intercept (we've already hard-coded it as an extra column)
lm1 = lm(y~x-1)

summary(lm1)
betacovlm = vcov(lm1)
sqrt(diag(betacovlm))

epsilon <- y-x%*%betahat
sigma_hat <- sum(epsilon^2)/(length(y)-length(betahat))
xtx_inv <- solve(t(x) %*% x) 
Cov_beta <- sigma_hat*xtx_inv
std_err <-round(sqrt(diag(Cov_beta)),4)
names(std_err) <- c("Intercept",colnames(x)[-1])

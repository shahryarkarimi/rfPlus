library(rfPlus)
library(scatterplot3d)
set.seed(1)
n <- 200
x1 <- runif(n)
x2 <- runif(n)
y <- sin(2 * pi * x1) + 0.5 * x2 + rnorm(n, sd = 0.2)
d <- data.frame(y, x1, x2)
head(d, 10)
scatterplot3d(
  x = x1, y = x2, z = y,
  xlab = "x1", ylab = "x2", zlab = "y",
  pch = 16, color = "blue"
)
fit_rf <- rf(y ~ ., data = d, ntree = 50)
yhat_rf <- predict(fit_rf, newdata = d)
head(yhat_rf)
# yhat check (lambda = 0)
fit_a0 <- rfPlus(fit_rf, lambda = 0, alpha = 0)
fit_a1 <- rfPlus(fit_rf, lambda = 0, alpha = 1)
max(abs(predict(fit_a0, newdata = d) - yhat_rf))
max(abs(predict(fit_a1, newdata = d) - yhat_rf))
# lambdas
lambdas <- c(0, 0.01, 0.1, 1, 10, 100, 1000)
# ridge
summary_ridge <- data.frame()
fit_ridge <- list()
yhat_ridge <- matrix(0, nr = length(y), nc = length(lambdas))
cf_ridge <- NULL
for (i in 1:length(lambdas)) {
  fit_ridge[[i]] <- rfPlus(fit_rf, lambda = lambdas[i], alpha = 0)
  yhat_ridge[, i] <- predict(fit_ridge[[i]], newdata = d)
  cf_ridge <- cbind(cf_ridge, coef(fit_ridge[[i]], trees = 50)[[1]])
  summary_ridge <- rbind(
    summary_ridge,
    data.frame(
      lambda = lambdas[i],
      intercept = cf_ridge[1, i],
      psi_norm = sqrt(sum(cf_ridge[-1, i]^2)),
      nonzero = sum(cf_ridge[-1, i] != 0),
      pred_sd = sd(yhat_ridge[, i]),
      diff_rf = max(abs(yhat_ridge[, i] - yhat_rf))
    )
  )
}
print(summary_ridge, row.names = FALSE)
# lasso
summary_lasso <- data.frame()
fit_lasso <- list()
yhat_lasso <- matrix(0, nr = length(y), nc = length(lambdas))
cf_lasso <- NULL
for (i in 1:length(lambdas)) {
  fit_lasso[[i]] <- rfPlus(fit_rf, lambda = lambdas[i], alpha = 1)
  yhat_lasso[, i] <- predict(fit_lasso[[i]], newdata = d)
  cf_lasso <- cbind(cf_lasso, coef(fit_lasso[[i]], trees = 50)[[1]])
  summary_lasso <- rbind(
    summary_lasso,
    data.frame(
      lambda = lambdas[i],
      intercept = cf_lasso[1, i],
      psi_norm = sqrt(sum(cf_lasso[-1, i]^2)),
      nonzero = sum(cf_lasso[-1, i] != 0),
      pred_sd = sd(yhat_lasso[, i]),
      diff_rf = max(abs(yhat_lasso[, i] - yhat_rf))
    )
  )
}
print(summary_lasso, row.names = FALSE)
#### predictions vs y, Ridge ####
par(mfrow = c(3, 3))
for (i in 1:length(lambdas)) {
  plot(y, yhat_ridge[, 1], main = paste("Ridge, lambda =", lambdas[i]))
  abline(h = 0, col = 8)
  points(y, yhat_ridge[, i], col = i, pch = 20)
}
#### Plots: predictions vs y, Lasso ####
par(mfrow = c(3, 3))
for (i in 1:length(lambdas)) {
  plot(y, yhat_lasso[, 1], main = paste("Lasso, lambda =", lambdas[i]))
  abline(h = 0, col = 8)
  points(y, yhat_lasso[, i], col = i, pch = 20)
}

#### Plots: coefficient shrinkage, Ridge ####
par(mfrow = c(3, 3))
for (i in 1:length(lambdas)) {
  plot(cf_ridge[, 1], cf_ridge[, 1], main = paste("Ridge, lambda =", lambdas[i]))
  abline(h = 0, col = 8)
  points(cf_ridge[, 1], cf_ridge[, i], col = i, pch = 20)
}
#### Plots: coefficient shrinkage, Lasso ####
par(mfrow = c(3, 3))
for (i in 1:length(lambdas)) {
  plot(cf_lasso[, 1], cf_lasso[, 1], main = paste("Lasso, lambda =", lambdas[i]))
  abline(h = 0, col = 8)
  points(cf_lasso[, 1], cf_lasso[, i], col = i, pch = 20)
}

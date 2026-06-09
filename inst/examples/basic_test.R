#############################################################
# mtcars dataset
#############################################################
library(rfPlus)
set.seed(1)

dat <- mtcars
y_all <- dat$mpg
X_all <- subset(dat, select = -mpg)

### train / test
train_idx <- sample(1:nrow(dat), size = floor(0.7 * nrow(dat)))

train <- dat[train_idx, ]
test  <- dat[-train_idx, ]

y_train <- train$mpg
X_train <- subset(train, select = -mpg)

y_test <- test$mpg
X_test <- subset(test, select = -mpg)


### fitting rf
rf_model <- rf(mpg ~ ., data = train, ntree = 100)

class(rf_model)


### fitting rfPlus

rfp <- rfPlus(rf = rf_model, X = X_train, y = y_train)

cat("Class of rfp:", class(rfp), "\n")

inherits(rfp, "rfPlus")

# print(rfp)

## ------------------------------------------------------------------
## 3. Check Psi for a single tree
## ------------------------------------------------------------------

tp1 <- rfp$tree_info[[1]]

Psi_train_1 <- getPsi(tp1, data = X_train)
Psi_test_1  <- getPsi(tp1, data = X_test)

cat("\nDimensions of Psi (single tree):\n")
cat("  Train:", dim(Psi_train_1))
cat("  Test :", dim(Psi_test_1))

# nrow(Psi_train_1) == nrow(X_train)
# nrow(Psi_test_1) == nrow(X_test)


## ------------------------------------------------------------------
## 4. Check Psi for all trees
## ------------------------------------------------------------------

rft <- rfp$tree_info

Psi_train_full <- getPsi(rft, data = X_train)
Psi_test_full  <- getPsi(rft, data = X_test)

cat("\nDimensions of Psi (all trees):\n")
cat("  Train:", dim(Psi_train_full))
cat("  Test :", dim(Psi_test_full))

# nrow(Psi_train_full) == nrow(X_train)
# nrow(Psi_test_full) == nrow(X_test)

## ------------------------------------------------------------------
## 5. Predictions: RF vs OLS-RF surrogate
## ------------------------------------------------------------------

rf_pred_train  <- predict(rf_model, newdata = X_train)
ols_pred_train <- predict(rfp,    newdata = X_train)

rf_pred_test   <- predict(rf_model, newdata = X_test)
ols_pred_test  <- predict(rfp,    newdata = X_test)

## Compare on test set ----------------------------------------------

comp_test <- data.frame(
  y_true   = y_test,
  rf_pred  = rf_pred_test,
  ols_pred = ols_pred_test,
  diff     = rf_pred_test - ols_pred_test
)

cat("\nFirst 10 rows of comparison on TEST set:\n")
print(round(head(comp_test, 10), 10))

mse_rf_test  <- mean((y_test  - rf_pred_test)^2)
mse_ols_test <- mean((y_test  - ols_pred_test)^2)

cat(sprintf("\nTest MSEs:\n  Random Forest    = %.10f\n", mse_rf_test))
cat(sprintf("  OLS-RF surrogate = %.10f\n", mse_ols_test))

max_abs_diff_test <- max(abs(comp_test$diff))
cat(sprintf("  Max |RF - OLS-RF| on TEST = %.12f\n", max_abs_diff_test))

## Compare on training set ------------------------------------------

comp_train <- data.frame(
  y_true   = y_train,
  rf_pred  = rf_pred_train,
  ols_pred = ols_pred_train,
  diff     = rf_pred_train - ols_pred_train
)

mse_rf_train  <- mean((y_train - rf_pred_train)^2)
mse_ols_train <- mean((y_train - ols_pred_train)^2)

cat(sprintf("\nTrain MSEs:\n  Random Forest    = %.10f\n", mse_rf_train))
cat(sprintf("  OLS-RF surrogate = %.10f\n", mse_ols_train))

max_abs_diff_train <- max(abs(comp_train$diff))
cat(sprintf("  Max |RF - OLS-RF| on TRAIN = %.12f\n", max_abs_diff_train))

## ------------------------------------------------------------------
## 6. Simple pass/fail checks
## ------------------------------------------------------------------

tol <- 1e-8

if (max_abs_diff_train < tol && max_abs_diff_test < tol) {
  cat("\n>>> PASS: rfPlus predictions match randomForest predictions within tolerance.\n")
} else {
  cat("\n>>> WARNING: rfPlus predictions differ from randomForest by more than tolerance.\n")
}




#############################################################
# boston dataset
#############################################################
library(ISLR2)
library(rfPlus)

set.seed(1)

dat <- Boston
y_all <- dat$medv
X_all <- subset(dat, select = -medv)

train_idx <- sample(1:nrow(dat), size = floor(0.7 * nrow(dat)))

train <- dat[train_idx, ]
test <- dat[-train_idx, ]

y_train <- train$medv
X_train <- subset(train, select = -medv)

y_test <- test$medv
X_test <- subset(test, select = -medv)

##### fitting
rf_model <- rf(medv ~ ., data = train, ntree = 100)
class(rf_model)

rfp <- rfPlus(rf = rf_model, X = X_train, y = y_train)
class(rfp)
print(rfp)

#############################################################
rft <- rfp$tree_info  # list of treePlus objects

# extract first treePlus object from rfPlus fit
tp1 <- rfp$tree_info[[1]]
# tp1 <- getTreePlus(rf_model, idx = 1)

# Psi for tree 1 on train and test sets
Psi_train_1 <- getPsi(tp1, data = X_train)
Psi_test_1  <- getPsi(tp1, data = X_test)

dim(Psi_train_1)
dim(Psi_test_1)


Psi_train_full <- getPsi(rft, data = X_train)
Psi_test_full  <- getPsi(rft, data = X_test)

dim(Psi_train_full)
dim(Psi_test_full)

#############################################################
# predictions
#############################################################

rf_pred_test  <- predict(rf_model, newdata = X_test)
ols_pred_test <- predict(rfp,    newdata = X_test)

## compare RF vs rfPlus predictions
comp <- data.frame(
  y_true   = y_test,
  rf_pred  = rf_pred_test,
  ols_pred = ols_pred_test,
  diff     = rf_pred_test - ols_pred_test
)

round(head(comp, 10), 6)

mse_rf  <- mean((y_test - rf_pred_test)^2)
mse_ols <- mean((y_test - ols_pred_test)^2)

cat(sprintf("\nRandom Forest:\n  MSE = %.6f\n", mse_rf))
cat(sprintf("OLS-RF surrogate:\n  MSE = %.6f\n", mse_ols))


#############################################################
# new test
#############################################################

# Exact-equality should still hold at alpha = 0 (untouched path)
set.seed(1)
m <- rf(mpg ~ ., data = mtcars, ntree = 100)
rfp <- rfPlus(m)                                   # X, y defaulted from rf
stopifnot(max(abs(predict(m) - predict(rfp))) < 1e-10)

# Alpha > 0 should NOT collapse predictions toward 0 (the bug we just fixed)
rfp_a10 <- rfPlus(m, alpha = 10)
mean(predict(rfp_a10))    # should be near mean(mtcars$mpg), not near 0
mean(mtcars$mpg)          # ~ 20.09

#' @export
rf <- function(formula = NULL, data = NULL,
               x = NULL, y = NULL, ...) {

  using_formula <- !is.null(formula)
  using_xy <- !is.null(x) && !is.null(y)
  if ((using_formula + using_xy) != 1L) {
    stop("Use either model formula with the data OR x and y (but not both).")
  }

  if (using_formula) {
    if (!inherits(formula, "formula")) stop("Please provide a model formula, e.g. y ~ .")
    if (is.null(data)) stop("A dataset must be provided.")

    if (anyNA(data)) {
      stop("Missing values detected in `data` for variables in the formula. Please clean NAs before calling rf().")
    }

    rf <- randomForest(formula = formula, data = data, keep.inbag = TRUE, ...)
    rf$training_data <- data

  } else {
    X_df <- as.data.frame(x)

    if (anyNA(X_df) || anyNA(y)) {
      stop("Missing values detected in `x` or `y`. Please clean NAs before calling rf().")
    }

    rf <- randomForest(x = X_df, y = y, keep.inbag = TRUE, ...)

    # determine response name
    y_name <- deparse(substitute(y))
    if (!nzchar(y_name)) y_name <- "response"
    y_name <- make.names(y_name, unique = TRUE)

    rf$training_data <- data.frame(X_df, setNames(list(y), y_name), check.names = TRUE)
  }

  rf$call <- match.call()
  class(rf) <- c("rf", class(rf))
  rf
}

perm_data <- function(rf, idx = NULL) {
  if (!inherits(rf, "randomForest")) stop("`rf` must be a randomForest object.")
  if (is.null(rf$training_data)) stop("`rf$training_data` missing. Fit the model with rf().")

  training_data <- as.data.frame(rf$training_data)
  n <- nrow(training_data);
  p <- ncol(training_data)

  fullset <- seq_len(rf$ntree)
  if (is.null(idx)) {
    idx <- fullset
  } else {
    idx <- as.integer(idx)
    if (any(!(idx %in% fullset))) stop("`idx` must be in 1:", rf$ntree)
  }
  t <- length(idx)

  p_data <- array(NA_real_, dim = c(n, p, t),
                  dimnames = list(NULL, colnames(training_data), paste0("tree_", idx)))

  base_idx <- seq_len(n)

  for (j in seq_along(idx)) {
    k <- idx[j]
    w <- rf$inbag[, k]
    reps <- rep.int(base_idx, times = w)
    p_data[, , j] <- as.matrix(training_data[reps, ])
  }

  if (t==1) return(p_data[,,1]) else return(p_data)
}

node_memberships <- function(tree, data) {
  data <- as.data.frame(data)
  n <- nrow(data)
  J <- nrow(tree)

  split_var <- as.character(tree[["split var"]])
  split_thr <- as.numeric(tree[["split point"]])
  left_id   <- as.integer(tree[["left daughter"]])
  right_id  <- as.integer(tree[["right daughter"]])
  status    <- as.integer(tree[["status"]])

  memb <- matrix(FALSE, nrow = n, ncol = J)

  for (i in seq_len(n)) {
    node <- 1L
    repeat {
      memb[i, node] <- TRUE
      if (status[node] == -1L) break
      xval <- data[[ split_var[node] ]][i]
      thr  <- split_thr[node]
      go_left   <- (!is.na(xval)) && (xval <= thr)
      next_node <- if (go_left) left_id[node] else right_id[node]
      if (is.na(next_node) || next_node < 1L || next_node > J) break
      node <- next_node
    }
  }

  memb
}

#' @export
tree2Plus <- function(tree, data) {
  J <- nrow(tree)
  memb <- node_memberships(tree, data)

  treePlus <- tree
  treePlus$NL <- NA_real_
  treePlus$NR <- NA_real_

  status <- as.integer(treePlus[["status"]])
  Ls     <- as.integer(treePlus[["left daughter"]])
  Rs     <- as.integer(treePlus[["right daughter"]])

  for (nid in seq_len(J)) {
    if (status[nid] == -1L) next
    L <- Ls[nid]; R <- Rs[nid]
    if (is.na(L) || is.na(R) || L < 1L || R < 1L || L > J || R > J) next

    NL <- sum(memb[, L])
    NR <- sum(memb[, R])

    if (NL > 0 && NR > 0) {
      treePlus$NL[nid] <- as.numeric(NL)
      treePlus$NR[nid] <- as.numeric(NR)
    }
  }

  treePlus <- list(tree = treePlus, data = data)
  class(treePlus) <- c("treePlus", class(treePlus))

  treePlus
}

#' @export
getTree <- function(x, ...) {
  UseMethod("getTree")
}

#' @export
getTreePlus <- function(x, ...) {
  UseMethod("getTreePlus")
}

#' @export
#' @method getTree rf
getTree.rf <- function (rfobj, k = 1, labelVar = TRUE) {
  if (is.null(rfobj$forest)) {
    stop("No forest component in ", deparse(substitute(rfobj)))
  }
  if (k > rfobj$ntree) {
    stop("There are fewer than ", k, "trees in the forest")
  }
  if (rfobj$type == "regression") {
    tree <- cbind(rfobj$forest$leftDaughter[, k], rfobj$forest$rightDaughter[,
      k], rfobj$forest$bestvar[, k], rfobj$forest$xbestsplit[,
      k], rfobj$forest$nodestatus[, k], rfobj$forest$nodepred[,
      k])[1:rfobj$forest$ndbigtree[k], ]
  }
  else {
    tree <- cbind(rfobj$forest$treemap[, , k], rfobj$forest$bestvar[,
      k], rfobj$forest$xbestsplit[, k], rfobj$forest$nodestatus[,
      k], rfobj$forest$nodepred[, k])[1:rfobj$forest$ndbigtree[k],
      ]
  }
  dimnames(tree) <- list(1:nrow(tree), c("left daughter",
    "right daughter", "split var", "split point", "status",
    "prediction"))
  if (labelVar) {
    tree <- as.data.frame(tree)
    v <- tree[[3]]
    v[v == 0] <- NA
    tree[[3]] <- factor(rownames(rfobj$importance)[v])
    if (rfobj$type == "classification") {
      v <- tree[[6]]
      v[!v %in% 1:nlevels(rfobj$y)] <- NA
      tree[[6]] <- levels(rfobj$y)[v]
    }
  }
  tree
}
#' @export
#' @method getTreePlus rf
getTreePlus.rf <- function(rf, data = NULL, idx = NULL) {

  if (!inherits(rf, "randomForest"))
    stop("getTreePlus(): rf must be a randomForest model created by rf().")

  all_ids <- seq_len(rf$ntree)
  if (is.null(idx)) idx <- all_ids else idx <- as.integer(idx)
  if (any(!(idx %in% all_ids)))
    stop("getTreePlus(): idx must be in 1:", rf$ntree)

  if (is.null(data)) {
    A <- perm_data(rf) # A is either an array or matrix
  } else {
    A <- data
  }

  out <- vector("list", length(idx))

  for (j in seq_along(idx)) {
    k <- idx[j]

    tree_k <- getTree(rf, k, labelVar = TRUE)

    if (is.matrix(A)) {
      Xk <- A
    } else {
      Xk <- A[,,k]
    }

    out[[j]] <- tree2Plus(tree_k, Xk)
  }

  names(out) <- paste0("tree_", idx)

  if (length(out) == 1L) {
    return(out[[1L]])
  } else {
    class(out) <- "treePlus"
    return(out)
  }
}

#' @export
getPsi <- function(x, ...) {
  UseMethod("getPsi")
}

#' @export
#' @method getPsi treePlus
getPsi.treePlus <- function(tp, data = NULL, idx = NULL, ...) {

  # single treePlus object
  if (!is.null(tp$tree) && !is.null(tp$data)) {

    tree_df <- tp$tree


    if (is.null(data)) {
      data <- tp$data
    }
    data <- as.data.frame(data)
    n    <- nrow(data)

    memb <- node_memberships(tree_df, data)

    is_internal <- as.integer(tree_df[["status"]]) != -1L
    has_counts  <- !is.na(tree_df$NL) & !is.na(tree_df$NR)
    sel         <- which(is_internal & has_counts)

    if (length(sel) == 0L) {
      return(matrix(0.0, nrow = n, ncol = 0))
    }

    L  <- as.integer(tree_df[["left daughter"]][sel])
    R  <- as.integer(tree_df[["right daughter"]][sel])
    NL <- as.numeric(tree_df$NL[sel])
    NR <- as.numeric(tree_df$NR[sel])

    S   <- length(sel)
    Psi <- matrix(0.0, nrow = n, ncol = S)

    for (j in seq_len(S)) {
      denom <- sqrt(NL[j] * NR[j])
      if (!is.finite(denom) || denom == 0) next
      Psi[memb[, L[j]], j] <-  NR[j] / denom
      Psi[memb[, R[j]], j] <- -NL[j] / denom
    }

    colnames(Psi) <- paste0("s", seq_len(ncol(Psi)))
    return(Psi)
  }

  # list of treePlus objects
  tp_list <- tp

  if (is.null(data)) {
    stop("When `x` is a list of treePlus objects, you must supply `data`.")
  }
  data <- as.data.frame(data)
  n    <- nrow(data)

  if (!length(tp_list)) {
    return(matrix(0.0, nrow = n, ncol = 0))
  }


  if (!is.null(idx)) {
    idx <- as.integer(idx)
    if (any(idx < 1L | idx > length(tp_list))) {
      stop("`idx` must be between 1 and ", length(tp_list), call. = FALSE)
    }
    tp_list <- tp_list[idx]
  }

  Psi_list <- vector("list", length(tp_list))

  for (j in seq_along(tp_list)) {
    tp_j  <- tp_list[[j]]
    Psi_j <- getPsi(tp_j, data = data, ...)

    if (!is.null(Psi_j) && ncol(Psi_j) > 0L) {
      tree_name <- names(tp_list)[j]
      if (is.null(tree_name) || tree_name == "") tree_name <- paste0("tree_", j)
      colnames(Psi_j) <- paste0(tree_name, "_", colnames(Psi_j))
      Psi_list[[j]]   <- Psi_j
    } else {
      Psi_list[[j]] <- NULL
    }
  }

  Psi_list <- Filter(Negate(is.null), Psi_list)
  if (!length(Psi_list)) {
    return(matrix(0.0, nrow = n, ncol = 0))
  }

  do.call(cbind, Psi_list)
}

#' @export
#' @method getPsi rf
getPsi.rf <- function(rf, data = NULL, idx = NULL, ...) {

  if (is.null(data)) {
    if (is.null(rf$training_data)) {
      stop("`data` is NULL and `rf$training_data` is missing. ",
           "Fit the model with rf() or supply `data`.", call. = FALSE)
    }
    data <- rf$training_data
  }

  tp <- getTreePlus(rf, idx = idx)

  getPsi(tp, data = data, ...)
}

# model
#' @export
#' @export
rfPlus <- function(rf, X = NULL, y = NULL, lambda = 0, alpha = 0) {

  # ---- Internal solvers --------------------------------------------------

  # Closed-form weighted ridge in glmnet's lambda scale.
  # Solves  (X'WX + lambda * N_w * D) beta = X'Wy,
  # where N_w = sum(w) and D is a diagonal mask (1 = penalize, 0 = leave alone).
  # Used when alpha = 0 (any lambda) or when lambda = 0 (any alpha, = OLS).
  .ridge_solve <- function(X, y, w, lambda, penalize) {
    w    <- as.numeric(w)
    N_w  <- sum(w)
    D    <- diag(as.numeric(penalize), nrow = ncol(X))
    Gram <- t(X) %*% (X * w) + lambda * N_w * D
    XtWy <- t(X) %*% (y * w)
    beta <- tryCatch(
      solve(Gram, XtWy),
      error = function(e) {
        stop(
          "Per-tree ridge solve failed -- likely near-singular Psi'WPsi.\n",
          "  lambda = ", lambda, " (consider a larger value).\n",
          "  Underlying error: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
    as.vector(beta)
  }

  # Coordinate descent for weighted elastic net with unpenalized intercept.
  # Objective (matches glmnet's parameterization):
  #   (1 / (2 N_w)) * sum_i w_i (y_i - b0 - x_i' b)^2
  #     + lambda * [ (1 - alpha)/2 * ||b||_2^2  +  alpha * ||b||_1 ]
  #
  # Per-coordinate update (with partial residual r^(j) = y - b0 - X b + X_j b_j):
  #   z_j     = (1/N_w) * sum_i w_i * X_ij * r_i^(j)
  #   v_j     = (1/N_w) * sum_i w_i * X_ij^2
  #   b_j_new = soft(z_j, lambda * alpha) / ( v_j + lambda * (1 - alpha) )
  # where soft(z, g) = sign(z) * max(|z| - g, 0).
  #
  # Intercept update (unpenalized): b0 += (1/N_w) * sum_i w_i r_i.
  # Returns c(intercept, slopes) as a plain numeric vector.
  .enet_solve <- function(X, y, w, lambda, alpha) {

    max_iter <- 1000L
    tol      <- 1e-7

    X   <- as.matrix(X)
    n   <- nrow(X)
    p   <- ncol(X)
    w   <- as.numeric(w)
    N_w <- sum(w)

    # Precompute weighted column scales v_j once.
    v <- colSums(w * X * X) / N_w

    # Soft-thresholding operator.
    soft <- function(z, g) sign(z) * max(abs(z) - g, 0)

    # L1 threshold and L2 shrinkage factor (same across coordinates).
    g_l1 <- lambda * alpha
    g_l2 <- lambda * (1 - alpha)

    # Initialize: intercept = weighted mean of y, all slopes zero.
    b0 <- sum(w * y) / N_w
    b  <- numeric(p)
    r  <- y - b0          # current residuals (b = 0 initially)

    converged <- FALSE
    for (iter in seq_len(max_iter)) {
      max_delta <- 0

      # ---- Coordinate sweep over the slope coefficients ------------------
      for (j in seq_len(p)) {
        if (v[j] == 0) next            # dead column, leave b_j = 0

        # Add current b_j * X[, j] back to residuals -> partial residual.
        r_partial <- r + X[, j] * b[j]

        # Weighted partial covariance with column j.
        z_j <- sum(w * X[, j] * r_partial) / N_w

        # Soft-thresholded, ridge-shrunk update.
        new_bj <- soft(z_j, g_l1) / (v[j] + g_l2)

        delta <- new_bj - b[j]
        if (delta != 0) {
          # Restore residuals using the new b_j.
          r         <- r_partial - X[, j] * new_bj
          b[j]      <- new_bj
          max_delta <- max(max_delta, abs(delta))
        }
      }

      # ---- Intercept update (unpenalized; equiv. weighted mean of r) ----
      shift <- sum(w * r) / N_w
      if (shift != 0) {
        b0        <- b0 + shift
        r         <- r - shift
        max_delta <- max(max_delta, abs(shift))
      }

      if (max_delta < tol) { converged <- TRUE; break }
    }

    if (!converged) {
      warning("Elastic-net coordinate descent did not converge in ",
              max_iter, " iterations (max_delta = ",
              signif(max_delta, 3), ").", call. = FALSE)
    }

    c(b0, b)
  }

  # ---- Input checks ------------------------------------------------------
  if (!inherits(rf, "randomForest"))
    stop("`rf` must be a randomForest object (e.g., created by rf()).")
  if (!is.null(rf$type) && rf$type != "regression")
    stop("This implementation supports regression random forests only.")
  if (is.null(rf$inbag))
    stop("`rf$inbag` is NULL. Fit the random forest with keep.inbag = TRUE (use rf()).")
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0)
    stop("`lambda` must be a single non-negative numeric value.")
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1)
    stop("`alpha` must be a single numeric value in [0, 1] (0 = ridge, 1 = lasso).")

  # ---- Default X / y to the data stashed on the rf object ----------------
  if (is.null(X)) {
    if (is.null(rf$training_data))
      stop("X not supplied and rf$training_data is missing. Pass X explicitly.")
    X <- rf$training_data
    fm <- rf$call$formula
    if (!is.null(fm)) {
      resp <- tryCatch(all.vars(stats::as.formula(fm))[1],
                       error = function(e) NULL)
      if (!is.null(resp) && resp %in% colnames(X))
        X <- X[, setdiff(colnames(X), resp), drop = FALSE]
    }
  }
  X <- as.data.frame(X)

  if (is.null(y)) {
    if (is.null(rf$y))
      stop("y not supplied and rf$y is missing. Pass y explicitly.")
    y <- rf$y
  }
  if (!is.numeric(y))
    stop("`y` must be numeric for regression.")
  if (length(y) != nrow(X))
    stop("length(y) (", length(y), ") must match nrow(X) (", nrow(X), ").")
  if (nrow(rf$inbag) != nrow(X))
    stop("nrow(X) (", nrow(X),
         ") must match the data used to fit `rf` (", nrow(rf$inbag), ").")

  ntree     <- rf$ntree
  coef_list <- vector("list", ntree)

  # ---- treePlus objects for all trees (NULL data => bootstrap NL/NR) -----
  tp_obj  <- getTreePlus(rf)
  tp_list <- if (!is.null(tp_obj$tree) && !is.null(tp_obj$data)) {
    list(tree_1 = tp_obj)
  } else {
    tp_obj
  }

  # ---- Solver dispatch ---------------------------------------------------
  # Closed form covers: lambda = 0 (OLS, any alpha) and alpha = 0 (ridge, any lambda).
  # All other cases (alpha > 0 AND lambda > 0) use hand-coded coordinate descent.
  use_closed_form <- (alpha == 0) || (lambda == 0)

  # ---- Per-tree fit ------------------------------------------------------
  for (k in seq_len(ntree)) {
    tp_k  <- tp_list[[k]]
    Psi_k <- getPsi(tp_k, data = X)
    w     <- rf$inbag[, k]
    inbag <- w > 0

    if (!any(inbag))
      stop("No in-bag observations for tree ", k, " (all zero in rf$inbag).")

    if (ncol(Psi_k) == 0L) {
      # Stump-less tree: weighted mean of y on the in-bag rows.
      mu <- sum(w[inbag] * y[inbag]) / sum(w[inbag])
      coef_list[[k]] <- c(Intercept = mu)
      next
    }

    yd    <- y[inbag]
    wd    <- w[inbag]
    Psi_d <- Psi_k[inbag, , drop = FALSE]

    if (use_closed_form) {
      Xd       <- cbind(Intercept = 1, Psi_d)
      penalize <- c(FALSE, rep(TRUE, ncol(Psi_d)))
      beta     <- .ridge_solve(Xd, yd, wd, lambda, penalize)
      names(beta) <- colnames(Xd)
    } else {
      beta        <- .enet_solve(Psi_d, yd, wd, lambda = lambda, alpha = alpha)
      names(beta) <- c("Intercept", colnames(Psi_d))
    }

    coef_list[[k]] <- beta
  }

  structure(
    list(
      rf            = rf,
      X             = X,
      y             = y,
      lambda        = lambda,
      alpha         = alpha,
      ntree         = ntree,
      feature_names = colnames(X),
      tree_info     = tp_list,
      coef_list     = coef_list,
      call          = match.call()
    ),
    class = "rfPlus"
  )
}

#' @export
print.rfPlus <- function(rfPlus, ...) {
  pen_label <- if (rfPlus$alpha == 0) {
    "ridge"
  } else if (rfPlus$alpha == 1) {
    "lasso"
  } else {
    paste0("elastic net (alpha = ", rfPlus$alpha, ")")
  }
  cat("rfPlus model\n")
  cat("  penalty:        ", pen_label, "\n", sep = "")
  cat("  lambda:         ", rfPlus$lambda, "\n", sep = "")
  cat("  number of trees:", rfPlus$ntree, "\n")
  cat("  features:       ", paste(rfPlus$feature_names, collapse = ", "), "\n")
  invisible(rfPlus)
}

#' @export
coef.rfPlus <- function(rfPlus, trees = NULL, ...) {
  if (is.null(trees)) {
    trees <- seq_len(rfPlus$ntree)
  } else {
    trees <- as.integer(trees)
    if (any(trees < 1L | trees > rfPlus$ntree)) {
      stop("`trees` must be between 1 and ", rfPlus$ntree, call. = FALSE)
    }
  }
  rfPlus$coef_list[trees]
}

#' @export
predict.rfPlus <- function(rfPlus, newdata = NULL, trees = NULL, ...) {

  if (is.null(newdata)) newdata <- rfPlus$X
  newdata <- as.data.frame(newdata)

  if (is.null(trees)) {
    trees <- seq_len(rfPlus$ntree)
  } else {
    trees <- as.integer(trees)
    if (any(trees < 1L | trees > rfPlus$ntree)) {
      stop("`trees` must be between 1 and ", rfPlus$ntree, call. = FALSE)
    }
  }

  nt <- length(trees)
  n  <- nrow(newdata)

  per_tree <- matrix(NA_real_, nrow = n, ncol = nt)

  for (j in seq_along(trees)) {
    k <- trees[j]

    tree_k <- rfPlus$tree_info[[k]]
    Psi_k  <- getPsi(tree_k, data = newdata)
    coef_k <- rfPlus$coef_list[[k]]

    if (ncol(Psi_k) == 0L) {
      per_tree[, j] <- rep(coef_k[1L], n)
    } else {
      intercept <- coef_k[1L]
      beta      <- coef_k[-1L]
      per_tree[, j] <- as.numeric(intercept + Psi_k %*% beta)
    }
  }

  rowMeans(per_tree)
}

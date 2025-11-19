#' @import randomForest

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
    if (status[nid] == -1L) next  # skip leaves
    L <- Ls[nid]; R <- Rs[nid]
    if (is.na(L) || is.na(R) || L < 1L || R < 1L || L > J || R > J) next

    # counts over bootstrapped rows (duplicates explicit)
    NL <- sum(memb[, L])
    NR <- sum(memb[, R])

    if (NL > 0 && NR > 0) {
      treePlus$NL[nid] <- as.numeric(NL)
      treePlus$NR[nid] <- as.numeric(NR)
    }
  }
  class(treePlus) <- c("treePlus", class(tree))
  treePlus <- list(tree=treePlus, data=data)
}

getTreePlus <- function(rf, data = NULL, idx = NULL) {

  if (!inherits(rf, "randomForest"))
    stop("getRfPlus(): rf must be a randomForest model created by rf().")

  all_ids <- seq_len(rf$ntree)
  if (is.null(idx)) idx <- all_ids else idx <- as.integer(idx)
  if (any(!(idx %in% all_ids)))
    stop("getRfPlus(): idx must be in 1:", rf$ntree)

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

build_Psi <- function(x, ...) {
  UseMethod("build_Psi")
}

build_Psi.treePlus <- function(treePlus, data) {
  data <- as.data.frame(data)
  n <- nrow(data)

  memb <- node_memberships(treePlus, data)


  is_internal <- as.integer(treePlus[["status"]]) != -1L
  has_counts  <- !is.na(treePlus$NL) & !is.na(treePlus$NR)
  idx <- which(is_internal & has_counts)

  if (length(idx) == 0L) {
    return(matrix(0.0, nrow = n, ncol = 0))
  }

  # Children and counts for the selected splits
  L  <- as.integer(treePlus[["left daughter"]][idx])
  R  <- as.integer(treePlus[["right daughter"]][idx])
  NL <- as.numeric(treePlus$NL[idx])
  NR <- as.numeric(treePlus$NR[idx])

  S <- length(idx)
  Psi <- matrix(0.0, nrow = n, ncol = S)

  for (j in seq_len(S)) {
    denom <- sqrt(NL[j] * NR[j])
    if (!is.finite(denom) || denom == 0) next
    Psi[memb[, L[j]], j] <-  NR[j] / denom   # + on left child
    Psi[memb[, R[j]], j] <- -NL[j] / denom   # - on right child
  }

  colnames(Psi) <- paste0("s", seq_len(ncol(Psi)))
  Psi
}

# Task for Shahryar
build_Psi.RftPlus <- function(RftPlus_obj, data, idx = NULL, ...) {

  data <- as.data.frame(data)
  n <- nrow(data)

  if (!length(RftPlus_obj)) {
    return(matrix(0.0, nrow = n, ncol = 0))
  }

  if (!is.null(idx)) {
    idx <- as.integer(idx)
    if (any(idx < 1L | idx > length(RftPlus_obj))) {
      stop("`idx` must be between 1 and ", length(RftPlus_obj), call. = FALSE)
    }
    RftPlus_obj <- RftPlus_obj[idx]
  }

  Psi_list <- vector("list", length(RftPlus_obj))

  for (j in seq_along(RftPlus_obj)) {
    tree_j <- RftPlus_obj[[j]]

    Psi_j  <- build_Psi(tree_j, data, ...)

    if (!is.null(Psi_j) && ncol(Psi_j) > 0L) {
      tree_name <- names(RftPlus_obj)[j]
      if (is.null(tree_name) || tree_name == "") tree_name <- paste0("tree_", j)
      colnames(Psi_j) <- paste0(tree_name, "_", colnames(Psi_j))
      Psi_list[[j]] <- Psi_j
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

# model
rfPlus <- function(rf, X, y) {
  # weighted normal equations helper
  .normal_eqs <- function(X, y, w) {
    w <- as.numeric(w)
    XtW <- t(X) %*% (X * w)
    XtWy <- t(X) %*% (y * w)
    solve(XtW, XtWy)
  }

  # --- basic checks -------------------------------------------------------
  if (!inherits(rf, "randomForest"))
    stop("`rf` must be a randomForest object (e.g., created by rf()).")

  if (!is.null(rf$type) && rf$type != "regression")
    stop("This implementation supports regression random forests only.")

  if (is.null(rf$inbag))
    stop("`rf$inbag` is NULL. Fit the random forest with keep.inbag = TRUE (use rf()).")

  X <- as.data.frame(X)
  if (!is.numeric(y))
    stop("`y` must be numeric for regression.")

  if (length(y) != nrow(X))
    stop("Length of y must match nrow(X).")

  if (nrow(rf$inbag) != nrow(X))
    stop("Number of rows in X must match number of rows used to fit `rf`.")

  ntree <- rf$ntree
  coef_list <- vector("list", ntree)

  RfPlus_trees <- getRfPlus(rf)

  if (!is.list(RfPlus_trees) || inherits(RfPlus_trees, "treePlus")) {
    RfPlus_trees <- list(RfPlus_trees)
    names(RfPlus_trees) <- paste0("tree_", seq_len(ntree))
  }

  for (k in seq_len(ntree)) {
    tree_k <- RfPlus_trees[[k]]

    Psi <- build_Psi(tree_k, X)

    w <- rf$inbag[, k]
    inbag <- w > 0

    if (ncol(Psi) == 0L) {
      if (!any(inbag))
        stop("No in-bag observations for tree ", k, " (all zero in `rf$inbag`).")
      mu <- sum(w[inbag] * y[inbag]) / sum(w[inbag])
      coef_list[[k]] <- c(mu)
    } else {
      Xd <- cbind(Intercept = 1, Psi)[inbag, , drop = FALSE]
      yd <- y[inbag]
      wd <- w[inbag]
      coef_list[[k]] <- as.vector(.normal_eqs(Xd, yd, wd))
    }
  }

  structure(
    list(
      rf            = rf,
      X             = X,
      y             = y,
      ntree         = ntree,
      feature_names = colnames(X),
      tree_info     = RfPlus_trees,
      coef_list     = coef_list,
      call          = match.call()
    ),
    class = "rfPlus"
  )
}

print.rfPlus <- function(rfPlus, ...) {
  cat("RfPlus model\n")
  cat("number of trees:", rfPlus$ntree, "\n")
  cat("features:", paste(rfPlus$feature_names, collapse = ", "), "\n")
  invisible(rfPlus)
}

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

  # matrix of per-tree predictions
  per_tree <- matrix(NA_real_, nrow = n, ncol = nt)

  for (j in seq_along(trees)) {
    k <- trees[j]

    tree_k <- rfPlus$tree_info[[k]]
    Psi_k  <- build_Psi(tree_k, newdata)

    coef_k <- rfPlus$coef_list[[k]]

    if (ncol(Psi_k) == 0L) {
      # intercept-only tree
      per_tree[, j] <- rep(coef_k[1L], n)
    } else {
      intercept <- coef_k[1L]
      beta      <- coef_k[-1L]
      per_tree[, j] <- as.numeric(intercept + Psi_k %*% beta)
    }
  }

  rowMeans(per_tree)
}



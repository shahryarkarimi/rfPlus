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
get.Tree <- function(x, ...) {
  UseMethod("get.Tree")
}

#' @export
get.TreePlus <- function(x, ...) {
  UseMethod("get.TreePlus")
}

#' @export
#' @method get.Tree rf
get.Tree.rf <- function (rfobj, k = 1, labelVar = TRUE) {
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
#' @method get.TreePlus rf
get.TreePlus.rf <- function(rf, data = NULL, idx = NULL) {

  if (!inherits(rf, "randomForest"))
    stop("get.TreePlus(): rf must be a randomForest model created by rf().")

  all_ids <- seq_len(rf$ntree)
  if (is.null(idx)) idx <- all_ids else idx <- as.integer(idx)
  if (any(!(idx %in% all_ids)))
    stop("get.TreePlus(): idx must be in 1:", rf$ntree)

  if (is.null(data)) {
    A <- perm_data(rf) # A is either an array or matrix
  } else {
    A <- data
  }

  out <- vector("list", length(idx))

  for (j in seq_along(idx)) {
    k <- idx[j]

    tree_k <- get.Tree(rf, k, labelVar = TRUE)

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

  tp <- get.TreePlus(rf, idx = idx)

  getPsi(tp, data = data, ...)
}

# model
#' @export
rfPlus <- function(rf, X, y) {
  # weighted normal equations helper
  .normal_eqs <- function(X, y, w) {
    w <- as.numeric(w)
    XtW  <- t(X) %*% (X * w)
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

  # --- get treePlus objects for all trees ---------------------------------
  tp_obj <- get.TreePlus(rf)  # all trees by default

  # Ensure we always work with a list of treePlus objects
  if (!is.null(tp_obj$tree) && !is.null(tp_obj$data)) {
    # single treePlus object
    tp_list <- list(tp_obj)
    names(tp_list) <- "tree_1"
  } else {
    # already a list of treePlus objects
    tp_list <- tp_obj
  }

  # --- per-tree OLS -------------------------------------------------------
  for (k in seq_len(ntree)) {
    tp_k <- tp_list[[k]]

    # Psi for tree k on the training X
    Psi_k <- getPsi(tp_k, data = X)


    w <- rf$inbag[, k]
    inbag <- w > 0

    if (ncol(Psi_k) == 0L) {

      if (!any(inbag))
        stop("No in-bag observations for tree ", k, " (all zero in `rf$inbag`).")

      mu <- sum(w[inbag] * y[inbag]) / sum(w[inbag])
      coef_list[[k]] <- c(mu)
    } else {
      Xd <- cbind(Intercept = 1, Psi_k)[inbag, , drop = FALSE]
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
      tree_info     = tp_list,
      coef_list     = coef_list,
      call          = match.call()
    ),
    class = "rfPlus"
  )
}

#' @export
print.rfPlus <- function(rfPlus, ...) {
  cat("RfPlus model\n")
  cat("number of trees:", rfPlus$ntree, "\n")
  cat("features:", paste(rfPlus$feature_names, collapse = ", "), "\n")
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

    # updated: use getPsi(), not build_Psi()
    Psi_k <- getPsi(tree_k, data = newdata)

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

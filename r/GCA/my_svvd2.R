hard_thresh_rows <- function(x, gamma){
  I.rows = c()
  for(i in 1:dim(x)[1]){
    if(sum(x[i, ]^2) > gamma^2){
      I.rows = c(I.rows, i)
    }
    else{
      x[i, ] = rep(0, length(x[i,]))
    }
  }
  return(list(I.rows = I.rows, x = x))
  
}


my.ssvd2 <-
  function (x, Sigma_u, Sigma_v,
            method = c("theory", "method"), alpha.method = 0.05,
            alpha.theory = 1.5, huber.beta = 0.95, sigma = NA, r = 1,
            gamma.u = sqrt(8), gamma.v = sqrt(8), dothres = "hard", tol = 1e-08,
            n.iter = 100, n.boot = 100, non.orth = FALSE, reps = 1)
    # the main function
  {
    ans.initial <- ssvd.initial2(x, Sigma_u, Sigma_v, method = method, alpha.method = alpha.method,
                                alpha.theory = alpha.theory, huber.beta = huber.beta,
                                sigma = sigma, r = r)
    my.ssvd.iter.thresh2(x, Sigma_u, Sigma_v,
                        method = method, u.old = ans.initial$u,
                        v.old = ans.initial$v, gamma.u = gamma.u, gamma.v = gamma.v,
                        dothres = dothres, r = r, tol = tol, n.iter = n.iter,
                        n.boot = n.boot, sigma = ans.initial$sigma.hat, 
                        non.orth = non.orth, reps = reps)
  }

ssvd.initial2 <-
  function (x, Sigma_u, Sigma_v,
            method = c("theory", "method"), alpha.method = 0.05,
            alpha.theory = 1.5, huber.beta = 0.95, sigma = NA, r = 1)
    # implement SSVD initial selection in both the methodology and theoretical paper
    # to use theoretical one, set method to "theory", and set alpha.theory
    # to use methodology, set method to the "method", and set alpha.method, huber.beta
    # sigma can be provided if known, otherwise will be estimated
    # r is the desired rank
  {
    x <- as.matrix(x)
    pu <- nrow(x)
    pv <- ncol(x)
    # estimate sigma
    if (is.na(sigma)) {
      sigma.hat <- mad(as.vector(x))
    }
    else {
      sigma.hat <- sigma
    }
    # huberize x^2
    if (method == "theory") {
      x.huber <- x^2
    }
    else {
      x.huber <- huberize(x, huber.beta = huber.beta)
    }
    # row and column selection
    rownorm2 <- apply(x.huber, 1, sum)
    colnorm2 <- apply(x.huber, 2, sum)
    I.row <- get.subset(rownorm2, method = method, alpha.method = alpha.method,
                        alpha.theory = alpha.theory, sigma = sigma.hat, df = pv)
    print(I.row)
    I.col <- get.subset(colnorm2, method = method, alpha.method = alpha.method,
                        alpha.theory = alpha.theory, sigma = sigma.hat, df = pu)
    # sanitary
    if (length(I.row) < r) {
      warning("SSVD.initial: Number of selected rows less than rank!")
      I.row <- (order(rownorm2, decreasing = TRUE))[1:min(r +
                                                            10, pu)]
    }
    if (length(I.col) < r) {
      warning("SSVD.initial: Number of selected cols less than rank!")
      I.col <- (order(colnorm2, decreasing = TRUE))[1:min(r +
                                                            10, pv)]
    }
    # SVD on selected submatrix
    sigma_u.sub = Sigma_u[I.row, I.row, drop = FALSE]
    sigma_v.sub = Sigma_v[I.col, I.col, drop = FALSE]
    x.sub = x[I.row, I.col, drop = FALSE] 
    cor.sub = sqrtm(sigma_u.sub)$Binv %*% x.sub %*% sqrtm(sigma_v.sub)$Binv 
    x.sub.svd <- svd(cor.sub, nu = r, nv = r)
    # expanding
    u.hat <- matrix(0, pu, r)
    v.hat <- matrix(0, pv, r)
    u.hat[I.row, ] <- sqrtm(sigma_u.sub)$Binv %*% x.sub.svd$u
    v.hat[I.col, ] <- sqrtm(sigma_v.sub)$Binv %*% x.sub.svd$v
    d.hat <- x.sub.svd$d[1:r]
    list(u = u.hat, v = v.hat, d = d.hat, sigma.hat = sigma.hat)
  }

my.ssvd.iter.thresh2 <-
  function (x, Sigma_u, Sigma_v,
            method = c("theory", "method"), u.old, v.old, gamma.u = sqrt(2),
            gamma.v = sqrt(2), dothres = "hard", r = ncol(u.old), tol = 1e-08,
            reps = 1,
            n.iter = 100, n.boot = 100, sigma = NA, non.orth = FALSE)
    # x is the observed matrix
    # u.old and v.old are the starting points
    # gamma.u and gamma.v are for theoretical computation:
    # threshold level is set to be gamma*sqrt(log(p))
    # n.iter: max number of iteration
    # n.boot: number of bootstrap
    # sigma: noise level
    # r: rank
    # tol: tolerance level for the convergence
    # non.orth: if set to be TRUE, then the last step does not involve orthoganalization
  {
    x <- as.matrix(x)
    pu <- nrow(x)
    pv <- ncol(x)
    puv <- max(pu, pv)
    # estimate sigma
    if (is.na(sigma)) {
      sigma.hat <- mad(as.vector(x))
    }
    else {
      sigma.hat <- sigma
    }
    # rescaling
    x.scaled <- x/sigma.hat
    # thresholding rule
    if (dothres == "hard") {
      thresh <- hard_thresh_rows 
    }
    else if (dothres == "soft") {
      thresh <- soft.thresh
    }
    else {
      warning("SSVD.iter.thresh: argument dothres not recognized! Use hard-thresholding as default")
      thresh <- hard_thresh_rows 
    }
    # initialization
    dist.u <- 1
    dist.v <- 1
    i.iter = 1
    u.cur <- u.old
    v.cur <- v.old
    
    
    
    while (i.iter <= n.iter & max(dist.u, dist.v) > tol) {
      u.old <- u.cur
      # multiplication
      sel.v <- apply(v.cur == 0, 1, all)
      
      u.cur <- x.scaled[,!sel.v, drop = FALSE] %*% v.cur[!sel.v,, drop = FALSE]
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Orthogonalize Y using economic QR decomposition: Y=QR
      #If q > 0 perfrom q subspace iterations
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if( reps > 1 ) {
        for( i in 1:reps) {
          u.cur <- qr.Q( qr(u.cur, complete = FALSE) , complete = FALSE )
          Z <- crossprod_help(x.scaled , u.cur )
          Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
          u.cur <- u.cur %*% Z
        }#End for
        remove(Z)
      }#End if
      
      #x.mul = x.scaled[,!sel.v, drop = FALSE]
      #if (reps > 1){
      #  ### increase signal
      #  svd_result <- svd(x.mul)
      #  x.mul <- svd_result$u %*% diag(svd_result$d^reps) %*% t(svd_result$v)
      #}
      
      #u.cur <- x.mul %*% v.cur[!sel.v,, drop = FALSE]
      
      # thresholding
      if (method == "theory") {
        cur <-  thresh(u.cur, gamma.u * sqrt(log(pu)))
        u.cur <- cur$x
        rows <- cur$I.rows
      }
      else {
        threshold.u <- error.est(x.scaled, u.old, v.cur,
                                 n.boot = n.boot)
        cur <- thresh(u.cur, threshold.u)
        u.cur <- cur$x
        rows <- cur$I.rows
      }
      if (all(u.cur == 0)) {
        # every entry is zero
        warning("SSVD.iter.thresh: Overthresh: all zero!")
        u.cur <- u.old
        dist.u <- 0
        dist.v <- 0
        break
      }
      else {
        # QR decomposition
        #u.cur <- qr.Q(qr(u.cur))
        #selected_col = which(apply(u.cur^2, 2, sum) >0)
        u.cur = solve(Sigma_u[rows, rows], u.cur)
        
        norm_u = t(u.cur) %*% Sigma_u %*% u.cur
        index_zero = which(diag(norm_u) < 1e-8)
        print(index_zero)
        for (ind in index_zero){
          norm_u[ind, ind] = 1
        }
        u.cur <- u.cur %*% (sqrtm(norm_u)$Binv)
        dist.u <- subsp.dist.orth(u.cur, u.old)
        #print(c("Dim u", dim(u.cur)))
      }
      v.old <- v.cur
      # multiplication
      sel.u <- apply(u.cur == 0, 1, all)
      
      t_x.mul = t(x.scaled[!sel.u,, drop = FALSE]) 
      if (reps > 1){
        svd_result <- svd(t_x.mul)
        x.mul <- svd_result$u %*% diag(svd_result$d^reps) %*% t(svd_result$v)
      }
      v.cur <- t_x.mul %*% u.cur[!sel.u,, drop = FALSE]
      # thresholding
      if (method == "theory") {
        cur <- thresh(v.cur, gamma.v * sqrt(log(pv)))
        v.cur <- cur$x
        cols <- cur$I.cols
      }
      else {
        threshold.v <- error.est(t(x.scaled), v.old, u.cur,
                                 n.boot = n.boot)
        cur <- thresh(v.cur, threshold.v)
        v.cur <- cur$x
        cols <- cur$I.cols
      }
      if (all(v.cur == 0)) {
        # every entry is zero
        warning("SSVD.iter.thresh: Overthresh: all zero!")
        v.cur <- v.old
        dist.u <- 0
        dist.v <- 0
        break
      }
      else {
        # QR decomposition
        #v.cur <- qr.Q(qr(v.cur))
        ### check if some rows are 0
        v.cur = solve(Sigma_v[cols, cols], v.cur)
        norm_v = t(v.cur) %*% Sigma_v %*% v.cur
        index_zero = which(diag(norm_v) < 1e-8)
        for (ind in index_zero){
          norm_v[ind, ind] = 1
        }
        v.cur <- v.cur %*% (sqrtm(norm_v)$Binv)
        dist.v <- subsp.dist.orth(v.cur, v.old)
        #print(c("Dim v", dim(v.cur)))
        #print(c("Dim u after v", dim(u.cur)))
      }
      i.iter <- i.iter + 1
      #print("iteration is")
      #print(i.iter)
      #print("dim U")
      #print(dim(u.cur))
    }
    
    if(non.orth == TRUE){
      u.old <- u.cur
      u.cur <- x.scaled %*% v.cur
      # thresholding
      if (method == "theory") {
        u.cur <- thresh(u.cur, gamma.u * sqrt(log(pu)))
      }
      else {
        threshold.u <- error.est(x.scaled, u.old, v.cur,
                                 n.boot = n.boot)
        u.cur <- thresh(u.cur, threshold.u)
      }
      u.cur <- apply(u.cur, 2, function(x){x/sqrt(sum(x^2))})
      v.old <- v.cur
      # multiplication
      v.cur <- t(x.scaled) %*% u.old
      # thresholding
      if (method == "theory") {
        v.cur <- thresh(v.cur, gamma.v * sqrt(log(pv)))
      }
      else {
        threshold.v <- error.est(t(x.scaled), v.old, u.cur,
                                 n.boot = n.boot)
        v.cur <- thresh(v.cur, threshold.v)
      }
      v.cur <- apply(v.cur, 2, function(x){x/sqrt(sum(x^2))})
    }
    d.cur <- diag(t(u.cur) %*% x %*% v.cur)
    # To make singular values positive
    u.cur <- u.cur %*% diag(sign(d.cur), r, r)
    if (i.iter == n.iter) {
      warning("increase n.iter")
    }
    if (non.orth == TRUE){
      return(list(u = u.cur, v = v.cur,
                  u.orth = u.old %*% diag(sign(apply(u.cur*u.old, 2, sum)), r, r),
                  v.orth = v.old %*% diag(sign(apply(v.cur*v.old, 2, sum)), r, r),
                  d = abs(d.cur), niter = i.iter - 1,
                  sigma.hat = sigma.hat, dist.u = dist.u, dist.v = dist.v))
    } else{
      return(list(u = u.cur, v = v.cur, d = abs(d.cur), niter = i.iter - 1,
                  sigma.hat = sigma.hat, dist.u = dist.u, dist.v = dist.v))
    }
  }


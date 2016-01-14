#' Simulate data
#'
#' @param n Sample size.
#' @param rho Pairwise correlation between covariates.
#' @param theta Standard deviation of the random error.
#' @param binary If TRUE, generate binary responses; otherwise, by default,
#' create continuous responses.
#' @return gen.data returns a list containning at least the following
#' components:
#' "X", a covariate matrix of n observations and p predictors;
#' "y", a univariate response;
#' "b.true", the actual coefficients for each predictor group.
#' @details
#' This function simulates data as presented in Liu (2015).
#' @references Liu, Y. (2015). Approaches to reduce and integrate data in
#' structured and high-dimensional regression problems in Genomics. Ph.D.
#' Dissertation, The Pennsylvania State University, University Park,
#' Department of Statistics.
#' @examples
#' data <- gen.data(n=100)
#' names(data)
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm

gen.data<-function(n, rho=0.5, theta=1, binary=FALSE)
{
    ps<-c(5, 10)
    ngrp<-length(ps)
    p<-sum(ps)

    b1<-c(1, -1, 0, 0, rep(0, ps[1]-4))
    b2<-c(1, 1, -1, -1, rep(0, ps[2]-4))

    b.true<-list(cbind(b1), cbind(b2))

    Gamma=diag(p)
    covx= matrix(rho,nrow=p,ncol=p)+diag(1-rho,p)
    X=mvrnorm(n,numeric(p),covx)
    epsilon=rnorm(n,0,1)

    V1=X[,1:ps[1]]
    V2=X[,(ps[1]+1):(ps[1]+ps[2])]

    if (binary){
        y = as.numeric(exp(0.8*V1%*%b1)+(2*V2%*%b2)+theta*epsilon > 1)
    }else{
        y = exp(0.8*V1%*%b1)+(2*V2%*%b2)+theta*epsilon
    }

    ans<-list(X=X, y=y, b.true=b.true)
    return(ans)
}

#' Center a vector
#'
#' @param v A vector.
#' @return A vector with mean zero.
#' @details
#' This function centers any vector and returns a vector with mean zero.
#' @examples
#' data <- gen.data(n=100)
#' y.centered <- center(data$y)
#' @export

center<-function(v){
    v - mean(v)
}

#' Covariance matrix
#'
#' @param X a n x p matrix of n observations and p predictors.
#' @return A p x p covariance matrix.
#' @details
#' This function returns A p x p covariance matrix for any n x p matrix.
#' @examples
#' data <- gen.data(n=100)
#' x.cov <- cov.x(data$X)
#' @export

cov.x<-function(X)
{
    Xc<-apply(X, 2, center)
    t(Xc) %*% Xc / nrow(Xc)
}

#' Normalize a vector
#'
#' @param v A vector.
#' @return A vector with norm 1.
#' @details
#' This function normalizes any non-zero vector and returns a vector with
#' the norm equal to 1.
#' @examples
#' data <- gen.data(n=100)
#' y.norm1 <- norm1(data$y)
#' @export

norm1<-function(v)
{
    sumv2<-sum(v^2)
    if(sumv2 == 0) sumv2<-1
    v/sqrt(sumv2)
}

#' Gram-Schmidt orthonormalization
#'
#' @param X a n x p matrix of n observations and p predictors.
#' @return A n x p matrix of n observations and p predictors.
#' @details
#' This function orthonormalizes any n x p matrix.
#' @examples
#' data <- gen.data(n=100)
#' x.orth <- orthnormal(data$X)
#' @export

orthnormal<-function(X)
{
    X<-as.matrix(X)
    n<-nrow(X)
    p<-ncol(X)

    W<-NULL
    if(p > 1) {
        W<-cbind(W, X[,1])
        for(k in 2:p) {
            gw<-rep(0, n)
            for(i in 1:(k-1)) {
                gki<-as.vector((t(W[,i]) %*% X[,k])/(t(W[,i]) %*% W[,i]))
                gw<-gw + gki * W[,i]
            }
            W<-cbind(W, X[,k] - gw)
        }
    } else {
        W<-cbind(W, X[,1])
    }

    W<-apply(W, 2, norm1)
    W
}

#' Power of a matrix
#'
#' @param X A p x p square matrix.
#' @param alpha A scaler determining the order of the power.
#' @return A p x p square matrix.
#' @details
#' This function calculates the power of a square matrix.
#' @examples
#' data <- gen.data(n=100)
#' cov.squared <- matpower(cov.x(data$X), 2)
#' @export

matpower = function(X,alpha){
    X = (X + t(X))/2
    tmp = eigen(X)
    return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
               t(tmp$vectors))}

#' Matrix standardization
#'
#' @param x A n x p matrix of n observations and p predictors.
#' @return A n x p matrix of n observations and p predictors.
#' @details
#' This function standardizes a matrix treating each row as a random vector
#' in an iid sample. It returns a n x p matrix with column-mean zero
#' and identity-covariance matrix.
#' @examples
#' data <- gen.data(n=100)
#' x.std <- standmat(data$X)
#' @export
#' @importFrom stats var

standmat = function(x){
    mu = apply(x,2,mean)
    sig = var(x)
    signrt = matpower(sig,-1/2)
    return(t(t(x) - mu)%*%signrt)
}

#' Subspace distance
#'
#' @param v1 A matrix, each column consists of a p-dimensional vector.
#' @param v2 A matrix, each column consists of a p-dimensional vector.
#' @return A scaler represents the distance between the two spaces spanned by
#' v1 and v2 respectively.
#' @details
#' This function computes the distances between two spaces using the formulation
#' in Li, Zha, Chiaromonte (2005), which is the Frobenius norm of the difference
#' between the two orthogonal projection matrices defined by v1 and v2.
#' @references Li, B., Zha, H., and Chiaromonte, F. (2005). Contour regression:
#' a general approach to dimension reduction. Annals of Statistics,
#' 33(4):1580-1616.
#' @examples
#' v1 <- c(1, 0, 0)
#' v2 <- c(0, 1, 0)
#' disvm(v1, v1)
#' disvm(v1, v2)
#' @export

disvm = function(v1,v2){
    v1 <- as.matrix(v1)
    v2 <- as.matrix(v2)
    if (dim(v1)[2] == 1){
        p1 = v1%*%t(v1)/c(t(v1)%*%v1)
    }else{
        p1 <- v1%*%matpower(t(v1)%*%v1,-1)%*%t(v1)
    }
    if (dim(v2)[2] == 1){
        p2 = v2%*%t(v2)/c(t(v2)%*%v2)
    }else{
        p2 <- v2%*%matpower(t(v2)%*%v2,-1)%*%t(v2)
    }
    d = sqrt(sum((p1-p2)*(p1-p2)))
    return(d)
}

#' Groupwise OLS (gOLS)
#'
#' @param X A covariate matrix of n observations and p predictors.
#' @param Y A univariate response.
#' @param groups A vector with the number of predictors in each group.
#' @param Gamma A p x p identity matrix.
#' @param dims A vector with the dimension (at most 1) for each predictor group.
#' @param tol A tolerence factor in the convergence criterium.
#' @param maxiter A maximum iteration number.
#' @return gOLS returns a list containning at least the following components:
#' "b_init", the initial estimatied directions for each group with its own
#' dimension using assembled OLS;
#' "b_est", the estimated directions for each group with its own dimension
#' using gOLS AFTER normalization;
#' "B", the estimated directions for each group using gOLS BEFORE normalization;
#' "t_est", the estimated p-dimensional directions for each group using gOLS
#' AFTER normalization;
#' "obj", the objective functions from each iteration;
#' "niter", the number of iterations.
#' @details
#' This function estimates directions for each predictor group using gOLS.
#' Predictors need to be organized in groups within the "X" matrix, as the
#' same order saved in "groups". We only allow continuous covariates
#' in the "X" matrix; while categorical covariates can be handled outside of
#' gOLS, e.g. structured OLS.
#' "b_est", "B", and "t_est" are basically the same estimation with different
#' layouts and/or normalizations.
#' @references Liu, Y., Chiaromonte, F., and Li, B. (2015). Structured Ordinary
#' Least Squares: a sufficient dimension reduction approach for regressions with
#'  partitioned predictors and heterogeneous units. Submitted.
#' @examples
#' data <- gen.data(n=1000, binary=FALSE) # generate data
#' dim(data$X) # covariate matrix of 1000 observations and 15 predictors
#' dim(data$y) # univariate response
#' groups <- c(5, 10) # two predictor groups and their numbers of predictors
#' Gamma <- diag(sum(groups)) # identity matrix
#' dims <- c(1,1) # dimension of each predictor group
#' est_gOLS <- gOLS(data$X,data$y,groups,Gamma,dims,tol=1e-5,maxiter=100)
#' names(est_gOLS)
#' @export
#' @importFrom MASS ginv
#' @importFrom Matrix bdiag
#' @importFrom stats cov

gOLS = function(X,Y,groups,Gamma,dims,tol=1e-5,maxiter=100){
    #input :X, Y, groups, Gamma, dims, tol, maxiter
    #output: b_est: estimated direction after gamma
    #        t_est: estimated directions for each group
    #        obj: objective functions from each iteration
    #        niter: number of iterations
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n=nrow(X)
    p=ncol(X)
    ngroup=length(groups)
    ggamma=NULL
    #create a list of matrices store group information (gamma's)
    for (g in 1:ngroup) {
        first=sum(groups[0:(g-1)])+1
        last=sum(groups[1:g])
        tmp= Gamma[,first:last]
        if (is.null(ncol(tmp))) tmp=matrix(tmp,nrow=p)
        ggamma=c(ggamma,list(tmp))
    }

    #Center X, estimate Sigma_x, Center Y
    Xc=apply(X,2,center)
    Sigma_x=cov(Xc)/n*(n-1)

    #   Sigma_inv=solve(Sigma_x)    #conventional inverse
    Sigma_inv=ginv(Sigma_x)   #Moore-Penrose generalized inverse

    eig=eigen(Sigma_x)
    Sigma_half=eig$vectors%*%diag(sqrt(round(eig$values,10)))%*%t(eig$vectors)
    inveig=eigen(Sigma_inv)
    Sigmainv_half=
        inveig$vectors%*%diag(sqrt(round(inveig$values,10)))%*%t(inveig$vectors)

    Yc=apply(Y,2,center)

    #Estimate U and V
    U <- Sigma_inv %*% cov(Xc, Yc)
    # U <- Sigma_inv %*% t(Xc) %*% Yc / n
    V <- Sigma_half %*% U

    #Obtain initial estimator of (b1,...,bg)
    b_init=NULL
    for (g in 1:ngroup){
        Xg=Xc%*%ggamma[[g]]
        out.ols=ginv(t(Xg)%*%Xg)%*%t(Xg)%*%Yc
        dg=dims[g]
        bg_init=as.matrix(out.ols)
        b_init=c(b_init,list(bg_init))
    }
    b_est=b_init
    #estimate B based on initial esimates
    B=b_est[[1]]
    if (ngroup>1) {
        for (g in 2:ngroup) { B=bdiag(B,b_est[[g]])}
    }
    obj_diff=1    #initial value for the difference of objective function
    obj_rec=NULL  #recording the ojective function from each iteration
    k=0

    #start the iteration, until convergence criterions met
    while (obj_diff>tol && k<maxiter) {
        k=k+1
        #Obtain estimate of C
        #estimate A
        A <- as.matrix(Sigma_half %*% B)
        #esitmate C
        C <- as.matrix(ginv(t(A)%*%A)%*%t(A)%*%V)

        #estimate V2
        V2=NULL
        for (g in 1:ngroup) {
            first=sum(dims[0:(g-1)])+1
            last=sum(dims[1:g])
            fgs=C[first:last,]
            if (is.null(ncol(fgs))) fgs=matrix(fgs,ncol=1)
            V2=cbind(V2, kronecker(t(fgs),Sigma_half%*%ggamma[[g]]))
        }

        #estimate b1,..bg
        vec_b=as.matrix(ginv(t(V2)%*%V2)%*%t(V2)%*%V)

        b_est=NULL
        for (g in 1:ngroup){
            pg=groups[g]
            dg=dims[g]
            first=sum((groups*dims)[0:(g-1)])+1
            last=sum((groups*dims)[1:g])
            vec_bg=vec_b[first:last]
            bg_est=matrix(vec_bg,nrow=pg,ncol=dg,byrow=FALSE)
            b_est=c(b_est,list(bg_est))
        }

        #update B matrix
        B=b_est[[1]]
        if (ngroup>1) {
            for (g in 2:ngroup) { B=bdiag(B,b_est[[g]])}
        }

        #Calculate objective function
        objfn=sum((Sigma_half%*%U-Sigma_half%*%Gamma%*%B%*%C)^2)
        obj_rec=c(obj_rec,objfn)
        if (k==1) {obj_diff=1} else {
            obj_diff=abs(obj_rec[k-1]-objfn)
        }
    }

    #orthnormal b_est
    for (g in 1:ngroup){
        b_est[[g]]=orthnormal(b_est[[g]])
    }

    #output the results
    t_est=NULL
    for ( g in 1:ngroup){
        tg=ggamma[[g]]%*%b_est[[g]]
        t_est=c(t_est,list(tg))
    }

    ans=list(b_init=b_init,b_est=b_est,t_est=t_est,obj=obj_rec,niter=k,B=B)
    return(ans)
}

#' Groupwise SIR (gSIR) for binary response
#'
#' @param X A covariate matrix of n observations and p predictors.
#' @param Y A binary response.
#' @param groups A vector with the number of predictors in each group.
#' @param Gamma A p x p identity matrix.
#' @param dims A vector with the dimension (at most 1) for each predictor group.
#' @param tol A tolerence factor in the convergence criterium.
#' @param maxiter A maximum iteration number.
#' @return gSIR returns a list containning at least the following components:
#' "b_init", the initial estimatied directions for each group with its own
#' dimension using assembled SIR;
#' "b_est", the estimated directions for each group with its own dimension
#' using gSIR AFTER normalization;
#' "B", the estimated directions for each group using gSIR BEFORE normalization;
#' "t_est", the estimated p-dimensional directions for each group using gSIR
#' AFTER normalization;
#' "obj", the objective functions from each iteration;
#' "niter", the number of iterations.
#' @details
#' This function estimates directions for each predictor group using gSIR.
#' Predictors need to be organized in groups within the "X" matrix, as the
#' same order saved in "groups". We only allow continuous covariates
#' in the "X" matrix; while categorical covariates can be handled outside of
#' gSIR, e.g. structured SIR.
#' "b_est", "B", and "t_est" are basically the same estimation with different
#' layouts and/or normalizations.
#' @references Guo, Z., Li, L., Lu, W., and Li, B. (2014). Groupwise dimension
#' reduction via envelope method. Journal of the American Statistical
#' Association, accepted.
#' @examples
#' data <- gen.data(n=1000, binary=TRUE) # generate data
#' dim(data$X) # covariate matrix of 1000 observations and 15 predictors
#' length(data$y) # binary response
#' groups <- c(5, 10) # two predictor groups and their numbers of predictors
#' Gamma <- diag(sum(groups)) # identity matrix
#' dims <- c(1,1) # dimension of each predictor group
#' est_gSIR<-gSIR(data$X,data$y,groups,Gamma,dims,tol=1e-5,maxiter=100)
#' names(est_gSIR)
#' @export
#' @importFrom MASS ginv
#' @importFrom Matrix bdiag
#' @importFrom stats cov

gSIR = function(X,Y,groups,Gamma,dims,tol=1e-5,maxiter=100){
    #input :X, Y, groups, Gamma, dims, tol, maxiter
    #output: b_est: estimated direction after gamma
    #        t_est: estimated directions for each group
    #        obj: objective functions from each iteration
    #        niter: number of iterations
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n=nrow(X)
    p=ncol(X)
    ngroup=length(groups)
    ggamma=NULL
    #create a list of matrices store group information (gamma's)
    for (g in 1:ngroup) {
        first=sum(groups[0:(g-1)])+1
        last=sum(groups[1:g])
        tmp= Gamma[,first:last]
        if (is.null(ncol(tmp))) tmp=matrix(tmp,nrow=p)
        ggamma=c(ggamma,list(tmp))
    }

    #Center X, estimate Sigma_x, Center Y
    Xc=apply(X,2,center)
    Yc <- as.factor(Y)
    Y1 <- which(Yc == levels(Yc)[1])
    Y2 <- which(Yc == levels(Yc)[2])

    Sigma_x=cov(Xc)/n*(n-1)

    #   Sigma_inv=solve(Sigma_x)    #conventional inverse
    Sigma_inv=ginv(Sigma_x)   #Moore-Penrose generalized inverse

    eig=eigen(Sigma_x)
    Sigma_half=eig$vectors%*%diag(sqrt(round(eig$values,10)))%*%t(eig$vectors)
    inveig=eigen(Sigma_inv)
    Sigmainv_half=
        inveig$vectors%*%diag(sqrt(round(inveig$values,10)))%*%t(inveig$vectors)

    #Estimate U and V
    U <- Sigma_inv %*% (apply(Xc[Y2,], 2, mean) - apply(Xc[Y1,], 2, mean))
    #   U <- Sigma_inv %*% cov(Xc, Yc)
    #   U <- Sigma_inv %*% t(Xc) %*% Yc / n
    V <- Sigma_half %*% U

    #Obtain initial estimator of (b1,...,bg)
    b_init=NULL
    for (g in 1:ngroup){
        Xg=Xc%*%ggamma[[g]]
        out.sir=
            ginv(cov(Xg)/n*(n-1))%*%
            (apply(Xg[Y2,], 2, mean) - apply(Xg[Y1,], 2, mean))
        dg=dims[g]
        bg_init=as.matrix(out.sir)
        b_init=c(b_init,list(bg_init))
    }
    b_est=b_init
    #estimate B based on initial esimates
    B=b_est[[1]]
    if (ngroup>1) {
        for (g in 2:ngroup) { B=bdiag(B,b_est[[g]])}
    }
    obj_diff=1    #initial value for the difference of objective function
    obj_rec=NULL  #recording the ojective function from each iteration
    k=0

    #start the iteration, until convergence criterions met
    while (obj_diff>tol && k<maxiter) {
        k=k+1
        #Obtain estimate of C
        #estimate A
        A <- as.matrix(Sigma_half %*% B)
        #esitmate C
        C <- as.matrix(ginv(t(A)%*%A)%*%t(A)%*%V)

        #estimate V2
        V2=NULL
        for (g in 1:ngroup) {
            first=sum(dims[0:(g-1)])+1
            last=sum(dims[1:g])
            fgs=C[first:last,]
            if (is.null(ncol(fgs))) fgs=matrix(fgs,ncol=1)
            V2=cbind(V2, kronecker(t(fgs),Sigma_half%*%ggamma[[g]]))
        }

        #estimate b1,..bg
        vec_b=as.matrix(ginv(t(V2)%*%V2)%*%t(V2)%*%V)

        b_est=NULL
        for (g in 1:ngroup){
            pg=groups[g]
            dg=dims[g]
            first=sum((groups*dims)[0:(g-1)])+1
            last=sum((groups*dims)[1:g])
            vec_bg=vec_b[first:last]
            bg_est=matrix(vec_bg,nrow=pg,ncol=dg,byrow=FALSE)
            b_est=c(b_est,list(bg_est))
        }

        #update B matrix
        B=b_est[[1]]
        if (ngroup>1) {
            for (g in 2:ngroup) { B=bdiag(B,b_est[[g]])}
        }

        #Calculate objective function
        objfn=sum((Sigma_half%*%U-Sigma_half%*%Gamma%*%B%*%C)^2)
        obj_rec=c(obj_rec,objfn)
        if (k==1) {obj_diff=1} else {
            obj_diff=abs(obj_rec[k-1]-objfn)
        }
    }

    #orthnormal b_est
    for (g in 1:ngroup){
        b_est[[g]]=orthnormal(b_est[[g]])
    }

    #output the results
    t_est=NULL
    for ( g in 1:ngroup){
        tg=ggamma[[g]]%*%b_est[[g]]
        t_est=c(t_est,list(tg))
    }

    ans=list(b_init=b_init,b_est=b_est,t_est=t_est,obj=obj_rec,niter=k,B=B)
    return(ans)
}

#' Groupwise OLS (gOLS) BIC criterion to estimate dimensions with
#' eigen-decomposition
#'
#' @param X A covariate matrix of n observations and p predictors.
#' @param y A univariate response.
#' @param groups A vector with the number of predictors in each group.
#' @param Gamma A p x p identity matrix.
#' @return gOLS.comp.d returns a list containning at least the following
#' components:
#' "d", the estimated dimension (at most 1) for each predictor group;
#' "crit", the BIC criterion from each iteration.
#' @details
#' This function estimates dimension for each predictor group using
#' eigen-decomposition. Predictors need to be organized in groups within the
#' "X" matrix, as the same order saved in "groups". We only allow continuous
#' covariates in the "X" matrix; while categorical covariates can be handled
#' outside of gOLS, e.g. structured OLS.
#' @references Liu, Y., Chiaromonte, F., and Li, B. (2015). Structured Ordinary
#' Least Squares: a sufficient dimension reduction approach for regressions with
#'  partitioned predictors and heterogeneous units. Submitted.
#' @examples
#' data <- gen.data(n=1000, binary=FALSE) # generate data
#' dim(data$X) # covariate matrix of 1000 observations and 15 predictors
#' dim(data$y) # univariate response
#' groups <- c(5, 10) # two predictor groups and their numbers of predictors
#' Gamma <- diag(sum(groups)) # identity matrix
#' dim_gOLS<-gOLS.comp.d(data$X,data$y,groups,Gamma)
#' names(dim_gOLS)
#' @export
#' @importFrom stats sd

gOLS.comp.d<-function(X, y, groups, Gamma)
{
    n<-nrow(X)
    ngrp<-length(groups)

    dsw<-rep(1, ngrp)
    d<-rep(0, ngrp)

    Z <- standmat(X)
    Y <- (y-mean(y))/sd(y)
    out.d<-gOLS(Z,Y,groups,Gamma,dsw,tol=1e-6,maxiter=500)

    M <- t(as.matrix(out.d$B)) %*% (as.matrix(out.d$B))
    vecM <- diag(M)
    orderM <- order(vecM, decreasing = TRUE)

    crit <- -1/(n^(1/8)*log(n))
    for (i in 1:length(vecM)){
        crit.d <- sum(vecM[orderM][1:i]) - (i+1)/(n^(1/8)*log(n))
        crit<-c(crit, crit.d)
    }

    if (which(crit == max(crit)) > 1){
        d[orderM[1:(which(crit == max(crit))-1)]] <- 1
    }

    ans<-list(d=d, crit=crit)
    return(ans)
}

#' Groupwise SIR (gSIR) BIC criterion to estimate dimensions with
#' eigen-decomposition (binary response)
#'
#' @param X A covariate matrix of n observations and p predictors.
#' @param y A binary response.
#' @param groups A vector with the number of predictors in each group.
#' @param Gamma A p x p identity matrix.
#' @return gSIR.comp.d returns a list containning at least the following
#' components:
#' "d", the estimated dimension (at most 1) for each predictor group;
#' "crit", the BIC criterion from each iteration.
#' @details
#' This function estimates dimension for each predictor group using
#' eigen-decomposition. Predictors need to be organized in groups within the
#' "X" matrix, as the same order saved in "groups". We only allow continuous
#' covariates in the "X" matrix; while categorical covariates can be handled
#' outside of gSIR, e.g. structured SIR.
#' @references Liu, Y. (2015). Approaches to reduce and integrate data in
#' structured and high-dimensional regression problems in Genomics. Ph.D.
#' Dissertation, The Pennsylvania State University, University Park,
#' Department of Statistics.
#' @examples
#' data <- gen.data(n=1000, binary=TRUE) # generate data
#' dim(data$X) # covariate matrix of 1000 observations and 15 predictors
#' length(data$y) # univariate response
#' groups <- c(5, 10) # two predictor groups and their numbers of predictors
#' Gamma <- diag(sum(groups)) # identity matrix
#' dim_gSIR<-gSIR.comp.d(data$X,data$y,groups,Gamma)
#' names(dim_gSIR)
#' @export

gSIR.comp.d<-function(X, y, groups, Gamma)
{
    n<-nrow(X)
    ngrp<-length(groups)

    dsw<-rep(1, ngrp)
    d<-rep(0, ngrp)

    Z <- standmat(X)
    Y <- y
    out.d<-gSIR(Z,Y,groups,Gamma,dsw,tol=1e-6,maxiter=500)

    M <- t(as.matrix(out.d$B)) %*% (as.matrix(out.d$B))
    vecM <- diag(M)
    orderM <- order(vecM, decreasing = TRUE)

    crit <- -1/(n^(1/16)*log(n))
    for (i in 1:length(vecM)){
        crit.d <- sum(vecM[orderM][1:i]) - (i+1)/(n^(1/16)*log(n))
        crit<-c(crit, crit.d)
    }

    if (which(crit == max(crit)) > 1){
        d[orderM[1:(which(crit == max(crit))-1)]] <- 1
    }

    ans<-list(d=d, crit=crit)
    return(ans)
}

#' Structured OLS (sOLS) outer level BIC criterion to estimate dimension with
#' eigen-decomposition
#'
#' @param X A matrix containing directions estimated from all subpopulations.
#' @param sizes A vector with the sample sizes of all subpopulation.
#' @return sOLS.comp.d returns the dimension estimated across subpopulations.
#' @details
#' This function estimates dimension across the subpopulations using
#' eigen-decomposition. The order of the subpopulations in the "sizes" vector
#' should match the one in the "X" matrix.
#' @references Liu, Y., Chiaromonte, F., and Li, B. (2015). Structured Ordinary
#' Least Squares: a sufficient dimension reduction approach for regressions with
#'  partitioned predictors and heterogeneous units. Submitted.
#' @examples
#' v1 <- c(1, 1, 0, 0)
#' v2 <- c(0, 1, 1, 0)
#' v3 <- c(0, 0, 1, 1)
#' v4 <- c(1, 1, 1, 1)
#' m1 <- cbind(v1, v2)
#' sizes1 <- c(100, 200)
#' sOLS.comp.d(m1, sizes1)
#' m2 <- cbind(v1, v2, v3)
#' sizes2 <- c(100, 200, 500)
#' sOLS.comp.d(m2, sizes2)
#' m3 <- cbind(v1, v3, v4)
#' sizes3 <- c(100, 500, 1000)
#' sOLS.comp.d(m3, sizes3)
#' @export

sOLS.comp.d<-function(X, sizes)
{
    MM <- X %*% t(X)
    crit <- NULL
    for (w in 1:min(dim(X))){
        crit.d <- sum(eigen(MM)$values[1:w]) - w/(min(sizes)^(1/8))
        crit<-c(crit, crit.d)
    }
    ans <- which(crit == max(crit))
    return(ans)
}

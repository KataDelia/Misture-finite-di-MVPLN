customvga <- function(dataset,
                         membership = "none",
                         gmin,
                         gmax,
                         initMethod = "kmeans",
                         nInitIterations = 0,
                         normalize = "Yes") {
    
    ptm <- base::proc.time()
    
    # Performing checks
    if (typeof(unlist(dataset)) != "double" & typeof(unlist(dataset)) != "integer") {
        stop("dataset should be a list of count matrices.");}
    
    if (any((unlist(dataset) %% 1 == 0) == FALSE)) {
        stop("dataset should be a list of count matrices.")
    }
    
    if (is.list(dataset) != TRUE) {
        stop("dataset needs to be a list of matrices.")
    }
    
    if(is.numeric(gmin) != TRUE || is.numeric(gmax) != TRUE) {
        stop("Class of gmin and gmin should be numeric.")
    }
    
    if (gmax < gmin) {
        stop("gmax cannot be less than gmin.")
    }
    
    if (gmax > length(dataset)) {
        stop("gmax cannot be larger than nrow(dataset).")
    }
    
    
    if(all(membership != "none") && is.numeric(membership) != TRUE) {
        stop("membership should be a numeric vector containing the
      cluster membership. Otherwise, leave as 'none'.")
    }
    
    # Checking if missing membership values
    # First check for the case in which G = 1, otherwise check
    #   if missing cluster memberships
    if(all(membership != "none") && length(unique(membership)) != 1) {
        if(all(membership != "none") &&
           all((diff(sort(unique(membership))) == 1) != TRUE) ) {
            stop("Cluster memberships in the membership vector
        are missing a cluster, e.g. 1, 3, 4, 5, 6 is missing cluster 2.")
        }
    }
    
    if(all(membership != "none") && length(membership) != length(dataset)) {
        stop("membership should be a numeric vector, where length(membership)
      should equal the number of observations. Otherwise, leave as 'none'.")
    }
    
    if (is.character(initMethod) == TRUE) {
        if(initMethod != "none" & initMethod != "kmeans" & initMethod != "random" & initMethod != "medoids" & initMethod != "medoids" & initMethod != "clara" & initMethod != "fanny") {
            stop("initMethod should of class character, specifying
      either: none, kmeans, random, medoids, clara, or fanny.")
        }
    } else if (is.character(initMethod) != TRUE) {
        stop("initMethod should of class character, specifying
      either: none, kmeans, random, medoids, clara, or fanny.")
    }
    
    if (is.numeric(nInitIterations) != TRUE) {
        stop("nInitIterations should be positive integer or zero, specifying
      the number of initialization runs to be considered.")
    }
    
    if (is.character(normalize) != TRUE) {
        stop("normalize should be a string of class character specifying
      if normalization should be performed.")
    }
    
    if (normalize != "Yes" && normalize != "No") {
        stop("normalize should be a string indicating Yes or No, specifying
       if normalization should be performed.")
    }
    
    
    n <- length(dataset)
    p <- ncol(dataset[[1]])
    r <- nrow(dataset[[1]])
    d <- p * r
    
    
    # Changing dataset into a n x rp dataset
    TwoDdataset <- matrix(NA, ncol = d, nrow = n)
    sampleMatrix <- matrix(c(1:d), nrow = r, byrow = TRUE)
    for (u in 1:n) {
        for (e in 1:p) {
            for (s in 1:r) {
                TwoDdataset[u, sampleMatrix[s, e]] <- dataset[[u]][s, e]
            }
        }
    }
    
    # Check if entire row is zero
    if(any(rowSums(TwoDdataset) == 0)) {
        TwoDdataset <- TwoDdataset[- which(rowSums(TwoDdataset) == 0), ]
        n <- nrow(TwoDdataset)
        membership <- membership[- which(rowSums(TwoDdataset) == 0)]
    }
    
    if(all(is.na(membership) == TRUE)) {
        membership <- "Not provided" }
    
    # Calculating normalization factors
    if(normalize == "Yes") {
        normFactors <- as.vector(edgeR::calcNormFactors(as.matrix(TwoDdataset),
                                                        method = "TMM"))
        
    } else if(normalize == "No") {
        normFactors <- rep(0, d)
    } else {
        stop("Argument normalize should be 'Yes' or 'No' ")
    }
    
    parallelFA <- function(G, dataset,
                           TwoDdataset,
                           r, p, d, n,
                           normFactors,
                           nInitIterations,
                           initMethod) {
        
        # arranging normalization factors
        libMat <- matrix(normFactors, n, d, byrow = T)
        libMatList <- list()
        for (i in 1:n) {
            libMatList[[i]] <- t(matrix(libMat[i, ], nrow = p))
        }
        
        if (initMethod == "none") {
            
            # Initialization
            mu <- omega <- phi <- list() # mu is M;
            delta <- kappa <- sigma <- isigma <- iphi <- iomega <- list()
            # delta  is variational parameter Delta
            # kappa is variational parameter kappa
            # sigma is Psi in the math - kronecker of omega and Phi
            # isigma is inverse of Psi
            # iphi is inverse of Phi
            # iomega is inverse of omega
            
            m <- S <- list()
            # m is vectorized xi
            # S is Psi
            
            # Other intermediate items initialized
            Sk <- array(0, c(d, d, G) )
            start <- GX <- dGX <- zS <- list()
            
            iKappa <- iDelta <- startList <- list()
            # iKappa is inverse of kappa
            # iDelta is inverse of Delta
            
            kMeansResults <- kmeans(log(TwoDdataset+1),
                                    centers = G,
                                    nstart = 50,
                                    iter.max = 20)$cluster
            zValue <- mclust::unmap(kMeansResults) ### Starting value for Z
            piG <- colSums(zValue) / n
            
            for (g in 1:G) {
                obs <- which(zValue[, g] == 1)
                mu[[g]] <- colMeans(log(TwoDdataset[obs, ] + 1 / 6))
                sigma[[g]] <- var(log(TwoDdataset[obs, ] + 1 / 6))
                isigma[[g]] <- solve(sigma[[g]])
                phi[[g]] <- diag(r) * sqrt(min(diag(var(log(TwoDdataset[obs, ] + 1 / 6)))))
                omega[[g]] <- diag(p) * sqrt(min(diag(var(log(TwoDdataset[obs, ] + 1 / 6)))))
                iphi[[g]] <- solve(phi[[g]])
                iomega[[g]] <- solve(omega[[g]])
            }
            
            for (g in 1:G) {
                start[[g]] <- log(TwoDdataset + 1 / 6) ###Starting value for M
                m[[g]] <- log(TwoDdataset + 1 / 6)
                S[[g]] <- list()
                delta[[g]] <- list()
                kappa[[g]] <- list()
                startList[[g]] <- list()
                for (i in 1:n) {
                    startList[[g]][[i]] <- log(dataset[[i]])
                    delta[[g]][[i]] <- diag(r) * 0.001
                    kappa[[g]][[i]] <- diag(p) * 0.001
                    S[[g]][[i]] <- delta[[g]][[i]] %x% kappa[[g]][[i]]
                }
            }
        } else {
            
            # Initialize based on specified initMethod method
            outputInitialization <- list()
            checklogL <- vector()
            for (initIt in seq_along(1:nInitIterations)) {
                set.seed(initIt)
                outputInitialization[[initIt]] <- initializationRunVGA(
                    G = G,
                    dataset = dataset,
                    TwoDdataset = TwoDdataset,
                    r = r,
                    p = p,
                    d = d,
                    n = n,
                    normFactors = normFactors,
                    initMethod = initMethod)
                checklogL[initIt] <- outputInitialization[[initIt]]$finalLogLik
            }
            
            # select init run with highest logL
            maxRun <- which.max(checklogL)
            
            mu <- outputInitialization[[maxRun]]$mu
            omega <- outputInitialization[[maxRun]]$omega
            phi <- outputInitialization[[maxRun]]$phi
            delta <- outputInitialization[[maxRun]]$delta
            kappa <- outputInitialization[[maxRun]]$kappa
            sigma <- outputInitialization[[maxRun]]$sigma
            isigma <- outputInitialization[[maxRun]]$isigma
            iphi <- outputInitialization[[maxRun]]$iphi
            iomega <- outputInitialization[[maxRun]]$iomega
            m <- outputInitialization[[maxRun]]$m
            S <- outputInitialization[[maxRun]]$S
            Sk <- outputInitialization[[maxRun]]$Sk
            start <- outputInitialization[[maxRun]]$start
            GX <- outputInitialization[[maxRun]]$GX
            dGX <- outputInitialization[[maxRun]]$dGX
            zS <- outputInitialization[[maxRun]]$zS
            iKappa <- outputInitialization[[maxRun]]$iKappa
            iDelta <- outputInitialization[[maxRun]]$iDelta
            startList <- outputInitialization[[maxRun]]$startList
            zValue <- outputInitialization[[maxRun]]$zValue
            piG <- outputInitialization[[maxRun]]$piG
        }
        
        # start clustering after initialization
        it <- 1
        aloglik <- loglik <- NULL
        checks <- aloglik[c(1:5)] <- 0
        itMax <- 200
        
        while (checks == 0) {
            
            for (g in 1:G) {
                GX[[g]] <- dGX[[g]] <- zS[[g]] <- list()
                iDelta[[g]] <- iKappa[[g]] <- list()
                deltaO <- delta[[g]]
                kappaO <- kappa[[g]]
                
                for (i in 1:n) {
                    iDelta[[g]][[i]] <- diag(c(t(diag(kappa[[g]][[i]])) %*%
                                                   t(exp(log(libMatList[[i]]) +
                                                             startList[[g]][[i]] +
                                                             0.5 * diag(delta[[g]][[i]]) %*%
                                                             t(diag(kappa[[g]][[i]])))))) +
                        iphi[[g]] * sum(diag(iomega[[g]] %*%
                                                 kappa[[g]][[i]]))
                    delta[[g]][[i]] <- p * solve(iDelta[[g]][[i]])
                    
                    
                    iKappa[[g]][[i]] <- diag(c(t(diag(delta[[g]][[i]])) %*%
                                                   (exp(log(libMatList[[i]]) +
                                                            startList[[g]][[i]] +
                                                            t(0.5 * diag(kappa[[g]][[i]]) %*%
                                                                  t(diag(delta[[g]][[i]]))))))) +
                        iomega[[g]] * sum(diag(iphi[[g]] %*%
                                                   delta[[g]][[i]]))
                    
                    kappa[[g]][[i]] <- solve(iKappa[[g]][[i]])
                    
                    
                    S[[g]][[i]] <- delta[[g]][[i]] %x% kappa[[g]][[i]]
                    zS[[g]][[i]] <- zValue[i, g] * S[[g]][[i]]
                    GX[[g]][[i]] <- TwoDdataset[i, ] -
                        exp(start[[g]][i, ] +
                                log(libMat[i,]) +
                                0.5 * diag(S[[g]][[i]])) -
                        (isigma[[g]]) %*% (start[[g]][i, ] - mu[[g]])
                    m[[g]][i, ] <- start[[g]][i, ] + S[[g]][[i]] %*% GX[[g]][[i]]
                    startList[[g]][[i]] <- t(matrix(m[[g]][i, ], nrow = p))
                }
                start[[g]] <- m[[g]]
                
                mu[[g]] <- colSums(zValue[, g] * m[[g]]) / sum(zValue[, g]) # this is xi
                
                # Updating Sample covariance
                muMat <- t(matrix(mu[[g]], nrow = p))
                phiM <- list()
                for (i in 1:n) {
                    phiM[[i]] <- zValue[i,g]*(startList[[g]][[i]] - muMat) %*%
                        iomega[[g]] %*% t(startList[[g]][[i]] - muMat) +
                        zValue[i, g] * delta[[g]][[i]] * sum(diag(iomega[[g]] %*%
                                                                      kappa[[g]][[i]]))
                }
                phi[[g]] <- Reduce("+", phiM) / sum(zValue[, g] * p)
                iphi[[g]] <- solve(phi[[g]])
                
                omegaM <- list()
                for (i in 1:n) {
                    omegaM[[i]] <- zValue[i, g] * t(startList[[g]][[i]] - muMat) %*%
                        iphi[[g]] %*% (startList[[g]][[i]] - muMat) +
                        zValue[i, g] * kappa[[g]][[i]] * sum(diag(iphi[[g]] %*%
                                                                      delta[[g]][[i]]))
                }
                omega[[g]] <- Reduce("+", omegaM) / sum(zValue[, g] * r)
                iomega[[g]] <- solve(omega[[g]])
                sigma[[g]] <- phi[[g]] %x% omega[[g]]
                isigma[[g]] <- iphi[[g]] %x% iomega[[g]]
            }
            
            piG <- colSums(zValue) / n
            # Internal functions
            funFive <- function(x, y = isigma[[g]]) {
                sum(diag(x %*% y))
            }
            
            FMatrix <- matrix(NA, ncol = G, nrow = n)
            
            for (g in 1:G) {
                two <- rowSums(exp(m[[g]] + log(libMat) +
                                       0.5 * matrix(unlist(lapply(S[[g]], diag)),
                                                    ncol = d, byrow = TRUE)))
                five <- 0.5 * unlist(lapply(S[[g]], funFive))
                six <- 0.5 * log(unlist(lapply(S[[g]], det)))
                FMatrix[, g] <- piG[g] * exp(rowSums(m[[g]] * TwoDdataset) -
                                                 two - rowSums(lfactorial(TwoDdataset)) +
                                                 rowSums(log(libMat) * TwoDdataset) - 0.5 *
                                                 mahalanobis(m[[g]], center = mu[[g]], cov = isigma[[g]],
                                                             inverted = TRUE) - five + six + 0.5 *
                                                 log(det(isigma[[g]])) - d/2)
            }
            
            # loglik[it] <- sum(log(rowSums(FMatrix)))
            # if FMatrix values is 0, then log(0) = -Inf and logL error
            rowSumsFMatrix <- rowSums(FMatrix)
            if(any(rowSumsFMatrix == 0) == TRUE) {
                zeroRows <- which(rowSumsFMatrix == 0)
                rowSumsFMatrix[zeroRows] <- min(rowSumsFMatrix[-zeroRows])
            }
            loglik[it] <- sum(log(rowSumsFMatrix))
            
            
            zValue <- FMatrix / rowSumsFMatrix
            
            if (it <= 5) {
                zValue[zValue == "NaN"] <- 0
            }
            
            if (it > 5) {
                # Aitkaine's stopping criterion
                # Print for check; commented out
                # cat("\n At it:", it ,"loglik[it - 1]:",
                # loglik[it - 1], "- loglik[it - 2]:", loglik[it - 2],
                # "is", loglik[it - 1] - loglik[it - 2], "\n")
                if ((loglik[it - 1] - loglik[it - 2]) == 0) checks <- 1 else {
                    a <- (loglik[it] - loglik[it - 1]) / (loglik[it - 1] - loglik[it - 2])
                    addTo <- (1 / (1 - a) * (loglik[it] - loglik[it - 1]))
                    aloglik[it] <- loglik[it - 1] + addTo
                    if (abs(aloglik[it] - aloglik[it - 1]) < 0.05) {
                        checks <- 1
                    } else {
                        checks <- checks
                    }
                }
            }
            
            it <- it + 1
            if (it == itMax) {
                checks <- 1
            }
            finalPhi <- finalOmega <- list()
            for (g in 1:G) {
                finalPhi[[g]] <- phi[[g]] / diag(phi[[g]])[1]
                finalOmega[[g]] <- omega[[g]] * diag(phi[[g]])[1]
            }
        }
        
        programclust <- mclust::map(zValue)
        
        FinalGResults <- list(mu = matrix(unlist(mu),
                                          ncol = p,
                                          byrow = TRUE),
                              sigma = sigma,
                              phi = finalPhi,
                              omega = finalOmega,
                              probaPost = zValue,
                              loglikelihood = loglik,
                              proportion = piG,
                              clusterlabels = programclust,
                              iterations = it)
        
        class(FinalGResults) <- "FinalGResults"
        
        return(FinalGResults)
    }
    
    parallelFAOutput <- list()
    for(g in seq_along(1:(gmax - gmin + 1))) {
        
        if(length(1:(gmax - gmin + 1)) == gmax) {
            clustersize <- g
        } else if(length(1:(gmax - gmin + 1)) < gmax) {
            clustersize <- seq(gmin, gmax, 1)[g]
        }
        
        parallelFAOutput[[g]] <- parallelFA(
            G = clustersize,
            dataset = dataset,
            TwoDdataset = TwoDdataset,
            r = r,
            p = p,
            d = d,
            n = n,
            normFactors = normFactors,
            nInitIterations = nInitIterations,
            initMethod = initMethod)
    }
    
    
    # cluster data result extracting
    BIC <- ICL <- AIC <- AIC3 <- k <- ll <- vector()
    for(g in seq_along(1:(gmax - gmin + 1))) {
        #print(g)
        if(length(1:(gmax - gmin + 1)) == gmax) {
            clustersize <- g
        } else if(length(1:(gmax - gmin + 1)) < gmax) {
            clustersize <- seq(gmin, gmax, 1)[g]
        }
        
        # save the final log-likelihood
        ll[g] <- unlist(tail(parallelFAOutput[[g]]$loglik, n = 1))
        
        k[g] <- calcParameters(g = clustersize,
                               r = r,
                               p = p)
        
        # starting model selection
        if (g == max(1:(gmax - gmin + 1))) {
            bic <- BICFunction(ll = ll,
                               k = k,
                               n = n,
                               run = parallelFAOutput,
                               gmin = gmin,
                               gmax = gmax,
                               parallel = FALSE)
            
            icl <- ICLFunction(bIc = bic,
                               gmin = gmin,
                               gmax = gmax,
                               run = parallelFAOutput,
                               parallel = FALSE)
            
            aic <- AICFunction(ll = ll,
                               k = k,
                               run = parallelFAOutput,
                               gmin = gmin,
                               gmax = gmax,
                               parallel = FALSE)
            
            aic3 <- AIC3Function(ll = ll,
                                 k = k,
                                 run = parallelFAOutput,
                                 gmin = gmin,
                                 gmax = gmax,
                                 parallel = FALSE)
        }
    }
    
    final <- base::proc.time() - ptm
    
    RESULTS <- list(dataset = dataset,
                    nUnits = n,
                    nVariables = p,
                    nOccassions = r,
                    normFactors = normFactors,
                    gmin = gmin,
                    gmax = gmax,
                    initalizationMethod = initMethod,
                    allResults = parallelFAOutput,
                    loglikelihood = ll,
                    nParameters = k,
                    trueLabels = membership,
                    ICLAll = icl,
                    BICAll = bic,
                    AICAll = aic,
                    AIC3All = aic3,
                    totalTime = final)
    
    class(RESULTS) <- "mvplnVGA"
    return(RESULTS)
    
}

initializationRunVGA <- function(G,
                                 dataset,
                                 TwoDdataset,
                                 r, p, d, n,
                                 normFactors,
                                 initMethod) {
    
    # arranging normalization factors
    libMat <- matrix(normFactors, n, d, byrow = T)
    libMatList <- list()
    for (i in 1:n) {
        libMatList[[i]] <- t(matrix(libMat[i, ], nrow = p))
    }
    
    # Initialization
    mu <- omega <- phi <- list() # mu is M;
    delta <- kappa <- sigma <- isigma <-
        iphi <- iomega <- list()
    # delta  is variational parameter Delta
    # kappa is variational parameter kappa
    # sigma is Psi in the math - kronecker of omega and Phi
    # isigma is inverse of Psi
    # iphi is inverse of Phi
    # iomega is inverse of omega
    
    m <- S <- list()
    # m is vectorized xi
    # S is Psi
    
    # Other intermediate items initialized
    Sk <- array(0, c(d, d, G) )
    start <- GX <- dGX <- zS <- list()
    
    iKappa <- iDelta <- startList <- list()
    # iKappa is inverse of kappa
    # iDelta is inverse of Delta
    
    
    if (initMethod == "kmeans") {
        # cat("\n initMethod == kmeans \n")
        zValue <- mclust::unmap(stats::kmeans(x = log(TwoDdataset + 1 / 3),
                                              centers = G)$cluster)
        
    } else if (initMethod == "random") {
        # cat("\n initMethod == random \n")
        if(G == 1) { # generating z if g = 1
            zValue <- as.matrix(rep.int(1, times = n),
                                ncol = G,
                                nrow = n)
        } else { # generating z if g>1
            zConv = 0
            while(! zConv) { # ensure that dimension of z is same as G (i.e.
                # if one column contains all 0s, then generate z again)
                zValue <- t(stats::rmultinom(n = n,
                                             size = 1,
                                             prob = rep(1 / G, G)))
                if(length(which(colSums(zValue) > 0)) == G) {
                    zConv <- 1
                }
            }
        }
    } else if (initMethod == "medoids") {
        # cat("\n initMethod == medoids \n")
        zValue <- mclust::unmap(cluster::pam(x = log(TwoDdataset + 1 / 3),
                                             k = G)$cluster)
    } else if (initMethod == "clara") {
        # cat("\n initMethod == clara \n")
        zValue <- mclust::unmap(cluster::clara(x = log(TwoDdataset + 1 / 3),
                                               k = G)$cluster)
    } else if (initMethod == "fanny") {
        # cat("\n initMethod == fanny \n")
        zValue <- mclust::unmap(cluster::fanny(x = log(TwoDdataset + 1 / 3),
                                               k = G)$cluster)
    }
    
    piG <- colSums(zValue) / n
    
    for (g in 1:G) {
        obs <- which(zValue[, g] == 1)
        mu[[g]] <- colMeans(log(TwoDdataset[obs, ] + 1 / 6))
        sigma[[g]] <- var(log(TwoDdataset[obs, ] + 1 / 6))
        isigma[[g]] <- solve(sigma[[g]])
        phi[[g]] <- diag(r) * sqrt(min(diag(var(log(TwoDdataset[obs, ] + 1 / 6)))))
        omega[[g]] <- diag(p) * sqrt(min(diag(var(log(TwoDdataset[obs, ] + 1 / 6)))))
        iphi[[g]] <- solve(phi[[g]])
        iomega[[g]] <- solve(omega[[g]])
    }
    
    for (g in 1:G) {
        start[[g]] <- log(TwoDdataset + 1/6) ###Starting value for M
        m[[g]] <- log(TwoDdataset + 1/6)
        S[[g]] <- list()
        delta[[g]] <- list()
        kappa[[g]] <- list()
        startList[[g]] <- list()
        for (i in 1:n) {
            startList[[g]][[i]] <- log(dataset[[i]])
            delta[[g]][[i]] <- diag(r) * 0.001
            kappa[[g]][[i]] <- diag(p) * 0.001
            S[[g]][[i]] <- delta[[g]][[i]] %x% kappa[[g]][[i]]
        }
    }
    
    
    it <- 1
    aloglik <- loglik <- NULL
    checks <- aloglik[c(1:5)] <- 0
    itMax <- 200
    
    while (checks == 0) {
        
        for (g in 1:G) {
            GX[[g]] <- dGX[[g]] <- zS[[g]] <- list()
            iDelta[[g]] <- iKappa[[g]] <- list()
            deltaO <- delta[[g]]
            kappaO <- kappa[[g]]
            
            for (i in 1:n) {
                iDelta[[g]][[i]] <- diag(c(t(diag(kappa[[g]][[i]])) %*%
                                               t(exp(log(libMatList[[i]]) +
                                                         startList[[g]][[i]] +
                                                         0.5 * diag(delta[[g]][[i]]) %*%
                                                         t(diag(kappa[[g]][[i]])))))) +
                    iphi[[g]] * sum(diag(iomega[[g]] %*%
                                             kappa[[g]][[i]]))
                delta[[g]][[i]] <- p * solve(iDelta[[g]][[i]])
                
                
                iKappa[[g]][[i]] <- diag(c(t(diag(delta[[g]][[i]])) %*%
                                               (exp(log(libMatList[[i]]) +
                                                        startList[[g]][[i]] +
                                                        t(0.5 * diag(kappa[[g]][[i]]) %*%
                                                              t(diag(delta[[g]][[i]]))))))) +
                    iomega[[g]] * sum(diag(iphi[[g]] %*%
                                               delta[[g]][[i]]))
                
                kappa[[g]][[i]] <- solve(iKappa[[g]][[i]])
                
                
                S[[g]][[i]] <- delta[[g]][[i]] %x% kappa[[g]][[i]]
                zS[[g]][[i]] <- zValue[i, g] * S[[g]][[i]]
                GX[[g]][[i]] <- TwoDdataset[i, ] -
                    exp(start[[g]][i, ] +
                            log(libMat[i,]) +
                            0.5 * diag(S[[g]][[i]])) -
                    (isigma[[g]]) %*% (start[[g]][i, ] - mu[[g]])
                m[[g]][i, ] <- start[[g]][i, ] + S[[g]][[i]] %*% GX[[g]][[i]]
                startList[[g]][[i]] <- t(matrix(m[[g]][i, ], nrow = p))
            }
            start[[g]] <- m[[g]]
            
            mu[[g]] <- colSums(zValue[, g] * m[[g]]) / sum(zValue[, g]) # this is xi
            
            # Updating Sample covariance
            muMat <- t(matrix(mu[[g]], nrow = p))
            phiM <- list()
            for (i in 1:n) {
                phiM[[i]] <- zValue[i,g]*(startList[[g]][[i]] - muMat) %*%
                    iomega[[g]] %*% t(startList[[g]][[i]] - muMat) +
                    zValue[i, g] * delta[[g]][[i]] * sum(diag(iomega[[g]] %*%
                                                                  kappa[[g]][[i]]))
            }
            phi[[g]] <- Reduce("+", phiM) / sum(zValue[, g] * p)
            iphi[[g]] <- solve(phi[[g]])
            
            omegaM <- list()
            for (i in 1:n) {
                omegaM[[i]] <- zValue[i, g] * t(startList[[g]][[i]] - muMat) %*%
                    iphi[[g]] %*% (startList[[g]][[i]] - muMat) +
                    zValue[i, g] * kappa[[g]][[i]] * sum(diag(iphi[[g]] %*%
                                                                  delta[[g]][[i]]))
            }
            omega[[g]] <- Reduce("+", omegaM) / sum(zValue[, g] * r)
            iomega[[g]] <- solve(omega[[g]])
            sigma[[g]] <- phi[[g]] %x% omega[[g]]
            isigma[[g]] <- iphi[[g]] %x% iomega[[g]]
        }
        
        piG <- colSums(zValue) / n
        # Internal functions
        funFive <- function(x, y = isigma[[g]]) {
            sum(diag(x %*% y))
        }
        
        FMatrix <- matrix(NA, ncol = G, nrow = n)
        
        for (g in 1:G) {
            two <- rowSums(exp(m[[g]] + log(libMat) +
                                   0.5 * matrix(unlist(lapply(S[[g]], diag)),
                                                ncol = d, byrow = TRUE)))
            five <- 0.5 * unlist(lapply(S[[g]], funFive))
            six <- 0.5 * log(unlist(lapply(S[[g]], det)))
            FMatrix[, g] <- piG[g] * exp(rowSums(m[[g]] * TwoDdataset) -
                                             two - rowSums(lfactorial(TwoDdataset)) +
                                             rowSums(log(libMat) * TwoDdataset) - 0.5 *
                                             mahalanobis(m[[g]], center = mu[[g]], cov = isigma[[g]],
                                                         inverted = TRUE) - five + six + 0.5 *
                                             log(det(isigma[[g]])) - d/2)
        }
        
        # loglik[it] <- sum(log(rowSums(FMatrix)))
        # if FMatrix values is 0, then log(0) = -Inf and logL error
        rowSumsFMatrix <- rowSums(FMatrix)
        if(any(rowSumsFMatrix == 0) == TRUE) {
            zeroRows <- which(rowSumsFMatrix == 0)
            rowSumsFMatrix[zeroRows] <- min(rowSumsFMatrix[-zeroRows])
        }
        loglik[it] <- sum(log(rowSumsFMatrix))
        
        
        zValue <- FMatrix / rowSumsFMatrix
        if (it <= 5) {
            zValue[zValue == "NaN"] <- 0
        }
        
        # cat("\n Initialization loglik[it - 1] - loglik[it - 2]:", loglik[it - 1] - loglik[it - 2])
        if (it > 5) {
            # Aitkaine's stopping criterion
            if ((loglik[it - 1] - loglik[it - 2]) == 0) checks <- 1 else {
                a <- (loglik[it] - loglik[it - 1]) / (loglik[it - 1] - loglik[it - 2])
                addTo <- (1 / (1 - a) * (loglik[it] - loglik[it - 1]))
                aloglik[it] <- loglik[it - 1] + addTo
                if (abs(aloglik[it] - aloglik[it - 1]) < 0.05) {
                    checks <- 1
                } else {
                    checks <- checks
                }
            }
        }
        
        it <- it + 1
        if (it == itMax) {
            checks <- 1
        }
        finalPhi <- finalOmega <- list()
        for (g in 1:G) {
            finalPhi[[g]] <- phi[[g]] / diag(phi[[g]])[1]
            finalOmega[[g]] <- omega[[g]] * diag(phi[[g]])[1]
        }
    }
    
    programclust <- mclust::map(zValue)
    
    initRunOutput <- list(mu = mu,
                          sigma = sigma,
                          omega = omega,
                          phi = phi,
                          omega = omega,
                          delta = delta,
                          kappa = kappa,
                          isigma = isigma,
                          iphi = iphi,
                          iomega = iomega,
                          m = m,
                          S = S,
                          Sk = Sk,
                          start = start,
                          GX = GX,
                          dGX = dGX,
                          zS = zS,
                          iKappa = iKappa,
                          iDelta = iDelta,
                          startList = startList,
                          zValue = zValue,
                          loglik = loglik,
                          finalLogLik = loglik[it-1],
                          piG = piG,
                          clusterlabels = programclust,
                          iterations = it)
    
    class(initRunOutput) <- "initializationRunVGA"
    
    return(initRunOutput)
}


# [END]
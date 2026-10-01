trueG <- 2 # due cluster
truer <- 2 # due tempi
truep <- 3 # tre variabili
trueN <- 1000 # simulo 1000 unità

# Mu: matrice delle medie
trueM1 <- matrix(rep(6.20, times = 6),
                 ncol = truep,
                 nrow = truer,
                 byrow = TRUE); trueM1

trueM2 <- matrix(rep(1.50, times = 6),
                 ncol = truep,
                 nrow = truer,
                 byrow = TRUE); trueM2

trueMall <- rbind(trueM1, trueM2); trueMall

# Phi: matrice di covarianza per i tempi 
truePhi1 <- diag(c(1, 1)); truePhi1
truePhi2 <- diag(c(1, 0.70)); truePhi2

truePhiall <- rbind(truePhi1, truePhi2); truePhiall

# Omega
trueOmega1 <- diag(c(1.67, 1.46, 1.44)); trueOmega1
trueOmega2 <- diag(c(0.75, 0.82, 0.90)); trueOmega2

trueOmegaAll <- rbind(trueOmega1, trueOmega2)

sampleData1 <- mixMVPLN::mvplnDataGenerator(nOccasions = truer,
                                            nResponses = truep,
                                            nUnits = trueN,
                                            mixingProportions = c(0.6, 0.4),
                                            matrixMean = trueMall,
                                            phi = truePhiall,
                                            omega = trueOmegaAll)

sampleData1$dataset

# TMM (Trimmed Mean of M-values)

Risultati1 <- mixMVPLN::mvplnVGAclus(
  dataset = sampleData1$dataset,
  membership = sampleData1$truemembership,
  gmin = 1,
  gmax = 5,
  initMethod = "clara",
  nInitIterations = 3,
  normalize = "Yes")

str(Risultati1)

Risultati1$BICAll
Risultati1$ICLAll
Risultati1$AICAll
Risultati1$AIC3All

# grafici
par(mfrow=c(1,2))

graphics::matplot(Risultati1$loglikelihood, xlab = "Run",
                  ylab = "logL", type = c("b"), pch = 1, lty = 2, xaxt="n") 
graphics::axis(1, at = 1:5, las = 1)

ICvalues <- matrix(c(Risultati1$BICAll$allBICvalues,   
                     Risultati1$ICLAll$allICLvalues,
                     Risultati1$AICAll$allAICvalues,
                     Risultati1$AIC3All$allAIC3values),
                   ncol=4) 
graphics::matplot(ICvalues, xlab = "Run", ylab = "Information criteria value", 
                  type = c("b"), pch = 1, col = 1:4, xaxt="n", lwd = 3) 

graphics::axis(1, at = 1:5, las = 1)

graphics::legend("topright", 
                 legend = c("BIC", "ICL", "AIC", "AIC3"), 
                 col = 1:4, pch = 1, bty = "n")

mat <- Risultati1$allResults[[2]]$probaPost
mat[1,1] <- 0.3
mat[1,2] <- 0.7
mvplnClustVisuals <- mixMVPLN::mvplnVisualize(
  dataset = sampleData1$dataset,
  plots = 'bar',
  probabilities = mat,
  clusterMembershipVector = Risultati1$allResults[[2]]$clusterlabels)

mvplnClustVisuals

Risultati1MCMC <- mixMVPLN::mvplnMCMCclus(
  dataset = sampleData1$dataset,
  membership = sampleData1$truemembership,
  nIterations = 40,
  gmin = 1,
  gmax = 5,
  initMethod = "clara",
  nInitIterations = 3,
  normalize = "Yes")

Risultati1HYB <- mixMVPLN::mvplnHybriDclus(
  dataset = sampleData1$dataset,
  membership = sampleData1$truemembership,
  gmin = 1,
  gmax = 5,
  initMethod = "clara",
  nInitIterations = 2,
  normalize = "Yes")
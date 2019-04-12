### R code from vignette source 'CloneSeeker.Rnw'

###################################################
### code chunk number 1: lib
###################################################
library(CloneSeeker)


###################################################
### code chunk number 2: simTumor
###################################################
set.seed(21303) # for reproducibility
simTumor <- Tumor(c(5, 3, 2), rounds = 100, 
                  nu = 10, pcnv = 0.8, norm.contam = FALSE)


###################################################
### code chunk number 3: truePsi
###################################################
simTumor@psi


###################################################
### code chunk number 4: cloneMeta
###################################################
class(simTumor@clones)
length(simTumor@clones)


###################################################
### code chunk number 5: oneClone
###################################################
oneClone <- simTumor@clones[[1]]
class(oneClone)
length(oneClone)
names(oneClone)


###################################################
### code chunk number 6: oneCN
###################################################
dim(oneClone$cn)
summary(oneClone$cn)


###################################################
### code chunk number 7: oneMut
###################################################
dim(oneClone$seq)
oneClone$seq


###################################################
### code chunk number 8: simData
###################################################
simData <- generateTumorData(simTumor,
                             snps.seq = 10000,
                             snps.cgh = 600000,
                             mu = 70,
                             sigma.reads = 25,
                             sigma0.lrr = 0.15,
                             sigma0.baf = 0.03,
                             density.sigma = 0.1)


###################################################
### code chunk number 9: peekData
###################################################
class(simData)
length(simData)
names(simData)


###################################################
### code chunk number 10: simCNData
###################################################
cnDat <- simData$cn.data
dim(cnDat)
summary(cnDat)


###################################################
### code chunk number 11: simMutData
###################################################
dim(simData$seq.data)
seqDat <- simData$seq.data
somatic <- seqDat[seqDat$status=='somatic',]
dim(seqDat)
summary(seqDat)
table(seqDat$status)


###################################################
### code chunk number 12: psis
###################################################
psis <- generateSimplex(20, 5)
dim(psis)
head(psis)
tail(psis)


###################################################
### code chunk number 13: cnmodels
###################################################
cnmodels <- expand.grid(rep(list(0:5),5))
include <- sapply(1:nrow(cnmodels), function(i) {
  length(which(cnmodels[i,] >= 1))==5 | length(which(cnmodels[i,] <= 1)) == 5
})
cnmodels <- cnmodels[include,]


###################################################
### code chunk number 14: pars
###################################################
pars <- list(sigma0 = 1,    # SNP-wise standard deviation
             ktheta = 0.3,  # geometric prior parameter on number of clones
             theta = 0.9,   # geometric prior parameter on copy number changes
             mtheta = 0.9,  # gemoetric prior parameter on point mutations
             alpha = 0.5,   # parameter for a symmetric Dirichlet prior on psi
             thresh = 0.04, # smallest possible detectble clone
             cutoff = 100,  # filter out copy number segments supported by fewer SNPs
             Q = 100,       # number of new psi vectors resamples at each iteration
             iters = 4)     # number of iterations


###################################################
### code chunk number 15: realWork
###################################################
f <- system.file("auxiliary/stash.Rda", package="CloneSeeker")
if (file.exists(f)) {
  load(f)
} else {
  resCN <- seekClones(cndata = cnDat, vardata = NULL, cnmodels = cnmodels, psis, pars = pars)
  resMut <- seekClones(cndata = NULL, vardata = seqDat, cnmodels = cnmodels, psiset = psis, pars = pars)
  resBoth <- seekClones(cndata = cnDat, vardata = somatic, cnmodels = cnmodels, psis, pars = pars)
  save(resCN, resMut, resBoth, file = "../inst/auxiliary/stash.Rda")
}
rm(f)


###################################################
### code chunk number 16: resCN (eval = FALSE)
###################################################
## resCN <- seekClones(cndata = cnDat, vardata = NULL, 
##                     cnmodels = cnmodels, psiset = psis, pars = pars)


###################################################
### code chunk number 17: compare.psis
###################################################
resCN$psi
simTumor@psi


###################################################
### code chunk number 18: CNassign
###################################################
trueCN_Assignments <- t(sapply(1:nrow(resCN$filtered.data$cndata.filt),
function(i) {
  index <- rownames(simTumor@clones[[1]]$cn) == 
           rownames(resCN$filtered.data$cndata.filt)[i]
  sapply(1:length(simTumor@clones),function(j){
    simTumor@clones[[j]]$cn$A[index] + simTumor@clones[[j]]$cn$B[index]
  })
}))

inferredCN_Assignments <- (resCN$A+resCN$B)[,1:length(simTumor@clones)]
colnames(inferredCN_Assignments) <- colnames(trueCN_Assignments) <- 
  paste("C", 1:3)

data.frame(Truth = trueCN_Assignments,
           Infer = inferredCN_Assignments)


###################################################
### code chunk number 19: resMut (eval = FALSE)
###################################################
## resMut <- seekClones(cndata = NULL, vardata = seqDat, 
##                      cnmodels = cnmodels, psiset = psis, pars = pars)


###################################################
### code chunk number 20: resMutResults
###################################################
resMut$psi
simTumor@psi


###################################################
### code chunk number 21: trueMut (eval = FALSE)
###################################################
## trueMutAssignments <- sapply(1:length(simTumor@clones),function(i){
##   sapply(1:length(resMut$filtered.data$mutdata.filt$mut.id),function(j){
##     index <- which(simTumor@clones[[i]]$seq$mut.id==resMut$filtered.data$mutdata.filt$mut.id[j])
##     if(length(index)==0){
##       val <- 0
##     }else{
##       val <- simTumor@clones[[i]]$seq$mutated.copies[index]
##     }
##     val
##   })
## })
## inferredMutAssignments <- resMut$mutated[,1:length(simTumor@clones)]
## colnames(inferredMutAssignments) <- colnames(trueMutAssignments) <- c('C1','C2','C3')
## 
## data.frame(Truth = trueMutAssignments,
##            Infer = inferredMutAssignments)


###################################################
### code chunk number 22: resBoth (eval = FALSE)
###################################################
## resBoth <- seekClones(cndata = cnDat, vardata = somatic, 
##                       cnmodels = cnmodels, psiset = psis, pars = pars)


###################################################
### code chunk number 23: reBothResults
###################################################
resBoth$psi
simTumor@psi



hardthres = function(v, low=0.9, high=1.1){
  n = length(v)
  for (i in 1:n){ if (v[i]>low && v[i]<high) v[i] = 1 }
  v
}

falcon.getASCN.epsilon = function (readMatrix, rdep = NULL, tauhat = NULL, threshold = 0.15, 
          pOri = c(0.49, 0.51), error = 1e-05, maxIter = 1000) 
{
  AN = readMatrix$AN
  BN = readMatrix$BN
  AT = readMatrix$AT
  BT = readMatrix$BT
  if (is.null(rdep)) 
    rdep = median(AT + BT)/median(AN + BN)
  if (is.null(tauhat)) 
    tauhat = getChangepoints(readMatrix, pOri = pOri, error = error, 
                             maxIter = maxIter)
  N = length(AT)
  tau = sort(unique(c(1, tauhat, N)))
  K = length(tau) - 1
  pa = pb = rep(0, K)
  for (i in 1:K) {
    ids = tau[i]:(tau[i + 1] - 1)
    if (i == K) 
      ids = tau[i]:tau[i + 1]
    p = as.numeric(.Call("GetP", as.numeric(AT[ids]), as.numeric(BT[ids]), 
                         as.numeric(AN[ids]), as.numeric(BN[ids]), as.numeric(error), 
                         as.numeric(maxIter), as.numeric(pOri), PACKAGE = "falcon"))
    pa[i] = p[1]
    pb[i] = p[2]
    if (diff(p) < 0.1) {
      temp = as.numeric(.Call("LikH", as.numeric(AT[ids]), 
                              as.numeric(BT[ids]), as.numeric(AN[ids]), as.numeric(BN[ids]), 
                              as.numeric(p), PACKAGE = "falcon"))
      p2 = sum(AT[ids] + BT[ids])/sum(AT[ids] + BT[ids] + 
                                        AN[ids] + BN[ids])
      temp2 = as.numeric(.Call("Lik", as.numeric(AT[ids]), 
                               as.numeric(BT[ids]), as.numeric(AN[ids]), as.numeric(BN[ids]), 
                               as.numeric(rep(p2, 2)), PACKAGE = "falcon"))
      if (!is.na(temp)[1] && !is.na(temp[2]) && !is.na(temp2)) {
        bic = temp[1] - temp2 - temp[2]/2 + log(p2 * 
                                                  (1 - p2) * sum(AT[ids] + BT[ids] + AN[ids] + 
                                                                   BN[ids]))/2 + log(2 * pi)/2
      }
      else if (!is.na(temp)[1] && !is.na(temp2)) {
        bic = temp[1] - temp2 + log(p2 * (1 - p2) * 
                                      sum(AT[ids] + BT[ids] + AN[ids] + BN[ids]))/2 + 
          log(2 * pi)/2
      }
      if (bic < 0) {
        pa[i] = p2
        pb[i] = p2
      }
    }
  }
  rawcns1 = pa/(1 - pa)/rdep
  rawcns2 = pb/(1 - pb)/rdep
  cns1 = hardthres(rawcns1, low = 1 - threshold, high = 1 + 
                     threshold)
  cns2 = hardthres(rawcns2, low = 1 - threshold, high = 1 + 
                     threshold)
  Haplotype = list()
  for (i in 1:K) {
    # if (cns1[i] != cns2[i]) {
      ids = tau[i]:(tau[i + 1] - 1)
      if (i == K) 
        ids = tau[i]:tau[i + 1]
      gt = 1/(1 + (pb[i]/pa[i])^(AT[ids] - BT[ids]) * 
                ((1 - pb[i])/(1 - pa[i]))^(AN[ids] - BN[ids]))
      temp3 = rep("A", length(ids))
      temp3[which(gt > 0.5)] = "B"
      Haplotype[[i]] = temp3
    # }
  }
  return(list(tauhat = tauhat, ascn = rbind(cns1, cns2), Haplotype = Haplotype, 
              readMatrix = readMatrix))
}

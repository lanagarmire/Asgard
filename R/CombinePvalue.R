
##Combine P-value
CombineP = function (p){
  keep <- (p > 0) & (p <= 1)
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    warning("Must have at least two valid p values")
    res <- list(chisq = NA_real_, df = NA_integer_, p = NA_real_, 
      validp = p[keep])
  }
  else {
    lnp <- log(p[keep])
    chisq <- (-2) * sum(lnp)
    df <- 2 * length(lnp)
    if (length(lnp) != length(p)) {
      warning("Some studies omitted")
    }
    res <- pchisq(chisq,df, lower.tail = FALSE)
  }
  return(res)
}


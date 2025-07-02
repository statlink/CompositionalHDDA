cv.alfahdda <- function(ina, x, a = seq(-1, 1, by = 0.1), d_select = "both", 
               threshold = c(0.001, 0.005, 0.05, 1:9 * 0.1), folds = NULL, stratified = TRUE, 
               nfolds = 10, seed = NULL) {

  if ( min(x) == 0 )  a <- a[a > 0]
  la <- length(a)   ;   lt <- length(threshold)
  n <- dim(x)[1]

  if ( is.null(folds) )
  folds <- Compositional::makefolds(1:n, nfolds = nfolds, stratified = stratified, seed = seed)
  nfolds <- length(folds)

  if  ( d_select == "Cattell"  ) {
    config <- as.matrix( expand.grid(threshold = threshold, a = a) )
    p <- dim(config)[1]
    cv <- matrix(nrow = nfolds, ncol = p) 
  } else if ( d_select == "BIC" )  {
    cv <- matrix(nrow = nfolds, ncol = la)
  } else if ( d_select == "both" ) {
    config <- as.matrix( expand.grid(threshold = threshold, a = a) )
    p <- dim(config)[1]
    cv <- matrix(nrow = nfolds, ncol = p + la) 
  }
  
  for ( k in 1:nfolds ) {
    inatrain <- ina[ -folds[[ k ]] ]
    inatest <- ina[ folds[[ k ]]  ]
    xtrain <- x[-folds[[ k ]], ]
    xtest <- x[folds[[ k ]], ]
    
    for ( i in 1:la ) {
      ytrain <- Compositional::alfa(xtrain, a[i])$aff  ## apply the alpha-transformation
      ytest <- Compositional::alfa(xtest, a[i])$aff
      if ( d_select == "Cattell" ) {
        for ( j in 1:length(threshold) ) {
          mod <- HDclassif::hdda( data = ytrain, cls = inatrain, model = "all", d_select = "Cattell", 
                                  threshold = threshold[j] )
          est <- predict(mod, ytest)$class
          cv[k, (i - 1) * lt + j] <- mean(est == inatest)
        } 
      } else if ( d_select == "BIC" ) {
        mod <- HDclassif::hdda( data = ytrain, cls = inatrain, model = "all", d_select = "Cattell", 
                                threshold = threshold[j] )
        est <- predict(mod, ytest)$class
        cv[k, i] <- mean(est == inatest)
      } else if ( d_select == "both" ) {
        for ( j in 1:length(threshold) ) {
          mod <- HDclassif::hdda( data = ytrain, cls = inatrain, model = "all", d_select = "Cattell", 
                                  threshold = threshold[j] )
          est <- predict(mod, ytest)$class
          cv[k, (i - 1) * lt + j] <- mean(est == inatest)
        } 
        mod <- HDclassif::hdda( data = ytrain, cls = inatrain, model = "all", d_select = "BIC", 
                                threshold = threshold[j] )
        est <- predict(mod, ytest)$class
        cv[k, p + i] <- mean(est == inatest)
      }
      
    }  ##  end for ( i in 1:length(a) ) {
  }  ##  end  for ( k in 1:nfolds ) {
  
  if  ( d_select == "Cattell"  ) {
    res <- cbind(config, Rfast::colmeans(cv) )
    colnames(res)[3] <- "performance"
  } else if ( d_select == "BIC" )  {
    res <- cbind(a, Rfast::colmeans(cv) )
    colnames(res)[2] <- "performance"
  } else if ( d_select == "both" ) {
    config <- rbind(config, cbind(NA, a) )
    res <- cbind(config, Rfast::colmeans(cv) )
    colnames(res)[3] <- "performance"
  }
  
  res
}  

  
  
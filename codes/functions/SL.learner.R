SL.dbarts = function(Y, X, newX=NULL, family, ...) {
  require("dbarts")
  model = dbarts::bart(x.train = X, y.train = Y, x.test = newX,keeptrees = TRUE, verbose = F)
  if (is.null(newX)) {
    if (family$family == "gaussian") {
      pred = model$yhat.train.mean}  
    if (family$family == "binomial") { # probit
      pred = colMeans(pnorm(model$yhat.train))}
  }else{
    if (family$family == "gaussian") {
      pred = model$yhat.test.mean}  
    if (family$family == "binomial") {
      pred = colMeans(pnorm(model$yhat.test))}
  }
  fit = list(object = model)
  class(fit) = c("SL.dbarts")
  out = list(pred = pred, fit = fit)
  return(out)
}

predict.SL.dbarts = function(object, newdata, ...) {
  require("dbarts")
  pred =colMeans(predict(object$object, newdata = newdata))
  return(pred)
}



SL.ranger.1000 <- function(...,min.node.size=1,num.trees = 1000){
  SL.ranger(...,min.node.size=min.node.size,num.trees=num.trees)
}
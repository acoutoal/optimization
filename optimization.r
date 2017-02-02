#############################################################################
#
# Steepest descent, minimizes a function given data and the function's gradient
#
#############################################################################
steepest_descent <- function(fun=NA, grad=NA, data=NA,start=NA,eta=0.1, epsilon=1E-6, nmax=1000){
  
  converged  = F
  param_old  = start
  epsilon_hat= NA
  param_hat  = NA
  param_new  = NA
  G_hat      = NA
  L_hat      = NA
  param_hst  = matrix(nrow = nmax,ncol = length(start));
  m          = length(start)
  I          = diag(nrow=m,ncol=m)
  
  
  for( idx_iter in 1:nmax){
    
    G = grad(param_old,data)

    param_new = param_old - eta * I %*% G
    
    param_hst[idx_iter,] = param_new
    
    #Stopping criteria
    if( all(abs((param_new - param_old)/param_old)<epsilon) ){
      epsilon_hat = t(abs((param_new - param_old)/param_old))
      param_hat   = param_new
      G_hat = grad(param_hat,data)
      L_hat = fun(params = param_hat,data=data)
      converged=T
      break
    }
    
    param_old = param_new
    
  }
  return(list(g_hat=G_hat,l_hat=L_hat,conv=converged,iter=idx_iter, param_hat=param_hat, param_trace=param_hst, epsilon=epsilon_hat))
} # END NEWTON - RAPHSON


#############################################################################
#
# Newton raphson, minimizes a function given data, the function's gradient and hessian
#
#############################################################################
newton_raphson <- function(fun=NA, grad=NA,hess=NA,data=NA,start=NA,eta=0.1, epsilon=1E-6, nmax=1000){
  
  converged  = F
  param_old  = start
  epsilon_hat= NA
  param_hat  = NA
  param_new  = NA
  G_hat      = NA
  H_hat      = NA
  L_hat      = NA
  param_hst  = matrix(nrow = nmax,ncol = length(start));
  
  
  
  for( idx_iter in 1:nmax){
    
    G = grad(param_old,data)
    H = hess(param_old,data)
    
    param_new = param_old - eta * solve(H) %*% G
    
    param_hst[idx_iter,] = param_new
 
    #Stopping criteria
    if( all(abs((param_new - param_old)/param_old)<epsilon) ){
      epsilon_hat = t(abs((param_new - param_old)/param_old))
      param_hat   = param_new
      G_hat = grad(param_hat,data)
      H_hat = hess(param_hat,data)
      L_hat = fun(params = param_hat,data=data)
      converged=T
      break
    }
    
    param_old = param_new
    
  }
  return(list(g_hat=G_hat,h_hat=H_hat,l_hat=L_hat,conv=converged,iter=idx_iter, param_hat=param_hat, param_trace=param_hst, epsilon=epsilon_hat))
} # END NEWTON - RAPHSON

#############################################################################
#
# Gradient descent, minimizes a function given data and the function's gradient
#
#############################################################################
gradient_descent <- function(fun=NA, grad=NA, data=NA, start=NA, eta=0.1, epsilon=1E-6, nmax=1000){
  
  converged  = F
  param_old  = start
  epsilon_hat= NA
  param_hat  = NA
  param_new  = NA
  G_hat      = NA
  H_hat      = NA
  L_hat      = NA
  param_hst  = matrix(nrow = nmax,ncol = length(start));
  
  
  
  for( idx_iter in 1:nmax){
    
    G = grad(param_old,data)

    param_new = param_old - eta * G
    
    param_hst[idx_iter,] = param_new
  
    #Stopping criteria
    if( all(abs((param_new - param_old)/param_old)<epsilon) ){
      
      epsilon_hat = t(abs((param_new - param_old)/param_old))
      param_hat   = param_new
      G_hat       = grad(param_hat,data)
      L_hat       = fun(params = param_hat,data=data)
      converged   = T
      break
      
    }
    
    param_old = param_new
    
  }
  return(list(g_hat=G_hat,l_hat=L_hat,conv=converged,iter=idx_iter, param_hat=param_hat, param_trace=param_hst, epsilon=epsilon_hat))
} # END NEWTON - RAPHSON


#############################################################################
#
# Levenber marquardt least squares, minimizes a error function given data, a function and it's jacobian
#
#############################################################################
levenberg_marquardt <-function (fun,jacob,start,data,lambda=0.05,nmax=1000){
  converged  = F
  param_old  = start
  epsilon_hat= NA
  param_hat  = NA
  param_new  = NA
  G_hat      = NA
  H_hat      = NA
  L_hat      = NA
  param_hst  = matrix(nrow = nmax,ncol = length(start));
  m          = length(param_old)
  param_old = start
  
  for( idx_iter in 1:nmax){
    
    J            = jacob(param_old,data)
    y_hat        = fun(param_old,data)
    I            = diag(m)        #Levenberg regularization
    diag(I)      = diag(t(J)%*%J) #Marquardt regularization 
    
    delta = solve( t(J)%*%J + lambda * I, t(J) %*% (data$y - y_hat) )
    param_new = param_old + delta
    
    param_hst[idx_iter,] = param_new
    
    #Stopping criteria
    if( all(abs((param_new - param_old)/param_old)<epsilon) ){
      
      epsilon_hat = t(abs((param_new - param_old)/param_old))
      param_hat   = param_new
      G_hat       = jacob(param_hat,data)
      converged   = T
      break
      
    }

    param_old=param_new
    
  }
  
  return(list(g_hat=G_hat,conv=converged,iter=idx_iter, param_hat=param_hat, param_trace=param_hst, epsilon=epsilon_hat))
  
}

#############################################################################
#
# Gauss newton least squares, minimizes a error function given data, a function and it's jacobian
#
#############################################################################
gauss_newton <-function (fun,jacob,start,data,lambda=1,nmax=1000){
  converged  = F
  param_old  = start
  epsilon_hat= NA
  param_hat  = NA
  param_new  = NA
  G_hat      = NA
  H_hat      = NA
  L_hat      = NA
  param_hst  = matrix(nrow = nmax,ncol = length(start));
  m          = length(param_old)
  param_old = start
  
  for( idx_iter in 1:nmax){
    
    J            = jacob(param_old,data)
    y_hat        = fun(param_old,data)

    delta = solve( t(J)%*%J , t(J) %*% (data$y - y_hat) )
    param_new = param_old + lambda * delta
    
    param_hst[idx_iter,] = param_new
    
    #Stopping criteria
    if( all(abs((param_new - param_old)/param_old)<epsilon) ){
      
      epsilon_hat = t(abs((param_new - param_old)/param_old))
      param_hat   = param_new
      G_hat       = jacob(param_hat,data)
      converged   = T
      break
      
    }
    
    param_old=param_new
    
  }
  
  return(list(g_hat=G_hat,conv=converged,iter=idx_iter, param_hat=param_hat, param_trace=param_hst, epsilon=epsilon_hat))
  
}




linear_sqr_error <- function(params, data){
  e = t(data$x %*% params - data$y) %*% (data$x %*% params - data$y) / length(data$y)
  return(e)
}

gradient_line <- function(params, data){
  g = 2 * t(data$x) %*% (data$x %*% params - data$y) / length(data$y)
  return(g)
}

hessian_line <- function(params, data){
  h = 2 * t(data$x) %*% (data$x) / length(data$y)
  return(h)
}

jacobian_line <- function(params,data){
  return(data$x)
}

predict_line <- function(params, data){
  return( data$x %*% params)
}



#Unit tests for all functions
unit_test <- function(){
  epsilon = 1E-6
  nmax    = 1000
  
  ###########################################################################################################
  #Test case 1
  b = matrix(c(3, -3),nrow = 2)
  x = matrix(rnorm(2*nmax),ncol=2);
  y = x%*%b + rnorm(nmax)
  coef(lm(y ~ x))
  
  data = list(x=x,y=y)
  
  sol = gauss_newton       (fun = predict_line, jacob = jacobian_line, data=data, start=b*3)

  sol = levenberg_marquardt(fun = predict_line, jacob = jacobian_line, data=data, start=b*3)

  sol = gradient_descent   (fun = linear_sqr_error, grad = gradient_line, data = data, start = b*3)
  
  sol = newton_raphson     (fun = linear_sqr_error, grad = gradient_line, hess=hessian_line, data=data, start=b*3,eta = 0.5)

  sol = steepest_descent   (fun = linear_sqr_error, grad = gradient_line, data=data, start=b*3,eta = 0.5)
  
  
  print(t(sol$param_hat))
  print(t(b))
  
  stopifnot( sol$conv == T)                                   #convergence
  stopifnot( all(sol$g_hat<epsilon))                          #gradient should be close to zero
  stopifnot( all( abs(sol$param_hat - b) < 0.05*b ) ) #Less than 5% error from the original solution
  
  print("Unit test 1: PASSED")
  
  ###########################################################################################################
  #Test case 2
  a = 0.8
  b = 20
  p = rbeta(n = nmax,shape1 = a,shape2 = b)
  
  sol=newton_raphson(fun=loglik_beta, grad=gradient_beta,hess=hessian_beta,data=p,start=c(a+2,b+5),eta=0.05, epsilon=1E-12, nmax=10000)
  
  print(t(sol$param_hat))
  print(c(a,b))
  
  stopifnot( sol$conv == T)                                   #convergence
  stopifnot( all(sol$g_hat<epsilon))                          #gradient should be close to zero
  stopifnot( all(abs(sol$param_hat - c(a,b)) < 0.05*c(a,b)) ) #Less than 5% error from the original solution
  stopifnot( ( loglik_beta(params=c(a,b), data=p) - sum(log(dbeta(p,a,b))) ) < epsilon )
  
  print("Unit test 2: PASSED")
  
  ###########################################################################################################
  #Test case N
  #sol=gradient_descent(fun=loglik_beta, grad=gradient_beta,data=p,start=c(a+1,b+1),eta=0.05, epsilon=1E-12, nmax=10000)
  
}



adjust_proposal = function( current_sd,  accept_rate,  nprops,  target_ac){
  adj = (1.0/(1.0*nprops))*(accept_rate - target_ac);
  new_sd = current_sd*exp(adj);
  return(new_sd);
}

omega_sampler = function(delta, a, b){
  
  a_star = a + sum(delta)
  b_star = b + sum(delta) + length(delta)
  
  return(rbeta(1, a_star, b_star))
}


rPoisGamma = function(n, beta1, beta2, gamma1, gamma2, phi, burnin=1000){
  
  p=ncol( matrix(beta2, nrow=1))-1
  r=ncol( matrix(beta1,nrow=1)) -1
  
  s=ncol( matrix(gamma1, nrow=1))
  k = ncol( matrix(gamma2, nrow=1))-1
  
  u<-max(p,k,r,s)
  
  N = n + burnin
  Y1<- array(1,N)
  Y2<- array(1,N)
  
  mu1<-array(1,c(N,1))
  mu2<-array(1,c(N,1))
  
  for(t in (u+1):N){
    
    # Conditional means
    mu1[t]<-exp(beta1%*%c(1,log(Y1[(t-1):(t-r)]))+
                  gamma1%*%log(Y2[(t-1):(t-s)]+1) )
    
    Y1[t]<-rgamma(n=1,shape=1/phi,scale=mu1[t]*phi)
    
    
    
    mu2[t]<-exp(beta2%*%c(1,log(Y2[(t-1):(t-p)]+1)) +
                  gamma2%*%log(Y1[t:(t-k)]) )
    
    
    Y2[t]<-rpois(n=1,lambda=mu2[t])
  } 
  
  
  Y = cbind(Y1, Y2)
  Y = Y[(burnin+1):N, ]
  
  return(Y)
}


rGeoGamma = function(n, beta1, beta2, gamma1, gamma2, phi, burnin=1000){
  
  p=ncol( matrix(beta2, nrow=1))-1
  r=ncol( matrix(beta1,nrow=1)) -1
  
  s=ncol( matrix(gamma1, nrow=1))
  k = ncol( matrix(gamma2, nrow=1))-1
  
  u<-max(p,k,r,s)
  
  N = n + burnin
  Y1<- array(1,N)
  Y2<- array(1,N)
  
  mu1<-array(1,c(N,1))
  mu2<-array(1,c(N,1))
  
  for(t in (u+1):N){
    
    # Conditional means
    mu1[t]<-exp(beta1%*%c(1,log(Y1[(t-1):(t-r)]))+
                  gamma1%*%log(Y2[(t-1):(t-s)]+1) )
    
    Y1[t]<-rgamma(n=1,shape=1/phi,scale=mu1[t]*phi)
    
    
    
    mu2[t]<- beta2%*%c(1,log(Y2[(t-1):(t-p)]+1)) +
                  gamma2%*%log(Y1[t:(t-k)]) 
    
    
    Y2[t]<- rgeom(n=1, exp(mu2[t])/(1 + exp(mu2[t])))
  } 
  
  
  Y = cbind(Y1, Y2)
  Y = Y[(burnin+1):N, ]
  
  return(Y)
}

rPoisPois = function(n, beta1, beta2, gamma1, gamma2, burnin=1000){
  
  p=ncol( matrix(beta2, nrow=1))-1
  r=ncol( matrix(beta1,nrow=1)) -1
  
  s=ncol( matrix(gamma1, nrow=1))
  k = ncol( matrix(gamma2, nrow=1))-1
  
  u<-max(p,k,r,s)
  
  N = n + burnin
  Y1<- array(1,N)
  Y2<- array(1,N)
  
  mu1<-array(1,c(N,1))
  mu2<-array(1,c(N,1))
  
  for(t in (u+1):N){
    
    # Conditional means
    mu1[t]<-exp(beta1%*%c(1,log(Y1[(t-1):(t-r)]+1))+
                  gamma1%*%log(Y2[(t-1):(t-s)]+1) )
    
    Y1[t]<- rpois(n=1,lambda=mu1[t])
    
    
    
    mu2[t]<-exp(beta2%*%c(1,log(Y2[(t-1):(t-p)]+1)) +
                  gamma2%*%log(Y1[t:(t-k)]+1) )
    
    
    Y2[t]<-rpois(n=1,lambda=mu2[t])
  } 
  
  
  Y = cbind(Y1, Y2)
  Y = Y[(burnin+1):N, ]
  
  return(Y)
}

rGammaPois = function(n, beta1, beta2, gamma1, gamma2, phi, burnin=1000){
  
  p=ncol( matrix(beta2, nrow=1))-1
  r=ncol( matrix(beta1,nrow=1)) -1
  
  s=ncol( matrix(gamma1, nrow=1))
  k = ncol( matrix(gamma2, nrow=1))-1
  
  u<-max(p,k,r,s)
  
  N = n + burnin
  Y1<- array(1,N)
  Y2<- array(1,N)
  
  mu1<-array(1,c(N,1))
  mu2<-array(1,c(N,1))
  
  for(t in (u+1):N){
    
    # Conditional means
    mu1[t]<-exp(beta1%*%c(1,log(Y1[(t-1):(t-r)]+1))+
                  gamma1%*%log(Y2[(t-1):(t-s)]) )
    
    Y1[t]<- rpois(n=1,lambda=mu1[t])
    
    
    mu2[t]<-exp(beta2%*%c(1,log(Y2[(t-1):(t-p)])) +
                  gamma2%*%log(Y1[t:(t-k)]+1) )
    
    
    Y2[t]<-rgamma(n=1,shape=1/phi,scale=mu2[t]*phi)
  } 
  
  
  Y = cbind(Y1, Y2)
  Y = Y[(burnin+1):N, ]
  
  return(Y)
}

rGeoGeo = function(n, beta1, beta2, gamma1, gamma2, burnin=1000){
  
  p=ncol( matrix(beta2, nrow=1))-1
  r=ncol( matrix(beta1,nrow=1)) -1
  
  s=ncol( matrix(gamma1, nrow=1))
  k = ncol( matrix(gamma2, nrow=1))-1
  
  u<-max(p,k,r,s)
  
  N = n + burnin
  Y1<- array(1,N)
  Y2<- array(1,N)
  
  mu1<-array(1,c(N,1))
  mu2<-array(1,c(N,1))
  
  for(t in (u+1):N){
    
    # Conditional means
    mu1[t]<-beta1%*%c(1,log(Y1[(t-1):(t-r)]+1))+
                  gamma1%*%log(Y2[(t-1):(t-s)]+1) 
    
    Y1[t]<- rgeom(n=1, exp(mu1[t])/(1 + exp(mu1[t])))
    
    
    mu2[t]<-beta2%*%c(1,log(Y2[(t-1):(t-p)]+1)) +
                  gamma2%*%log(Y1[t:(t-k)]+1) 
    
    
    Y2[t]<-rgeom(n=1, exp(mu2[t])/(1 + exp(mu2[t])))
  } 
  
  
  Y = cbind(Y1, Y2)
  Y = Y[(burnin+1):N, ]
  
  return(Y)
}


op_ll1 = function(par, r, s, p, k, Y, delta_gammas){
  
  beta1 = par[1:(r+1)]
  gamma1 = par[(r+2):(r+s+1)]    
  phi = par[length(par)]
  
  beta2 = c(0.5, rep(0, p))
  gamma2 = c(0, rep(0, k+1))
  
  val =ll_Poisson_Gamma( beta1,  beta2,  gamma1, gamma2, phi, Y, delta_gammas)[1]
  return(val)
}



op_ll1_pois_pois = function(par, r, s, p, k, Y, delta_gammas){
  
  beta1 = par[1:(r+1)]
  gamma1 = par[(r+2):(r+s+1)]    
  
  beta2 = c(0.5, rep(0, p))
  gamma2 = c(0, rep(0, k+1))
  
  val =ll_Pois_Pois( beta1,  beta2,  gamma1, gamma2, Y, delta_gammas)[1]
  return(val)
}

op_ll1_geogeo = function(par, r, s, p, k, Y, delta_gammas){
  
  beta1 = par[1:(r+1)]
  gamma1 = par[(r+2):(r+s+1)]    
  
  beta2 = c(0.5, rep(0, p))
  gamma2 = c(0, rep(0, k+1))
  
  val =ll_Geo_Geo( beta1,  beta2,  gamma1, gamma2, Y, delta_gammas)[1]
  return(val)
}

op_ll1_geogamma = function(par, r, s, p, k, Y, delta_gammas){
  
  beta1 = par[1:(r+1)]
  gamma1 = par[(r+2):(r+s+1)]    
  phi = par[length(par)]
  
  beta2 = c(0.5, rep(0, p))
  gamma2 = c(0, rep(0, k+1))
  
  val =ll_Geo_Gamma( beta1,  beta2,  gamma1, gamma2, phi, Y, delta_gammas)[1]
  return(val)
}


op_ll2 = function(par,  r, p, Y, delta_betas, delta_gammas){
  
  
  s=length(delta_gammas[[1]])
  k=length(delta_gammas[[2]])-1
  beta1 =c(1,rep(0, r))
  gamma1 = rep(0,s)
  phi=1
  beta2 = par[1:(p+1)]
  gamma2 = par[(p+2):(p+2+k)]    
  
  
  val =ll_Poisson_Gamma( beta1,  beta2,  gamma1, gamma2, phi, Y, delta_gammas)[2]
  return(val)
}



op_ll2_geogeo = function(par,  r, p, Y, delta_betas, delta_gammas){
  
  
  s=length(delta_gammas[[1]])
  k=length(delta_gammas[[2]])-1
  beta1 =c(1,rep(0, r))
  gamma1 = rep(0,s)
  beta2 = par[1:(p+1)]
  gamma2 = par[(p+2):(p+2+k)]    
  
  
  val =ll_Geo_Geo( beta1,  beta2,  gamma1, gamma2, Y, delta_gammas)[2]
  return(val)
}


op_ll2_pois_pois = function(par,  r, p, Y, delta_betas, delta_gammas){
  
  
  s=length(delta_gammas[[1]])
  k=length(delta_gammas[[2]])-1
  beta1 =c(1,rep(0, r))
  gamma1 = rep(0,s)
  
  beta2 = par[1:(p+1)]
  gamma2 = par[(p+2):(p+2+k)]    
  
  
  val =ll_Pois_Pois( beta1,  beta2,  gamma1, gamma2, Y, delta_gammas)[2]
  return(val)
}


op_ll2_geogamma = function(par,  r, p, Y, delta_betas, delta_gammas){
  
  
  s=length(delta_gammas[[1]])
  k=length(delta_gammas[[2]])-1
  beta1 =c(1,rep(0, r))
  gamma1 = rep(0,s)
  phi=1
  beta2 = par[1:(p+1)]
  gamma2 = par[(p+2):(p+2+k)]    
  
  
  val =ll_Geo_Gamma( beta1,  beta2,  gamma1, gamma2, phi, Y, delta_gammas)[2]
  return(val)
}


fit_mcmc_pois_gamma = function(r,p,s,k,Y, N, priors, e, inits){
  
  delta_gammas = list();
  delta_gammas[[1]] = rep(1,s) 
  delta_gammas[[2]] = rep(1,k+1)
  
  if(missing(e)){
    e =0
  }
  if(missing(inits)){
  

  op=optim(par=c(1, rep(0, r), rep(0, s), 1),fn= op_ll1,  Y=Y,delta_gammas=delta_gammas,
           r=r,s=s,k=k,p=p,
           control = list(fnscale=-1), method = 'BFGS')
  
  
  beta1 = op$par[1:(r+1)] +rnorm(r+1, 0, e)
  gamma1 = op$par[(r+2):(r+s+1)] +rnorm(s, 0, e)
  phi = max( c(0.1, op$par[length(op$par)] +rnorm(1, 0, e)))
  
  
  op=optim(par=c(1, rep(0.01, p), rep(0, k+1)),fn= op_ll2,  Y=Y,delta_gammas=delta_gammas,
           r=r,p=p,
           control = list(fnscale=-1), method = 'BFGS')
  
  beta2 = op$par[1:(p+1)] +rnorm(p+1, 0, e)
  gamma2 = op$par[(p+2):(p+2+k)] +rnorm(k+1, 0, e)
  } else{
    beta1 = inits$beta1
    beta2 = inits$beta2
    gamma1 = inits$gamma1
    gamma2 = inits$gamma2
    phi = inits$phi
  }

  
  beta1_chain = matrix(NA, ncol = r+1, nrow = N)
  beta2_chain = matrix(NA, ncol = p+1, nrow = N)
  gamma1_chain = matrix(NA, ncol = s, nrow = N)
  gamma2_chain = matrix(NA, ncol = k+1, nrow = N)
  delta_g1 = numeric(s)
  delta_g2 = numeric(k+1)
  phi_chain = ll_chain = numeric(N)
  
  
  sigma_beta1 = rep(0.5, r+1)
  sigma_beta2 = rep(0.5, p+1)
  sigma_gamma1 =rep(0.5, s)
  sigma_gamma2 =rep(0.5, k+1)
  sigma_phi = 0.5
  
  
  accept_beta1 = rep(0, r+1)
  accept_beta2 = rep(0, p+1)
  accept_gamma1 = rep(0, s)
  accept_gamma2 = rep(0, k+1)
  accept_phi =0
  
  #### prior hyperparameters
  
  a= priors$a
  b = priors$b
  
  prior_sd_b1 = rep(priors$sd_b1, r+1)
  prior_sd_g1 = rep(priors$sd_g1, s)
  prior_sd_b2 = rep(priors$sd_b2, p+1)
  prior_sd_g2 = rep(priors$sd_g2, k+1)
  
  
  prior_rate_phi =priors$rate_phi
  
  start =Sys.time()
  iter =1
  while(iter <=N){
    
    
    ##### First series -----------------------------------------------------
    
    #Betas
    for(i1 in 1:(r+1)){
      
      mh_beta = mh_beta1( i1-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, sigma_beta1, prior_sd_b1);
      
      beta1 = mh_beta$param
      
      accept_beta1[i1] = accept_beta1[ i1] + mh_beta$accepted
      sigma_beta1[i1] =adjust_proposal( sigma_beta1[i1],   accept_beta1[ i1]/iter,  iter,  0.44)
    }
    
    
    #Omega1: Gamma
    omega_g1= omega_sampler(delta_gammas[[1]], a, b)
    
    #Deltas: Gamma
    for(g_index in 1:s){
      
      delta_gammas[[1]][g_index] = delta_g1_cpp(g_index-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, omega_g1)
      
      gamma1[g_index] = ifelse( delta_gammas[[1]][g_index] ==0, 0,  gamma1[g_index])
    }
    delta_g1= delta_g1+ delta_gammas[[1]]
    
    
    #Gamma
    include1 = 1:s
   # include1 = which(delta_gammas[[1]]==1)
   if(sum(delta_gammas[[1]])>0){
      for(i1 in 1:length(include1)){
        
        mh_gamma = mh_gamma1(include1[i1]-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, sigma_gamma1, prior_sd_g1);
        
        
        gamma1 = mh_gamma$param
        
        accept_gamma1[ include1[i1]] = accept_gamma1[ include1[i1]] + mh_gamma$accepted
        sigma_gamma1[  include1[i1]] =adjust_proposal( sigma_gamma1[ include1[i1]],   accept_gamma1[ include1[i1]]/iter,  iter,  0.44)
        
      }
    }
    
    ##### Second series -----------------------------------------------------
    
    
    #Betas
    
    for(i2 in 1:(p+1)){
      mh_beta =  mh_beta2( i2-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, sigma_beta2, prior_sd_b2);
      beta2 = mh_beta$param
      
      accept_beta2[i2] = accept_beta2[i2] + mh_beta$accepted
      
      sigma_beta2[ i2] =adjust_proposal( sigma_beta2[i2],   accept_beta2[i2]/iter,  iter,  0.44)
    }
    
    
    #Omega: Gamma2
    omega_g2= omega_sampler(delta_gammas[[2]], a, b)
    
    #Deltas: Gamma2
    for(g_index in 1:(k+1)){
      
      delta_gammas[[2]][g_index] = delta_g2_cpp(g_index-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, omega_g2)
      gamma2[g_index] = ifelse( delta_gammas[[2]][g_index] ==0, 0,  gamma2[g_index])
    }
    
    delta_g2= delta_g2+ delta_gammas[[2]]
    
    #Gammas
    include2 = which(delta_gammas[[2]]==1)
    include2 = 1:(k+1)
    if(sum(delta_gammas[[2]])>0){
      for(i2 in 1:length(include2)){
        
        mh_gamma =mh_gamma2(include2[i2]-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, sigma_gamma2, prior_sd_g2);
        
        gamma2 = mh_gamma$param
        
        accept_gamma2[include2[i2]] = accept_gamma2[include2[i2]] + mh_gamma$accepted
        sigma_gamma2[ include2[i2]] =adjust_proposal( sigma_gamma2[include2[i2]],   accept_gamma2[include2[i2]]/iter,  iter,  0.44)
        
     }
    }
    
    
    #mh_p = MH_phi(Y, beta1, gamma1, phi, sigma_phi, delta_betas, delta_gammas)
    
    mh_p = mh_phi( beta1, beta2,  gamma1, gamma2,  phi,  Y,  delta_gammas, sigma_phi, prior_rate_phi)
    phi = mh_p$param
    accept_phi = accept_phi + mh_p$accepted
    sigma_phi =adjust_proposal( sigma_phi,   accept_phi/iter,  iter,  0.44)
    ll_chain[iter] = mh_p$ll
    
    
    #Store 
    beta1_chain[iter,] = beta1
    beta2_chain[iter,] = beta2
    gamma1_chain[iter,] = gamma1
    gamma2_chain[iter,] = gamma2
    phi_chain[iter] = phi
    
    iter = iter +1
   
    if (iter %% (N / 10) == 0) {
      pct <- iter / N
      bar_width <- 50
      filled <- round(pct * bar_width)
      
      bar <- paste0(strrep("=", filled), strrep(" ", bar_width - filled))
      
      # print progress bar only
      cat(sprintf("\r[%s] %3.0f%%", bar, 100 * pct))
      flush.console()
      
      # newline at the end
      if (iter == N) cat("\nDone!\n")
    }
    
    
  }
  end=Sys.time()
  
  output = list("beta1" = beta1_chain,
                "beta2" = beta2_chain,
                "gamma1" = gamma1_chain,
                "gamma2" = gamma2_chain,
                "phi" = phi_chain, "delta_gamma1" = delta_g1/N,
                "delta_gamma2" = delta_g2/N, "time" = as.numeric(end-start,units="secs"),
                "loglik" = ll_chain)
  
  
  return(output)
  
}



fit_mcmc_geogamma = function(r,p,s,k,Y, N, priors, e, inits){
  
  delta_gammas = list();
  delta_gammas[[1]] = rep(1,s) 
  delta_gammas[[2]] = rep(1,k+1)
  
  if(missing(e)){
    e =0
  }
  if(missing(inits)){
    
    
    op=optim(par=c(1, rep(0, r), rep(0, s), 1),fn= op_ll1_geogamma,  Y=Y,delta_gammas=delta_gammas,
             r=r,s=s,k=k,p=p,
             control = list(fnscale=-1), method = 'BFGS')
    
    
    beta1 = op$par[1:(r+1)] +rnorm(r+1, 0, e)
    gamma1 = op$par[(r+2):(r+s+1)] +rnorm(s, 0, e)
    phi = max( c(0.1, op$par[length(op$par)] +rnorm(1, 0, e)))
    
    
    op=optim(par=c(1, rep(0, p), rep(0, k+1)),fn= op_ll2_geogamma,  Y=Y,delta_gammas=delta_gammas,
             r=r,p=p,
             control = list(fnscale=-1), method = 'BFGS')
    
    beta2 = op$par[1:(p+1)] +rnorm(p+1, 0, e)
    gamma2 = op$par[(p+2):(p+2+k)] +rnorm(k+1, 0, e)
  } else{
    beta1 = inits$beta1
    beta2 = inits$beta2
    gamma1 = inits$gamma1
    gamma2 = inits$gamma2
    phi = inits$phi
  }
  
  beta1_chain = matrix(NA, ncol = r+1, nrow = N)
  beta2_chain = matrix(NA, ncol = p+1, nrow = N)
  gamma1_chain = matrix(NA, ncol = s, nrow = N)
  gamma2_chain = matrix(NA, ncol = k+1, nrow = N)
  delta_g1 = numeric(s)
  delta_g2 = numeric(k+1)
  phi_chain = ll_chain = numeric(N)
  
  
  sigma_beta1 = rep(0.5, r+1)
  sigma_beta2 = rep(0.5, p+1)
  sigma_gamma1 =rep(0.5, s)
  sigma_gamma2 =rep(0.5, k+1)
  sigma_phi = 0.5
  
  
  accept_beta1 = rep(0, r+1)
  accept_beta2 = rep(0, p+1)
  accept_gamma1 = rep(0, s)
  accept_gamma2 = rep(0, k+1)
  accept_phi =0
  
  #### prior hyperparameters
  
  a= priors$a
  b = priors$b
  
  prior_sd_b1 = rep(priors$sd_b1, r+1)
  prior_sd_g1 = rep(priors$sd_g1, s)
  prior_sd_b2 = rep(priors$sd_b2, p+1)
  prior_sd_g2 = rep(priors$sd_g2, k+1)
  
  
  prior_rate_phi =priors$rate_phi
  
  start =Sys.time()
  iter =1
  while(iter <=N){
    
    
    ##### First series -----------------------------------------------------
    
    #Betas
    for(i1 in 1:(r+1)){
      
      mh_beta = mh_beta1_geogamma( i1-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, sigma_beta1, prior_sd_b1);
      
      beta1 = mh_beta$param
      
      accept_beta1[i1] = accept_beta1[ i1] + mh_beta$accepted
      sigma_beta1[i1] =adjust_proposal( sigma_beta1[i1],   accept_beta1[ i1]/iter,  iter,  0.44)
    }
    
    
    #Omega1: Gamma
    omega_g1= omega_sampler(delta_gammas[[1]], a, b)
    
    #Deltas: Gamma
    for(g_index in 1:s){
      
      delta_gammas[[1]][g_index] = delta_g1_cpp_geogamma(g_index-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, omega_g1)
      
      gamma1[g_index] = ifelse( delta_gammas[[1]][g_index] ==0, 0,  gamma1[g_index])
    }
    delta_g1= delta_g1+ delta_gammas[[1]]
    
    
    #Gamma
    include1 = 1:s
    # include1 = which(delta_gammas[[1]]==1)
    if(sum(delta_gammas[[1]])>0){
      for(i1 in 1:length(include1)){
        
        mh_gamma = mh_gamma1_geogamma(include1[i1]-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, sigma_gamma1, prior_sd_g1);
        
        
        gamma1 = mh_gamma$param
        
        accept_gamma1[ include1[i1]] = accept_gamma1[ include1[i1]] + mh_gamma$accepted
        sigma_gamma1[  include1[i1]] =adjust_proposal( sigma_gamma1[ include1[i1]],   accept_gamma1[ include1[i1]]/iter,  iter,  0.44)
        
      }
    }
    
    ##### Second series -----------------------------------------------------
    
    
    #Betas
    
    for(i2 in 1:(p+1)){
      mh_beta =  mh_beta2_geogamma( i2-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, sigma_beta2, prior_sd_b2);
      beta2 = mh_beta$param
      
      accept_beta2[i2] = accept_beta2[i2] + mh_beta$accepted
      
      sigma_beta2[ i2] =adjust_proposal( sigma_beta2[i2],   accept_beta2[i2]/iter,  iter,  0.44)
    }
    
    
    #Omega: Gamma2
    omega_g2= omega_sampler(delta_gammas[[2]], a, b)
    
    #Deltas: Gamma2
    for(g_index in 1:(k+1)){
      
      delta_gammas[[2]][g_index] = delta_g2_cpp_geogamma(g_index-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, omega_g2)
      gamma2[g_index] = ifelse( delta_gammas[[2]][g_index] ==0, 0,  gamma2[g_index])
    }
    
    delta_g2= delta_g2+ delta_gammas[[2]]
    
    #Gammas
    include2 = which(delta_gammas[[2]]==1)
    include2 = 1:(k+1)
    if(sum(delta_gammas[[2]])>0){
      for(i2 in 1:length(include2)){
        
        mh_gamma =mh_gamma2_geogamma(include2[i2]-1, beta1, beta2, gamma1, gamma2, phi, Y, delta_gammas, sigma_gamma2, prior_sd_g2);
        
        gamma2 = mh_gamma$param
        
        accept_gamma2[include2[i2]] = accept_gamma2[include2[i2]] + mh_gamma$accepted
        sigma_gamma2[ include2[i2]] =adjust_proposal( sigma_gamma2[include2[i2]],   accept_gamma2[include2[i2]]/iter,  iter,  0.44)
        
      }
    }
    
    
    #mh_p = MH_phi(Y, beta1, gamma1, phi, sigma_phi, delta_betas, delta_gammas)
    
    mh_p = mh_phi_geogamma( beta1, beta2,  gamma1, gamma2,  phi,  Y,  delta_gammas, sigma_phi, prior_rate_phi)
    phi = mh_p$param
    accept_phi = accept_phi + mh_p$accepted
    sigma_phi =adjust_proposal( sigma_phi,   accept_phi/iter,  iter,  0.44)
    ll_chain[iter] = mh_p$ll
    
    
    #Store 
    beta1_chain[iter,] = beta1
    beta2_chain[iter,] = beta2
    gamma1_chain[iter,] = gamma1
    gamma2_chain[iter,] = gamma2
    phi_chain[iter] = phi
    
    iter = iter +1
    
    
    if (iter %% (N / 10) == 0) {
      pct <- iter / N
      bar_width <- 50
      filled <- round(pct * bar_width)
      
      bar <- paste0(strrep("=", filled), strrep(" ", bar_width - filled))
      
      # print progress bar only
      cat(sprintf("\r[%s] %3.0f%%", bar, 100 * pct))
      flush.console()
      
      # newline at the end
      if (iter == N) cat("\nDone!\n")
    }
    
    
  }
  end=Sys.time()
  
  output = list("beta1" = beta1_chain,
                "beta2" = beta2_chain,
                "gamma1" = gamma1_chain,
                "gamma2" = gamma2_chain,
                "phi" = phi_chain, "delta_gamma1" = delta_g1/N,
                "delta_gamma2" = delta_g2/N, "time" = as.numeric(end-start,units="secs"),
                "loglik" = ll_chain)
  
  
  return(output)
  
}







fit_mcmc_geogeo = function(r,p,s,k,Y, N, priors, e, inits){
  
  delta_gammas = list();
  delta_gammas[[1]] = rep(1,s) 
  delta_gammas[[2]] = rep(1,k+1)
  
  if(missing(e)){
    e =0
  }
  if(missing(inits)){
    
    
    op=optim(par=c(1, rep(0, r), rep(0, s)),fn= op_ll1_geogeo,  Y=Y,delta_gammas=delta_gammas,
             r=r,s=s,k=k,p=p,
             control = list(fnscale=-1), method = 'BFGS')

    
    beta1 = op$par[1:(r+1)] +rnorm(r+1, 0, e)
    gamma1 = op$par[(r+2):(r+s+1)] +rnorm(s, 0, e)
    
    
    op=optim(par=c(1, rep(0, p), rep(0, k+1)),fn= op_ll2_geogeo,  Y=Y,delta_gammas=delta_gammas,
             r=r,p=p,
             control = list(fnscale=-1), method = 'BFGS')
    
    beta2 = op$par[1:(p+1)] +rnorm(p+1, 0, e)
    gamma2 = op$par[(p+2):(p+2+k)] +rnorm(k+1, 0, e)
  } else{
    beta1 = inits$beta1
    beta2 = inits$beta2
    gamma1 = inits$gamma1
    gamma2 = inits$gamma2

  }
  
  beta1_chain = matrix(NA, ncol = r+1, nrow = N)
  beta2_chain = matrix(NA, ncol = p+1, nrow = N)
  gamma1_chain = matrix(NA, ncol = s, nrow = N)
  gamma2_chain = matrix(NA, ncol = k+1, nrow = N)
  delta_g1 = numeric(s)
  delta_g2 = numeric(k+1)
  ll_chain = numeric(N)
  
  
  sigma_beta1 = rep(0.5, r+1)
  sigma_beta2 = rep(0.5, p+1)
  sigma_gamma1 =rep(0.5, s)
  sigma_gamma2 =rep(0.5, k+1)
  
  accept_beta1 = rep(0, r+1)
  accept_beta2 = rep(0, p+1)
  accept_gamma1 = rep(0, s)
  accept_gamma2 = rep(0, k+1)
  

  #### prior hyperparameters
  
  a= priors$a
  b = priors$b
  
  prior_sd_b1 = rep(priors$sd_b1, r+1)
  prior_sd_g1 = rep(priors$sd_g1, s)
  prior_sd_b2 = rep(priors$sd_b2, p+1)
  prior_sd_g2 = rep(priors$sd_g2, k+1)
  
  
  start =Sys.time()
  iter =1
  while(iter <=N){
    
    
    ##### First series -----------------------------------------------------
    
    #Betas
    for(i1 in 1:(r+1)){
      
      mh_beta = mh_beta1_geo_geo( i1-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, sigma_beta1, prior_sd_b1);
      
      beta1 = mh_beta$param
      
      accept_beta1[i1] = accept_beta1[ i1] + mh_beta$accepted
      sigma_beta1[i1] =adjust_proposal( sigma_beta1[i1],   accept_beta1[ i1]/iter,  iter,  0.44)
    }
    
    
    #Omega1: Gamma
    omega_g1= omega_sampler(delta_gammas[[1]], a, b)
    
    #Deltas: Gamma
    for(g_index in 1:s){
      
      delta_gammas[[1]][g_index] = delta_g1_cpp_geogeo(g_index-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, omega_g1)
      
      gamma1[g_index] = ifelse( delta_gammas[[1]][g_index] ==0, 0,  gamma1[g_index])
    }
    delta_g1= delta_g1+ delta_gammas[[1]]
    
    
    #Gamma
    include1 = 1:s
    # include1 = which(delta_gammas[[1]]==1)
    if(sum(delta_gammas[[1]])>0){
      for(i1 in 1:length(include1)){
        
        mh_gamma = mh_gamma1_geo_geo(include1[i1]-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, sigma_gamma1, prior_sd_g1);
        
        
        gamma1 = mh_gamma$param
        
        accept_gamma1[ include1[i1]] = accept_gamma1[ include1[i1]] + mh_gamma$accepted
        sigma_gamma1[  include1[i1]] =adjust_proposal( sigma_gamma1[ include1[i1]],   accept_gamma1[ include1[i1]]/iter,  iter,  0.44)
        
      }
    }
    
    ##### Second series -----------------------------------------------------
    
    
    #Betas
    
    for(i2 in 1:(p+1)){
      mh_beta =  mh_beta2_geo_geo( i2-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, sigma_beta2, prior_sd_b2);
      beta2 = mh_beta$param
      
      accept_beta2[i2] = accept_beta2[i2] + mh_beta$accepted
      
      sigma_beta2[ i2] =adjust_proposal( sigma_beta2[i2],   accept_beta2[i2]/iter,  iter,  0.44)
    }
    
    
    #Omega: Gamma2
    omega_g2= omega_sampler(delta_gammas[[2]], a, b)
    
    #Deltas: Gamma2
    for(g_index in 1:(k+1)){
      
      delta_gammas[[2]][g_index] = delta_g2_cpp_geogeo(g_index-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, omega_g2)
      gamma2[g_index] = ifelse( delta_gammas[[2]][g_index] ==0, 0,  gamma2[g_index])
    }
    
    delta_g2= delta_g2+ delta_gammas[[2]]
    
    #Gammas
    include2 = which(delta_gammas[[2]]==1)
    include2 = 1:(k+1)
    if(sum(delta_gammas[[2]])>0){
      for(i2 in 1:length(include2)){
        
        mh_gamma =mh_gamma2_geo_geo(include2[i2]-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, sigma_gamma2, prior_sd_g2);
        
        gamma2 = mh_gamma$param
        
        accept_gamma2[include2[i2]] = accept_gamma2[include2[i2]] + mh_gamma$accepted
        sigma_gamma2[ include2[i2]] =adjust_proposal( sigma_gamma2[include2[i2]],   accept_gamma2[include2[i2]]/iter,  iter,  0.44)
        
      }
    }
    

    ll_chain[iter] = sum(ll_Geo_Geo(beta1, beta2, gamma1, gamma2, Y, delta_gammas))
    
    
    #Store 
    beta1_chain[iter,] = beta1
    beta2_chain[iter,] = beta2
    gamma1_chain[iter,] = gamma1
    gamma2_chain[iter,] = gamma2
    
    iter = iter +1
    
    
    if (iter %% (N / 10) == 0) {
      pct <- iter / N
      bar_width <- 50
      filled <- round(pct * bar_width)
      
      bar <- paste0(strrep("=", filled), strrep(" ", bar_width - filled))
      
      # print progress bar only
      cat(sprintf("\r[%s] %3.0f%%", bar, 100 * pct))
      flush.console()
      
      # newline at the end
      if (iter == N) cat("\nDone!\n")
    }
    
    
  }
  end=Sys.time()
  
  output = list("beta1" = beta1_chain,
                "beta2" = beta2_chain,
                "gamma1" = gamma1_chain,
                "gamma2" = gamma2_chain,
                "delta_gamma1" = delta_g1/N,
                "delta_gamma2" = delta_g2/N, "time" = as.numeric(end-start,units="secs"),
                "loglik" = ll_chain)
  
  
  return(output)
  
}



fit_mcmc_poispois = function(r,p,s,k,Y, N, priors, e, inits){
  
  delta_gammas = list();
  delta_gammas[[1]] = rep(1,s) 
  delta_gammas[[2]] = rep(1,k+1)
  
  
  if(missing(e)){
    e =0
  }
  if(missing(inits)){
    
    
    op=optim(par=c(1, rep(0, r), rep(0, s)),fn= op_ll1_pois_pois,  Y=Y,delta_gammas=delta_gammas,
             r=r,s=s,k=k,p=p,
             control = list(fnscale=-1), method = 'BFGS')
    
    beta1 = op$par[1:(r+1)] +rnorm(r+1, 0, e)
    gamma1 = op$par[(r+2):(r+s+1)] +rnorm(s, 0, e)
    
    
    op=optim(par=c(0.1, rep(0, p), rep(0, k+1)),fn= op_ll2_pois_pois,  Y=Y,delta_gammas=delta_gammas,
             r=r,p=p,
             control = list(fnscale=-1), method = 'BFGS')
    
    beta2 = op$par[1:(p+1)] +rnorm(p+1, 0, e)
    gamma2 = op$par[(p+2):(p+2+k)] +rnorm(k+1, 0, e)
  } else{
    beta1 = inits$beta1
    beta2 = inits$beta2
    gamma1 = inits$gamma1
    gamma2 = inits$gamma2

  }
  
  beta1_chain = matrix(NA, ncol = r+1, nrow = N)
  beta2_chain = matrix(NA, ncol = p+1, nrow = N)
  gamma1_chain = matrix(NA, ncol = s, nrow = N)
  gamma2_chain = matrix(NA, ncol = k+1, nrow = N)
  delta_g1 = numeric(s)
  delta_g2 = numeric(k+1)
  ll_chain = numeric(N)
  
  
  sigma_beta1 = rep(0.5, r+1)
  sigma_beta2 = rep(0.5, p+1)
  sigma_gamma1 =rep(0.5, s)
  sigma_gamma2 =rep(0.5, k+1)

  
  accept_beta1 = rep(0, r+1)
  accept_beta2 = rep(0, p+1)
  accept_gamma1 = rep(0, s)
  accept_gamma2 = rep(0, k+1)
  
  #### prior hyperparameters
  
  a= priors$a
  b = priors$b
  
  prior_sd_b1 = rep(priors$sd_b1, r+1)
  prior_sd_g1 = rep(priors$sd_g1, s)
  prior_sd_b2 = rep(priors$sd_b2, p+1)
  prior_sd_g2 = rep(priors$sd_g2, k+1)

  
  start =Sys.time()
  iter =1
  while(iter <=N){
    
    
    ##### First series -----------------------------------------------------
    
    #Betas
    for(i1 in 1:(r+1)){
      
      mh_beta = mh_beta1_pois_pois( i1-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, sigma_beta1, prior_sd_b1);
      
      beta1 = mh_beta$param
      
      accept_beta1[i1] = accept_beta1[ i1] + mh_beta$accepted
      sigma_beta1[i1] =adjust_proposal( sigma_beta1[i1],   accept_beta1[ i1]/iter,  iter,  0.44)
    }
    
    
    #Omega1: Gamma
    omega_g1= omega_sampler(delta_gammas[[1]], a, b)
    
    #Deltas: Gamma
    for(g_index in 1:s){
      
      delta_gammas[[1]][g_index] = delta_g1_cpp_poispois(g_index-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, omega_g1)
      
      gamma1[g_index] = ifelse( delta_gammas[[1]][g_index] ==0, 0,  gamma1[g_index])
    }
    delta_g1= delta_g1+ delta_gammas[[1]]
    
    
    #Gamma
    include1 = 1:s
    # include1 = which(delta_gammas[[1]]==1)
    if(sum(delta_gammas[[1]])>0){
      for(i1 in 1:length(include1)){
        
        mh_gamma = mh_gamma1_pois_pois(include1[i1]-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, sigma_gamma1, prior_sd_g1);
        
        
        gamma1 = mh_gamma$param
        
        accept_gamma1[ include1[i1]] = accept_gamma1[ include1[i1]] + mh_gamma$accepted
        sigma_gamma1[  include1[i1]] =adjust_proposal( sigma_gamma1[ include1[i1]],   accept_gamma1[ include1[i1]]/iter,  iter,  0.44)
        
      }
    }
    
    ##### Second series -----------------------------------------------------
    
    
    #Betas
    
    for(i2 in 1:(p+1)){
      mh_beta =  mh_beta2_pois_pois( i2-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, sigma_beta2, prior_sd_b2);
      beta2 = mh_beta$param
      
      accept_beta2[i2] = accept_beta2[i2] + mh_beta$accepted
      
      sigma_beta2[ i2] =adjust_proposal( sigma_beta2[i2],   accept_beta2[i2]/iter,  iter,  0.44)
    }
    
    
    #Omega: Gamma2
    omega_g2= omega_sampler(delta_gammas[[2]], a, b)
    
    #Deltas: Gamma2
    for(g_index in 1:(k+1)){
      
      delta_gammas[[2]][g_index] = delta_g2_cpp_poispois(g_index-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, omega_g2)
      gamma2[g_index] = ifelse( delta_gammas[[2]][g_index] ==0, 0,  gamma2[g_index])
    }
    
    delta_g2= delta_g2+ delta_gammas[[2]]
    
    #Gammas
    include2 = which(delta_gammas[[2]]==1)
    include2 = 1:(k+1)
    if(sum(delta_gammas[[2]])>0){
      for(i2 in 1:length(include2)){
        
        mh_gamma =mh_gamma2_pois_pois(include2[i2]-1, beta1, beta2, gamma1, gamma2, Y, delta_gammas, sigma_gamma2, prior_sd_g2);
        
        gamma2 = mh_gamma$param
        
        accept_gamma2[include2[i2]] = accept_gamma2[include2[i2]] + mh_gamma$accepted
        sigma_gamma2[ include2[i2]] =adjust_proposal( sigma_gamma2[include2[i2]],   accept_gamma2[include2[i2]]/iter,  iter,  0.44)
        
      }
    }
    

    ll_chain[iter] =  sum(ll_Pois_Pois(beta1, beta2, gamma1, gamma2, Y, delta_gammas))
    
    
    #Store 
    beta1_chain[iter,] = beta1
    beta2_chain[iter,] = beta2
    gamma1_chain[iter,] = gamma1
    gamma2_chain[iter,] = gamma2
    
    iter = iter +1
    
    if (iter %% (N / 10) == 0) {
      pct <- iter / N
      bar_width <- 50
      filled <- round(pct * bar_width)
      
      bar <- paste0(strrep("=", filled), strrep(" ", bar_width - filled))
      
      # print progress bar only
      cat(sprintf("\r[%s] %3.0f%%", bar, 100 * pct))
      flush.console()
      
      # newline at the end
      if (iter == N) cat("\nDone!\n")
    }
    
    
  }
  end=Sys.time()
  
  output = list("beta1" = beta1_chain,
                "beta2" = beta2_chain,
                "gamma1" = gamma1_chain,
                "gamma2" = gamma2_chain,
                "delta_gamma1" = delta_g1/N,
                "delta_gamma2" = delta_g2/N, 
                "time" = as.numeric(end-start,units="secs"),
                "loglik" = ll_chain)
  
  
  return(output)
  
}


fit_mcmc <- function(r, p, s, k, Y, N, priors, dist1, dist2){

  # Helper functions for checking data types
  is_positive_continuous <- function(x) all(x > 0 & x %% 1 != 0)
  is_integer <- function(x) all(x %% 1 == 0)
  
  # Check Y depending on distributions
  if (dist1 == "gamma" && dist2 == "poisson") {
    if (!is_positive_continuous(Y[,1])) stop("Y[,1] must be positive continuous for Gamma")
    if (!is_integer(Y[,2])) stop("Y[,2] must be integer for Poisson")
    fit_fun <- fit_mcmc_pois_gamma
  } else if (dist1 == "gamma" && dist2 == "geo") {
    if (!is_positive_continuous(Y[,1])) stop("Y[,1] must be positive continuous for Gamma")
    if (!is_integer(Y[,2])) stop("Y[,2] must be integer for Geometric")
    fit_fun <- fit_mcmc_geogamma
  } else if (dist1 == "geo" && dist2 == "geo") {
    if (!is_integer(Y[,1])) stop("Y[,1] must be integer for Geometric")
    if (!is_integer(Y[,2])) stop("Y[,2] must be integer for Geometric")
    fit_fun <- fit_mcmc_geogeo
  } else if (dist1 == "poisson" && dist2 == "poisson") {
    if (!is_integer(Y[,1])) stop("Y[,1] must be integer for Poisson")
    if (!is_integer(Y[,2])) stop("Y[,2] must be integer for Poisson")
    fit_fun <- fit_mcmc_poispois
  } else {
    stop("Invalid combination of dist1 and dist2")
  }
  
  # Call the appropriate fitting function
  fit <- fit_fun(r = r, p = p, s = s, k = k, Y = Y, N = N, priors = priors, e = 0)
  
  return(fit)
}




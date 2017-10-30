.COMPASS.discrete <- function(n_s, n_u, categories, iterations, replications,
  verbose=TRUE, ...) {

  vmessage <- function(...)
    if (verbose) message(...) else invisible(NULL)

  ## Initial parameters
  vmessage("Initializing parameters...")
  N <- iterations ## The number of iterations to run the model
  ttt <- 2000 ## The step size in mode search (fixed)
  SS <- 1L

  N_s <- rowSums(n_s)
  N_u <- rowSums(n_u)

  I <- nrow(n_s)
  K <- ncol(n_u)
  K1 <- K - 1L

  indi = array(1, dim=c(I,K)) # 0 indicate that gamma_ik=0
  for (k in 1:K1) {
    l2 = which((n_s[,k]/N_s)-(n_u[,k]/N_u) <= 0)
    indi[l2, k] <- 0;
  }
  indi[,K] <- rowSums(indi[,1:K1, drop=FALSE])
  indi <- matrix(as.integer(indi), nrow=I)

  #############################################
  mk = array(as.integer(0), dim = c(1,K1));
  Istar = 0;
  mKstar = 0;

  gamma =array(as.integer(0), dim=c(I,K,N));

  alpha_u = array(0, dim=c(N,K));
  alpha_s = array(0, dim=c(N,K));

  varp_s1 = array(sqrt(3),dim=c(K,1)); # sqrt(var)
  varp_s2 = array(sqrt(10),dim=c(K,1)); # sqrt(var)
  varp_s2[K] = sqrt(20);
  varp_s1[K] = sqrt(10);

  pvar_s = array(0.8,dim=c(K,1));
  pvar_s[K] = 0.6;
  varp_u = array(sqrt(8),dim=c(K,1));

  pp = array(0.65, dim=c(I,1))
  pb1 <- clamp( 1.5 / median( indi[, K] ), 0, 0.9 )
  pb2 <- clamp( 5.0 / median( indi[, K] ), 0, 0.9 )
  lambda_s = rep(0,K);
  lambda_s[1:K1] = (10^-2)*max(N_s,N_u)
  lambda_s[K] = max(N_s,N_u)-sum(lambda_s[1:K1])
  lambda_u = lambda_s

  alpha_u[1,1:(K-1)] =10
  alpha_u[1,K] = 150

  alpha_s[1,1:(K-1)] = 10
  alpha_s[1,K] = 100

  #################### acceptance rate ###########################
  A_gm = array(as.integer(0), dim=c(I,N));
  A_alphau = array(as.integer(0), dim=c(K,N));
  A_alphas = array(as.integer(0), dim=c(K,N));

  vmessage("Computing initial parameter estimates...")
  for (tt in 2:N) {

    if (tt %% 1000 == 0) vmessage("Iteration ", tt, " of ", iterations, ".")

    # update alphau
    res2 <- .Call(C_updatealphau_noPu_Exp, alphaut = alpha_u[tt-1,],n_s = n_s,n_u=n_u, I=I, K=K, lambda_u = lambda_u, var_p = varp_u, ttt = ttt,gammat =gamma[,,tt-1])
    if(length(alpha_u[tt,])!=length(res2$alphau_tt)){
      vmessage("res2 alphau_tt length:",length(res2$alphau_tt),"\n")    
      vmessage("alpha_u[tt,] length: ",length(alpha_u[tt,]),"\n")
    }
    alpha_u[tt,] = res2$alphau_tt;
    A_alphau[,tt] = res2$Aalphau;

    #update gamma
    res1 <- .Call(C_updategammak_noPu, n_s = n_s,n_u=n_u,gammat = gamma[,,tt-1],I=I,K=K,SS = SS,alphau = alpha_u[tt,],alphas = alpha_s[tt-1,],alpha=1,mk=mk,Istar = Istar,
      mKstar = mKstar,pp=pp, pb1 = pb1, pb2 = pb2, indi=indi)
    gamma[,,tt] = res1$gamma_tt;
    if(length(A_gm[,tt]) != length(res1$Ag)){
	vmessage("res1 Ag length: ",length(res1$Ag),"\n")
	vmessage("A_gm[,tt] length: ",length(A_gm[,tt]),"\n")
    1}
    A_gm[,tt] = res1$Ag;
    Istar = res1$xIstar;
    mk = res1$xmk;
    mKstar = res1$mKstar;

    # update alphas
    res3 <- .Call(C_updatealphas_Exp, alphast = alpha_s[tt-1,], n_s = n_s,  K=K, I=I, lambda_s = lambda_s, gammat =gamma[,,tt], var_1 = varp_s1,var_2 = varp_s2,p_var = pvar_s,ttt = ttt)
    if(length(alpha_s[tt,])!=length(res3$alphas_tt)){
      vmessage("res3 alphas_tt length:",length(res3$alphas_tt),"\n")
      vmessage("alpha_s[tt,] length: ",length(alpha_s[tt,]),"\n")
    }
    if(length(A_alphas[,tt])!=length(res3$Aalphas)){
      vmessage("res3 Aalphas length:",length(res3$Aalphas),"\n")
      vmessage("A_alphas[,tt] length: ",length(A_alphas[,tt]),"\n")
    }


    alpha_s[tt,] = res3$alphas_tt;
    A_alphas[,tt] = res3$Aalphas;
    ####
    if (tt>4000&&(tt %% 4000)==0) {
      intr = (tt-1000+1):tt
      for (kk in 1:K) {
        Au = mean(A_alphau[kk,intr])
        if (Au<0.001) {varp_u[kk] = varp_u[kk]*sqrt(0.1);}
        else if (Au <0.05) {varp_u[kk] = varp_u[kk]*sqrt(0.5);}
        else if (Au < 0.2) {varp_u[kk] = varp_u[kk]*sqrt(0.9);}
        else if (Au>0.5) {varp_u[kk] = varp_u[kk]*sqrt(1.1);}
        else if (Au>0.75) {varp_u[kk] = varp_u[kk]*sqrt(2);}
        else if (Au>0.95) {varp_u[kk] = varp_u[kk]*sqrt(10);}
        As = mean(A_alphas[kk,intr])
        if (As<0.001) {varp_s1[kk] = varp_s1[kk]*sqrt(0.1);}
        else if (As <0.05) {varp_s1[kk] = varp_s1[kk]*sqrt(0.5);}
        else if (As < 0.2) {varp_s1[kk] = varp_s1[kk]*sqrt(0.9);}
        else if (As>0.5) {varp_s2[kk] = varp_s2[kk]*sqrt(1.1);}
        else if (As>0.75) {varp_s2[kk] = varp_s2[kk]*sqrt(2);}
        else if (As>0.95) {varp_s2[kk] = varp_s2[kk]*sqrt(10);}
      }
      for ( i in 1:I) {
        Agm = mean(A_gm[i,intr]);
        if (Agm<0.001) {pp[i] = min(0.9,pp[i]*1.4);}
        else if (Agm<0.05) {pp[i] = min(0.9,pp[i]*1.2);}
        else if (Agm<0.2) {pp[i] = min(0.9,pp[i]*1.1);}
        else if (Agm>0.6) {pp[i] = max(0.1,pp[i]*0.8);}
        else if (Agm>0.75) {pp[i] = max(0.1,pp[i]*0.5);}
        else if (Agm>0.95) {pp[i] = max(0.1,pp[i]*0.2);}
      }
    }
  }
  alpha_u[1,] = alpha_u[N,]
  gamma[,,1] = gamma[,,N]
  alpha_s[1,] = alpha_s[N,]

  sNN=replications ## number of 'replications' -- should be user defined
  vmessage("Fitting model with ", sNN, " replications.")
  for (stt in 1:sNN) {
    vmessage("Running replication ", stt, " of ", sNN, "...")
    for (tt in 2:N) {
      # update alphau
      res2 <- .Call(C_updatealphau_noPu_Exp, alphaut = alpha_u[tt-1,],n_s = n_s,n_u=n_u, I=I, K=K, lambda_u = lambda_u, var_p = varp_u, ttt = ttt,gammat =gamma[,,tt-1])
      alpha_u[tt,] = res2$alphau_tt;
      A_alphau[,tt] = res2$Aalphau;

      #update gamma
      res1 <- .Call(C_updategammak_noPu,n_s = n_s,n_u=n_u,gammat = gamma[,,tt-1],I=I,K=K,SS = SS,alphau = alpha_u[tt,],alphas = alpha_s[tt-1,],alpha=1,mk=mk,Istar = Istar,
        mKstar = mKstar,pp=pp, pb1 = pb1, pb2 = pb2, indi=indi)
      gamma[,,tt] = res1$gamma_tt;
      A_gm[,tt] = res1$Ag;
      Istar = res1$xIstar;
      mk = res1$xmk;
      mKstar = res1$mKstar;

      # update alphas
      res3 <- .Call(C_updatealphas_Exp, alphast = alpha_s[tt-1,], n_s = n_s,  K=K, I=I, lambda_s = lambda_s, gammat =gamma[,,tt], var_1 = varp_s1,var_2 = varp_s2,p_var = pvar_s,ttt = ttt)

    if(length(alpha_s[tt,])!=length(res3$alphas_tt)){
      vmessage("res3 alphas_tt length:",length(res3$alphas_tt),"\n")
      vmessage("alpha_s[tt,] length: ",length(alpha_s[tt,]),"\n")
    }
    if(length(A_alphas[,tt])!=length(res3$Aalphas)){
      vmessage("res3 Aalphas length:",length(res3$Aalphas),"\n")
      vmessage("A_alphas[,tt] length: ",length(A_alphas[,tt]),"\n")
    }
      alpha_s[tt,] = res3$alphas_tt;
      A_alphas[,tt] = res3$Aalphas;
      if (tt %% 1000 == 0) vmessage("Iteration ", tt, " of ", iterations, ".")
    }
    if (stt == sNN) {break}
    alpha_u[1,] = alpha_u[N,];
    gamma[,,1] = gamma[,,N];
    alpha_s[1,] = alpha_s[N,];
  }

  ######################################

  dimnames(gamma) <- list(rownames(n_s), NULL, NULL)
  Nburn=0;
  Mgamma = mat.or.vec(I,K);
  Nseq = seq(Nburn+1,N,by=1)
  for (ttt in Nseq) {
    Mgamma = Mgamma + gamma[,,ttt]; #thining
  }
  Mgamma = Mgamma/(N-Nburn);

  vmessage("Done!")

  cats <- categories[, -ncol(categories), drop=FALSE]
  subsets_df <- as.data.frame(cats)
  for (i in seq_along(subsets_df)) {
    tmp <- subsets_df[[i]]
    subsets_df[[i]] <- paste0(
      swap(tmp, c(0, 1), c("!", "")),
      colnames(subsets_df)[[i]]
    )
  }
  subsets <- do.call( function(...) paste(..., sep="&"), subsets_df )

  colnames(Mgamma) <- subsets

  ## set names on the output
  output <- list(
    alpha_s=alpha_s,
    A_alphas=rowMeans(A_alphas),
    alpha_u=alpha_u,
    A_alphau=rowMeans(A_alphau),
    gamma=gamma,
    mean_gamma=Mgamma,
    A_gamma=rowMeans(A_gm),
    categories=categories,
    model="discrete"
  )

  return(output)

}

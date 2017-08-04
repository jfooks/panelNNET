OLStrick_function_2 <- function(parlist, hidden_layers, y, fe_var, lam, parapen, treatment){
  
  # iter <- pnn$iter
  # maxit <- pnn$maxit 
  # stopcounter <- pnn$stopcounter 
  # maxstopcounter <- pnn$maxstopcounter
  
  # parlist <- pnn$parlist
  # hidden_layers <- pnn$hidden_layers
  # y = pnn$y
  # fe_var <- pnn$fe_var
  # lam <- pnn$lam
  # parapen <- pnn$parapen
  # treatment = pnn$treatment
  # mse <- pnn$mse
  
  constraint <- sum(c(parlist$beta_param*parapen, parlist$beta)^2)
  #getting implicit regressors depending on whether regression is panel
  if (is.null(treatment)){
    if (!is.null(fe_var)){
      Zdm <- demeanlist(hidden_layers[[length(hidden_layers)]], list(fe_var))
      targ <- demeanlist(y, list(fe_var))
    } else {
      Zdm <- hidden_layers[[length(hidden_layers)]]
      targ <- y
    }    
  } else {#HTE case with fe's
    if (!is.null(fe_var)){
      V <- hidden_layers[[length(hidden_layers)]]
      ints <- sweep(V[,grepl('nodes', colnames(V))], MARGIN = 1, STATS = treatment, "*")
      colnames(ints) <- paste0(colnames(ints),"_int")
      dmat <- data.frame(V, ints) 
      Zdm <- as.matrix(demeanlist(dmat, list(fe_var)))
      targ <- demeanlist(y, list(fe_var))      
    } else {#HTE case 
      V <- hidden_layers[[length(hidden_layers)]]
      ints <- sweep(V[,grepl('nodes', colnames(V))], MARGIN = 1, STATS = treatment, "*")
      colnames(ints) <- paste0(colnames(ints),"_int")
      dmat <- data.frame(V, ints) 
      targ <- y
    }
  }
  

  #setup_penalty_vector
  D <- rep(1, ncol(Zdm)) #D=numeric
  if (is.null(fe_var)){
    pp <- c(0, parapen) #never penalize the intercept
  } else {
    pp <- parapen #parapen
  }
  D[1:length(pp)] <- D[1:length(pp)]*pp 
  #incorporate_parapen_into_diagonal_of_covmat
  #set_to_zero -> treatment_columns_and_main_effects_of_HTE_models 
  #treatments and ME's are not penalized 
  #only_penalize_the_interacted_nodes

  if (!is.null(treatment)){
    tozero <- grepl('treatment', colnames(Zdm)) | (!grepl('nodes', colnames(Zdm)) & !grepl('int', colnames(Zdm)))
    D[tozero] <- 0
  }
  
  
  # if (lam != 0){
  #   #implicit_lambda_function
  #   f <- function(lam){
  #     bi <- glmnet(y = targ, x = Zdm, lambda = lam, alpha = 0, intercept = FALSE, penalty.factor = D, standardize = FALSE)
  #     bi <- as.matrix(coef(bi))[-1,]
  #     bi1 <- (crossprod(bi*D) - constraint)^2
  #     return(bi1)
  #   }
  
  if (lam != 0){
    #implicit_lambda_function
    f <- function(lam){
      #minimize _mse_given_constraint_over_betas_and_lam
      bi <- glmnet(y = targ, x = Zdm, lambda = lam, alpha = 0, intercept = FALSE, penalty.factor = D, standardize = FALSE)
      bi <- as.matrix(coef(bi))[-1,]
      bi <- (crossprod(bi*D) - constraint)^2
      bi <- matrix(rep(bi, (ncol(Zdm)*(nrow(Zdm)))), ncol(Zdm)) #turn_bi_to_scalar_for_conformibility
      bi <- (Zdm %*% bi)
      bi <- as.matrix(bi[1,])
      mean((targ - bi)^2) # targ_==_y
      #mse_1 <- mean((targ - Zdm %*% c(parlist$beta_param, parlist$beta))^2 
    }
    
    
    #optimize it
    o <- optim(par = lam, fn = f, method = 'Brent', lower = lam, upper = 1e9) #fn_=_function_to_be_minimized
    #new_lambda
    newlam <- o$par
  } else {
    newlam <- 0 # lam == 0
  }
  
  #New_top_level_params
  br <- glmnet(y = targ, x = Zdm, lambda = newlam, alpha = 0, intercept = FALSE, penalty.factor = D, standardize = FALSE)
  b <- as.matrix(coef(br))[-1,,drop = FALSE]
  parlist$beta_param <- b[grepl('param', rownames(b))]
  parlist$beta <- b[grepl('nodes', rownames(b))]
  return(parlist)
}

# pnn$mse
# mean((targ - Zdm %*% b)^2)

#mean((targ - Zdm %*% c(parlist$beta_param, parlist$beta))^2)
#mean((targ - Zdm %*% b)^2)




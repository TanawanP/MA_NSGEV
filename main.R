# packages
library(ismev)
library(SpatialExtremes)

# set working directory
setwd("/Users/thanawanp/Desktop/Model averaging using nonstationary extreme value models R code")

# source file
source("16NS_GEV.R")
source("RL_estimate.R")

# Model averaging functions
RL_CI_ENS_BM <- function(data, rp=100, set_t=20, w0="balance", condition=6, criterion="AIC"){
  # fit all 16 NS-GEV models
    model_fit <- fit_all(a=data) 
    sorted_mod <- model_fit[order(model_fit[,criterion], decreasing = FALSE), ]
    cri_min <- min(sorted_mod[criterion])
    inclu_mod_idx <- which((sorted_mod[criterion]-min(sorted_mod[criterion])) < condition)
    
    # additional checking
    inclu_mod <- best_ens_model(a.fit=sorted_mod, n.mod=length(inclu_mod_idx), cri=criterion)

    # set time series
    nt <- 1:(length(data)+set_t)
    
    # calculate BM Rl
    bm.fit <- sorted_mod[1,]
    if(bm.fit$Model=="GD000" || bm.fit$Model=="GEV000"){
      bm.rl <- get(paste0("rtlv.", bm.fit$Model))(x=bm.fit,rp=rp)
      bm.rl <- rep(bm.rl, last(nt))
    }else{
      bm.rl <- get(paste0("rtlv.", bm.fit$Model))(x=bm.fit,rp=rp,nt=nt)
    }
    
    # define weight for MA model
    if(w0=="balance"){
      w <- 1/nrow(inclu_mod)
    }else if(w0=="unbalance"){
      cri_val <- inclu_mod[criterion] %>% unlist(.)
      cri_val.diff <- cri_val - cri_min
      w <- exp(-cri_val.diff/2)/sum(exp(-cri_val.diff/2))
      names(w) <- inclu_mod$Model
    }
    
    # calculate MA RL
    ma.split_mod <- split(inclu_mod, as.factor(inclu_mod$Model))
    ma.rl_est <- lapply(X=ma.split_mod, nt=nt, rp=rp, FUN= get_rtlv)
    ma.rl.df <- do.call(cbind, ma.rl_est) 
    
    cal.match.name <- w[match(colnames(ma.rl.df), names(w))][col(ma.rl.df)]*ma.rl.df
    ma.rl <- apply(cal.match.name, MARGIN = 1, sum) 
    # output
    return(list(bm_name = bm.fit$Model,
                bm_fit = bm.fit,
                bm.rl = bm.rl, 
                ma.rl = ma.rl,
                ma.rl.df = ma.rl.df,
                ma_combined_model = inclu_mod,
                ma_weight = w,
                all_model_fit = model_fit
                ))
  }
  
# example data
Bangkok_tsr <- c(123.2,75.8,69.8,90.4,114.7, 93.3,124.2, 54.1 ,153.7,81.2, 98.4, 97.8,141.0, 84.0,
                 109.4,68.3,108.6, 52.9, 63.5,167.3, 84.1, 105.7,88.4, 99.1, 125.5,107.3,248.6,
                 156.7,142.9,96.8,143.9,123.1,126.2, 60.1, 185.9, 128.9, 93.5,78.4,128.4,114.5,
                 96.1,103.1,108.2,185.9,110.6,102.7,132.9,148.4, 70.1,216.8,115.8,157.4, 87.9,103.7,
                 97.2, 174.3, 141.5, 188.3, 110.9,  57.7, 118.0)

data <- Bangkok_tsr

# run function
A <- RL_CI_ENS_BM(data, rp=100, set_t=20, w0="unbalance", condition=6, criterion="AIC")

# the best model (BM) fit and return level (RL)
A$bm_name
A$bm_fit
A$bm.rl %>% plot()

# When the best fit model is nonsttaionay model, MA technique is recomended.
# the model averaging (MA) fit and RL
A$ma_combined_model
A$ma_weight

plot(A$bm.rl, type='l')
lines(A$ma.rl, col="red")
abline(v=length(data)) # present year of data

# if the combinded model in MA is too many, then the MA is not exist.
# condition = 4 is recomended
B <- RL_CI_ENS_BM(data, rp=100, set_t=20, w0="unbalance", condition=4, criterion="AIC")
B$ma.rl
B$ma_weight
B$ma.rl.df  

plot(B$bm.rl, type='l')
lines(B$ma.rl, col="red")
abline(v=length(data)) 



  

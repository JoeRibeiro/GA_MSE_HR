### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###
#stop('remember to uncomment each line if running for both one-way and random.')
#args_local <- c("n_workers=0", "n_iter=500", "yrs_hist=100", "yrs_proj=50", "fhist='one-way'", "stock_id=27", "OM=TRUE")

#it looks like i need to modify brps to change the recruitment model to 4 seasons.
# brps[[stock]] tells you that dimensions of season is 1, when it is 4 after running seasonalise
# 
# Can't run rosana's as I need the .rdata
# Can't run much of https://github.com/flr/doc/blob/9511edeec5f22333044fb1eaf55718b2804ec7d5/a_seasonal_operating_model.Rmd as it throws an error on the line fbar(eq)=refpts(eq)["msy","harvest"]%*%FLQuant....

# This is the FLasher I installed, the newer ones will not install, the next commit stops installation succeeding
#remotes::install_github("flr/FLasher@0358698f696b527a6d363edf383902ee2c8a0478")

setwd('C:\\Users\\JR13\\Documents\\LOCAL_NOT_ONEDRIVE\\GA_MSE_HR')

args_local <- c("n_workers=0", "n_iter=500", "yrs_hist=100", "yrs_proj=50", "fhist='random'", "stock_id=27", "OM=TRUE")

args <- commandArgs(TRUE)
if (exists(x = "args_local")) args <- append(args, args_local)
print("arguments passed on to this script:")
print(args)

### evaluate arguments passed to R
for (i in seq_along(args)) eval(parse(text = args[[i]]))

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
### use mse fork from shfischer/mse, branch mseDL2.0 
### remotes::install_github("shfischer/mse", ref = "mseDL2.0)
# remotes::install_github("flr/mse", ref = "2.2.1")
req_pckgs <- c("FLCore", "FLasher", "FLBRP", "mse", 
               "tidyr", "dplyr", "foreach", "doParallel")
for (i in req_pckgs) library(package = i, character.only = TRUE)

### load additional functions
source("funs.R")
source("season_functions.R")

### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###

if (isTRUE(n_workers > 1)) {
  ### start doParallel cluster
  cl <- makeCluster(n_workers)
  registerDoParallel(cl)
  cl_length <- length(cl)
  ### load packages and functions into parallel workers
  . <- foreach(i = seq(cl_length)) %dopar% {
    for (i in req_pckgs) library(package = i, character.only = TRUE,
                                 warn.conflicts = FALSE, verbose = FALSE,
                                 quietly = TRUE)
    source("funs.R", echo = FALSE)
  }
} else {
  cl <- FALSE
}

### ------------------------------------------------------------------------ ###
### fishing history dimensions ####
### ------------------------------------------------------------------------ ###

# n_iter <- 500
# yrs_hist <- 100
# yrs_proj <- 50

set.seed(2)

### ------------------------------------------------------------------------ ###
### with uniform distribution and random F trajectories ####
### ------------------------------------------------------------------------ ###
# fhist <- "random"#"one-way"
if (identical(fhist, "random")) {
  start <- rep(0, n_iter)
  middle <- runif(n = n_iter, min = 0, max = 1)
  end <- runif(n = n_iter, min = 0, max = 1)
  df <- t(sapply(seq(n_iter), 
                 function(x) {
                   c(approx(x = c(1, yrs_hist*4/2), 
                            y = c(start[x], middle[x]), 
                            n = yrs_hist*4/2)$y,
                     approx(x = c(yrs_hist*4/2, yrs_hist*4 + 1), 
                            y = c(middle[x], end[x]), 
                            n = (yrs_hist*4/2) + 1)$y[-1]) # added *4 here to correctly apply F across all the years / seasons otherwise it comes out in a weird saw-tooth (4 teeth across the 100 years)
                 }))
  
  f_array <- array(dim = c(yrs_hist*4, 3, n_iter),
                   dimnames = list(seq(yrs_hist*4), c("min","value","max"),
                                   iter = 1:n_iter))
  f_array[, "value", ] <- c(t(df))
}

### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###

### get lhist for stocks
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)

### BRPs from Fischer et al. (2020)
brps <- readRDS("input/brps.rds")

# # Changes to seasonalise a stock, here we modify the original brps... params
# brps[["sar"]]@params = FLPar(NA, dimnames=list(params=c("a","b"), season=1:4, iter=1))
# #S-R takes places in q4
# brps[["sar"]]@params[1,4]=51.2
# brps[["sar"]]@params[2,4]=90.9
# 
# # refpts
# original=brps[["sar"]]@refpts
# brps[["sar"]]@refpts = FLPar(NA, dimnames=list(refpt=c("virgin","msy","crash","f0.1","fmax"),quant=c("harvest","yield","rec","ssb","biomass","revenue","cost","profit"), season=1:4, iter=1))
# #S-R takes places in q4
# for(x in 1:7){
#   for(y in 1:5){
#     brps[["sar"]]@refpts[y,x] = original[y,x]
#   }
# }
# 
# # change fbar to have 4 dimensions in season
# replacement =FLQuant(dimnames=list(age="all",unit="unique",year=1:101,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=fbar(brps[["sar"]])
# fbar(brps[["sar"]]) = replacement
# 
# # change fbar.obs to have 4 dimensions in season
# replacement =FLQuant(dimnames=list(age="all",unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@fbar.obs
# brps[["sar"]]@fbar.obs = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age="all",unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@landings.obs
# brps[["sar"]]@landings.obs = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age="all",unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@discards.obs
# brps[["sar"]]@discards.obs = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age="all",unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@rec.obs
# brps[["sar"]]@rec.obs  = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age="all",unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@ssb.obs
# brps[["sar"]]@ssb.obs  = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age="all",unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@stock.obs
# brps[["sar"]]@stock.obs = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age="all",unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@profit.obs
# brps[["sar"]]@profit.obs = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@landings.sel
# brps[["sar"]]@landings.sel = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@discards.sel
# brps[["sar"]]@discards.sel = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@bycatch.harvest
# brps[["sar"]]@bycatch.harvest = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@stock.wt
# brps[["sar"]]@stock.wt = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@landings.wt
# brps[["sar"]]@landings.wt = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@discards.wt
# brps[["sar"]]@discards.wt = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@bycatch.wt
# brps[["sar"]]@bycatch.wt = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@m
# brps[["sar"]]@m = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@mat
# brps[["sar"]]@mat = replacement
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@harvest.spwn
# brps[["sar"]]@harvest.spwn = replacement
# 
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@m.spwn
# brps[["sar"]]@m.spwn = replacement
# 
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@availability
# brps[["sar"]]@availability = replacement
# 
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age=1:5,unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@price
# brps[["sar"]]@price = replacement
# 
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age="all",unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@vcost
# brps[["sar"]]@vcost = replacement
# 
# 
# # etc, etc...
# replacement =FLQuant(dimnames=list(age="all",unit="unique",year=1,season=1:4,area="unique",iter=1))
# replacement@.Data[,,,4,,]=brps[["sar"]]@fcost
# brps[["sar"]]@fcost = replacement
# 



### create FLStocks
stocks_subset <- stocks$stock[stock_id]#"bll"

if (exists("OM")) { 
  
  if (isTRUE(OM)) {
    
    stks_hist <- foreach(stock = stocks_subset, .errorhandling = "stop", 
                         .packages = c("FLCore", "FLasher", "FLBRP")) %dopar% {
      stk <- as(brps[[stock]], "FLStock")
      stk=FLCore::expand(stk,season=1:4)
      refpts <- refpts(brps[[stock]])
      stk <- qapply(stk, function(x) {#browser()
        dimnames(x)$year <- as.numeric(dimnames(x)$year) - 1; return(x)
      })
      stk <- stf(stk, yrs_hist + yrs_proj - dims(stk)$year + 1)
      stk <- propagate(stk, n_iter)
      stk=seasonalise(stk)
      ### create stock recruitment model
      stk_sr <- FLSR(params = params(brps[[stock]]), model = model(brps[[stock]]))
      stk_sr=FLCore::expand(stk_sr,season=1:4)
      
      ### create residuals for (historical) projection
      set.seed(0)
      residuals(stk_sr) <- rlnoise(dim(stk)[6], rec(stk) %=% 0, 
                                   sd = 0.6, b = 0)
      ### replicate residuals from GA paper
      set.seed(0)
      residuals(stk_sr)[, ac(0:150)] <- rlnoise(dim(stk)[6], 
                                                rec(stk)[, ac(0:150)] %=% 0,
                                                sd = 0.6, b = 0)
      ### replicate residuals from catch rule paper for historical period
      set.seed(0)
      residuals <- rlnoise(dim(stk)[6], (rec(stk) %=% 0)[, ac(1:100)], 
                           sd = 0.6, b = 0)
      residuals(stk_sr)[, ac(1:100)] <- residuals[, ac(1:100)]
      
      ### fishing history from previous paper
      if (isTRUE(fhist == "one-way")) {
        
        ### 0.5Fmsy until year 75, then increase to 0.8Fcrash
        fs <- rep(c(refpts["msy", "harvest"]) * 0.5, 74)
        f0 <- c(refpts["msy", "harvest"]) * 0.5
        fmax <- c(refpts["crash", "harvest"]) * 0.8
        rate <- exp((log(fmax) - log(f0)) / (25))
        fs <- c(fs, rate ^ (1:25) * f0)
        
        ### control object
        ctrl <- fwdControl(data.frame(year = 2:100, quantity = "f", val = fs))
        
      ### roller-coaster
      } else if (isTRUE(fhist == "roller-coaster")) {
        
        ### 0.5Fmsy until year 75, 
        ### increase to 0.8Fcrash in 10 years
        ### keep at 0.8Fcrash for 5 years
        ### reduce to Fmsy in last 5 years
        fs <- rep(c(refpts["msy", "harvest"]) * 0.5, 75)
        f0_up <- c(refpts["msy", "harvest"]) * 0.5
        fmax_up <- c(refpts["crash", "harvest"]) * 0.8
        yrs_up <- 15
        rate_up <- exp((log(fmax_up) - log(f0_up)) / yrs_up)
        yrs_down <- 6
        f0_down <- c(refpts["msy", "harvest"])
        rate_down <- exp((log(fmax_up) - log(f0_down)) / yrs_down)
        fs <- c(fs, rate_up ^ seq(yrs_up) * f0_up, rep(fmax_up, 3),
                rev(rate_down ^ seq(yrs_down) * f0_down))
        
        ### control object
        ctrl <- fwdControl(data.frame(season=c(rep(seq(1:4),yrs_hist)), year = c(rep(2:100,each=4)), quantity = "f", val = fs))

      ### random F trajectories
      } else if (isTRUE(fhist == "random")) {
        
        ### control object template
        ctrl <- fwdControl(data.frame(season=c(rep(seq(1:4),yrs_hist)), year = c(rep(1:yrs_hist,each=4)), quantity = c("f"), val = NA))
        #ctrl <- fwdControl(data.frame(season=c(rep(1:4,each=100)), year = c(rep(seq(1:100),4)), quantity = c("f"), val = NA))
        
        ### add iterations
        ctrl@iters <- f_array
        ### target * Fcrash
        
        ctrl@iters <- ctrl@iters *  c(refpts["crash", "harvest"]) * 1
      }
      
      ### project fishing history
      stk_stf <- FLasher::fwd(stk,  control = ctrl, sr = stk_sr)#, sr.residuals = residuals(stk_sr),
                     #sr.residuals.mult = TRUE, maxF = 5) 
      #plot(stk_stf, iter = 1:50)
      #plot(ssb(stk_stf), iter = 1:50)
      ### run a few times to get closer to target
      # for (i in 1:5) {
      #   stk_stf <- fwd(stk_stf, ctrl, sr = stk_sr,
      #                  sr.residuals.mult = TRUE, maxF = 4)
      # }
      
      ### save OM files
      name(stk_stf) <- stock
      path <- paste0("input/", n_iter, "_", yrs_proj, "/OM/", fhist, 
                     "/", stock, "/")
      dir.create(path, recursive = TRUE)
      
      ### stock & sr (full history)
      saveRDS(stk_stf, file = paste0(path, "stk.rds"))
      saveRDS(stk_sr, file = paste0(path, "sr.rds"))
      
      return(NULL)
      #return(list(stk = stk_stf, sr = stk_sr))
    }
  }
}

# source("funs_OM.R")
# debugonce(input_mp)
# input <- input_mp(stock = "pol", fhist = "one-way", n_iter = 500, n_yrs = 50,
#                   MP = "hr")
# res <- do.call(mp, input)
# 
# input2 <- input_mp(stock = "pol", fhist = "one-way", n_iter = 500, n_yrs = 50,
#                   MP = "hr", interval = 2)
# res2 <- do.call(mp, input2)

# debugonce(input_mp)
# input_list <- input_mp(stocks = "pol", fhist = "one-way", n_iter = 500, n_yrs = 50,
#                        MP = "hr")
# res_ <- do.call(mp, input_list[[1]])

### ------------------------------------------------------------------------ ###
### gc() ####
### ------------------------------------------------------------------------ ###

gc()
if (!isFALSE(cl)) clusterEvalQ(cl, {gc()})

installed <- as.data.frame(installed.packages())
saveRDS(installed,'C:\\Users\\JR13\\Documents\\LOCAL_NOT_ONEDRIVE\\GA_MSE_HR\\environmentNotRENV.rds')
write.csv(installed,file='C:\\Users\\JR13\\Documents\\LOCAL_NOT_ONEDRIVE\\GA_MSE_HR\\environmentNotRENV.csv')
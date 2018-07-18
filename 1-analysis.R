###################################################
# Brandon P.M. Edwards & Adam C. Smith
# Model Analysis
# analysis.R
# Created July 2018
# Last Updated July 2018
###################################################

###################################################
# Clear Memory
###################################################

remove(list = ls())

###################################################
# Set Constants
###################################################

dir <- getwd()
models <- c("FD", "Standard", "GAM", "GAMYE")

###################################################
# Read Data
###################################################

bugs <- read.csv("input/bugs_data.csv")

###################################################
# Cross Validation LPPD
###################################################

loocv <- NULL
for (model in models)
{
  print(paste("Loocv", model))
  # Build data frame of posteriors for each year
  # For Barn Swallow, since we saved 2000 steps of the posterior, we will end up with
  #   a data frame that is 100762 columns (representing each count for BARS) by 
  #   2000 rows (representing each saved step)
  post <- NULL
  bugs_index <- NULL
  for (year in c(1:51))
  {
    load(paste(dir, "/input/cv/", model, "/year", year, "removed.Rdata", sep = ""))
    if (year == 1)
    {
      post <- data.frame(jagsjob$sims.list$LambdaSubset)
    }
    else
    {
      post <- cbind(post, data.frame(jagsjob$sims.list$LambdaSubset))
    }
    
    if (model == "GAM")
    {
      year_data <- jagsjob$model$cluster1$data()
    }else{
      year_data <- jagsjob$model$data()
    }
    
    bugs_index <- c(bugs_index, year_data$I)
    
    remove(jagsjob)
  }
  #save(post, file = paste(model,"_posterior.Rdata", sep = ""))
  
  # create data frame to store all loocv values for each count and all models
  if(model == models[1]){
    loocv_mod.m <- data.frame(index = bugs_index)
  }
  
  #create empty vector to store all loocv_mod.v values for this model
  loocv_mod.v <- vector(length = ncol(post))
  # LOOCV is the log of the mean of the probability of the true count given each
  #  saved lambda (2000 per count), summed across all counts (101762 counts, i.e. ncol(post))
  for (i in 1:ncol(post))
  {
    loocv_mod.v[i] <- log(mean(dpois(bugs[bugs_index[i],"count"], post[,i])))
  }
  loocv_mod <- sum(loocv_mod.v)
  loocv <- c(loocv, loocv_mod)
  loocv_mod.m[,paste(model,"_loocv", sep = "")] <- loocv_mod.v #add vector to dataframe
}

save(loocv, file = "loocv.Rdata")
#save(loocv_mod.m, file = "loocv_pointwise.Rdata")

bugs_loocv <- merge(bugs, loocv_mod.m, by.x = "X", by.y = "index", all = T)
write.csv(bugs_loocv, file = "bugs_loocv.csv")

# Remove 99.9th percentile of data points to account for weirdness
q <- NULL
loocv_trim <- NULL
for (model in models)
{
  q99 <- quantile(bugs_loocv[,model], c(.001), na.rm = TRUE)
  q <- c(q, q99)
  trim <- bugs_loocv[which(bugs_loocv[,model] >= q99),]
  loocv_trim <- c(loocv_trim, sum(trim[,model], na.rm = TRUE))
}

###################################################
# WAIC
###################################################

WAIC <- NULL
pWAIC_point <- data.frame(index = bugs$X)
lppd_point <- data.frame(index = bugs$X)

for (model in models)
{
  print(paste("wAIC", model))
  WAIC_mod <- 0
  load(paste(dir, "/input/waic/", model, ".Rdata", sep = ""))
  
  lppd_v <- vector(length = nrow(bugs))
  pwaic_v <- vector(length = nrow(bugs))
  
  for(i in 1:nrow(bugs))
  {
    lppd_v[i] <- log(mean(dpois(bugs[i,"count"], jagsjob$sims.list$lambda[,i])))
    pwaic_v[i] <- var(log(dpois(bugs[i,"count"], jagsjob$sims.list$lambda[,i])))
  }
  WAIC_mod <- -2 * (sum(lppd_v) - sum(pwaic_v))
  WAIC <- c(WAIC, WAIC_mod)
  
  lppd_point[,paste(model, "_lppd", sep = "")] <- lppd_v
  pWAIC_point[,paste(model, "_pwaic", sep = "")] <- pwaic_v
  remove(jagsjob)
}

bugs_waic <- merge(bugs, lppd_point, by.x = "X", by.y = "index")
bugs_waic <- merge(bugs_waic, pWAIC_point, by.x = "X", by.y = "index")

write.csv(bugs_waic, file = "bugs_waic.csv")
save(WAIC, file = "WAIC.Rdata")

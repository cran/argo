library(xts)
library(glmnet)
load("data/test_dt1.Rdata")
set.seed(48425959)

ar_regtype <- list()
regtype <- "Lasso"
alpha <- 1

logit <- function(x) log(x) - log(1-x)

set.seed(48425959)
argo_orig <-
  ar_lasso_exogen(data = do.call("logit", list(GFT_rescale$CDC.data))[1:170,],
                  exogen = do.call("log", list(GFT_origscale[,-1]+0.5))[1:170,],
                  N_lag = 1:52, alpha=alpha, use_all_previous=use_all_previous,
                  N_training=N_training)
set.seed(48425959)
argo_pkg <-
  argo(data = do.call("logit", list(GFT_rescale$CDC.data))[1:170,],
       exogen = do.call("log", list(GFT_origscale[,-1]+0.5))[1:170,],
       N_lag = 1:52, alpha=alpha, use_all_previous=use_all_previous,
       N_training=N_training)

expect_equal(argo_orig$pred, argo_pkg$pred)

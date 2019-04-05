library(xts)
library(glmnet)
load("data/gc09.RData")
set.seed(48425959)
orig.lasso1 <- ar_lasso_exogen(data = do.call(linkfun, list(GFT_rescale$CDC.data))[1:130,],
                               exogen = do.call(linkfun, list(google_data))[1:130,],
                               N_lag = NULL, alpha=1, use_all_previous=TRUE,
                               N_training=N_training)
set.seed(48425959)
orig.lasso2 <- argo::argo(data = do.call(linkfun, list(GFT_rescale$CDC.data))[1:130,],
                          exogen = do.call(linkfun, list(google_data))[1:130,],
                          N_lag = NULL, alpha=1, N_training=104, use_all_previous=TRUE)

expect_equal(orig.lasso1$pred, orig.lasso2$pred)

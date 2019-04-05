library(xts)
library(glmnet)
load("data/test_dt1.Rdata")
set.seed(48425959)

ar_3 <-  list()
linkfun <- "identity"
regtype <- "AR3"
ar_3_orig <-
  AR_model(data = do.call(linkfun, list(GFT_rescale$CDC.data)), use_all_previous=use_all_previous,
           N_training=N_training)


regtype <- "AR3+GFT"
ar_3_gft_orig <-
  AR_model(data = do.call(linkfun, list(GFT_rescale$CDC.data)), use_all_previous=use_all_previous,
           N_training=N_training, exogen = GFT_xts$updated_gft/100)


ar3_argo <- argo(data = do.call(linkfun, list(GFT_rescale$CDC.data)), use_all_previous=use_all_previous,
                 N_training=N_training, alpha = NA, N_lag = 1:3)
ar3_gft_argo <- argo(data = do.call(linkfun, list(GFT_rescale$CDC.data)), use_all_previous=use_all_previous,
                 N_training=N_training, alpha = NA, N_lag = 1:3,exogen = GFT_xts$updated_gft/100)

expect_equal(ar_3_orig$pred, ar3_argo$pred)
expect_equal(ar_3_gft_orig$pred, ar3_gft_argo$pred)

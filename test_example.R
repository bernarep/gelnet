load("MC_my_data_test_short.RData")
source("bic_gelnet_rev_comp.R")
source("CR_gelasticnet_mse.R")
outtest1 = bic_gelnet(X = full_data[,,1,1], ebic_par = 0, algorithm=1, alpha_seq = 0.5, lambda_seq = seq(0.01,1,length.out=10), toler = 1e-5, maxit = 500)
outtest2 = bic_gelnet(X = full_data[,,1,1], ebic_par = 0, algorithm=2, alpha_seq = 0.5, lambda_seq = seq(0.01,1,length.out=10), toler = 1e-5, maxit = 500)
outtest3 = CRgelnetCV_fixed(X = full_data[,,1,1], kcv=5, lambdas=seq(0.01,1,length.out=10), alphas=0.5, thr=1e-5, maxit=500, sym_mtd = 1)
outtest1$optimal_model$opt_O
outtest2$optimal_model$opt_O
outtest3$optimal_model$opt_O


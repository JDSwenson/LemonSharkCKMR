rm(list=ls())
set.seed(20)

source("functions/Data_format.R")
source("functions/get_P_lemon.R")
source("functions/lemon_neg_log_like.R")

lemon_neg_log_lik(Pars=Pars,Negatives_Mother=Data_mom_no,Negatives_Father=Data_dad_no,Pairs_Mother=Data_mom_yes,Pairs_Father=Data_dad_yes,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

library(optimx)
CK_fit <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE,Negatives_Mother=Data_mom_no,Negatives_Father=Data_dad_no,Pairs_Mother=Data_mom_yes,Pairs_Father=Data_dad_yes,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

head(Samples)
Samples$Indiv_ID[which(Samples$Indiv_ID %in% Samples$Father)]
Samples$Indiv_ID[which(Samples$Indiv_ID %in% Samples$Mother)]

summary(CK_fit)
#compute variance covariance matrix
#Put this on hold - try to understand
D=diag(length(Pars))*c(exp(CK_fit$p1[2]),exp(CK_fit$p2[2])) #derivatives of transformations
VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
VC = (t(D)%*%VC_trans%*%D) #delta method
SE=sqrt(diag(VC))

cat(paste("Estimated mature female abundance: ",exp(CK_fit$p1),", SE = ",SE[1],"\n"))
cat(paste("Estimated mature male abundance: ",exp(CK_fit$p2),", SE = ",SE[2],"\n"))

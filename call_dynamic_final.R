rm(list=ls())
load("/Users/qinghe/Documents/project/zigp/ZINB/simulation/simdata_dynamic.RData")
source("~/Documents/project/zigp/ZINB/ZINB_NNGP_final.R")
source("~/Documents/project/zigp/ZINB/plots_ZINB_NNGP_MV_V3.R")
nsim = 20000
burn=2000
thin=5
M=19
save_ypred=TRUE
results = ZINB_NNGP(X=X, x=x, y=y, coords=coords, nis=NULL, Vs=Vs, Vt=Vt, Ds=Ds, 
                    Dt=Dt, M=M,nsim=nsim, burn=burn, thin=thin,save_ypred=save_ypred)
save_location ="/Users/qinghe/Documents/project/zigp/ZINB/simulation/output/" 
plot_ZINB_NNGP(results=results, coords=coords, 
               true.y = true.y, true.alpha=true.alpha, true.beta=true.beta,
               true.a=true.a,true.b=true.b,true.c=true.c,true.d=true.d,
               true.eps1s=true.eps1s,true.eps2s=true.eps2s,
               true.eps1t=true.eps1t,true.eps2t=true.eps2t,
               true.l1t=true.l1t,true.sigma1t=true.sigma1t,true.l2t=true.l2t,true.sigma2t=true.sigma2t,
               true.phi_bin=true.phi_bin,true.sigma1s=true.sigma1s,
               true.phi_nb=true.phi_nb,true.sigma2s=true.sigma2s,
               true.sigma_eps1s=true.sigma_eps1s,true.sigma_eps2s=true.sigma_eps2s,
               true.sigma_eps1t=true.sigma_eps1t,true.sigma_eps2t=true.sigma_eps2t,
               save_location=save_location,save_file_suffix = "sim3_1_dynamic_M1_V5_testalg2")

source("/Users/qinghe/Documents/project/zigp/ZINB/heatmap_V2.R")
save_location ="/Users/qinghe/Documents/project/zigp/ZINB/simulation/output/" 
save_file_suffix="sim_dynamic"
plotSpatialSitesEffect(coords, results, true.a=true.a, true.c=true.c,save_location,save_file_suffix)
plotTemporalEffect2(results=results, true.b=true.b, true.d=true.d,
                    save_location=save_location,save_file_suffix=save_file_suffix)
  

mY_tp=tapply(y, tp_seq, mean)

#############################################
mY_pred<-colMeans(results$Y_pred[3000:5000,])
mY_pred_tp=tapply(mY_pred, tp_seq, mean)
test_data <-
  data.frame(
    true = mY_tp,
    predicted = mY_pred_tp,
    time = 1:ncol(Vt))
test_data %>%
  gather(key,value, true, predicted) %>%
  dplyr::rename(variable=key) %>% 
  ggplot(aes(x=time, y=value, colour=variable)) +
  geom_line()+
  geom_point()+
  ggtitle("Time series plot of the response") +
  theme(plot.title = element_text(hjust = 0.5,size=9))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())

############################################################
Y_pred=results$Y_pred
Y_pred_tp=matrix(nrow=nrow(Y_pred),ncol=n_time_points)
for (i in 1:nrow(Y_pred_tp)) {
  Y_pred_tp[i,]=tapply(Y_pred[i,], tp_seq, mean)
}
Y_pred_tp_mean<-colMeans(Y_pred_tp)
Y_pred_tp_ci<-apply(Y_pred_tp,2,quantile,prob=c(0.025,0.975))
test_data2 <-
  data.frame(
    true = mY_tp,
    predicted = Y_pred_tp_mean,
    time = 1:ncol(Vt),
    y_pred_low = Y_pred_tp_ci[1,],
    y_pred_high = Y_pred_tp_ci[2,])
pdf(file = "/Users/qinghe/Documents/project/zigp/ZINB/simulation/output/y_pred_sim_dynamic.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4)
ggplot(test_data2, aes(x=time)) + 
  geom_line(aes(y = true), color = "darkred") + 
  geom_line(aes(y = predicted), color="steelblue", linetype="twodash") +
  geom_ribbon(aes(ymin = y_pred_low, ymax = y_pred_high), alpha = 0.2) +
  scale_y_continuous(limits = c(0, 3))
dev.off()
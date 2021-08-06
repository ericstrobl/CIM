
### RUN EXPERIMENTS

# define number of samples
samps = 2000

# instantiate lists to record outputs
res=lapply(1:nIter,function(x) vector("list", length = length(samps)))
res_graphs=vector("list", length = nIter)

for (t in 1:50){
  print(t)
  waves1 = list(w1=1:8,w2=9:16,w3=17:24)
  
  
  for (ss in seq_len(length(samps))){
    
    # generate a mixture of DAGs
    mixDAG = generate_mix_DAGs2(24,en=2,waves1)
    
    # sample from the mixture of DAGs
    synth_data = sample_mix_DAGs2(mixDAG,samps)
    res_graphs[[t]] = mixDAG
    
    # randomize the variable order
    resort_p = sample(c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3),24-length(mixDAG$L),replace=FALSE)
    waves = list(w1 = match(mixDAG$waves_L$w1,resort_p), w2 = match(mixDAG$waves_L$w2,resort_p),
                 w3 = match(mixDAG$waves_L$w3,resort_p))
    prior_know = NULL
    
    suffStat$data = synth_data[,resort_p];
    
    # perform skeleton discovery
    resSkel = skel_discovery_waves3(suffStat, GCM_wrap, alpha=0.01,
                                    p=ncol(suffStat$data), waves) 
    
    # run CIM
    start_time <- proc.time()
    cim_out = CIM(suffStat, GCM_wrap, alpha=0.01, p=ncol(suffStat$data), 
                  verbose=FALSE, waves=waves, prior_know = prior_know, resIP=resSkel)
    cim_out$G = cim_out$pofaag[order(resort_p),order(resort_p)]
    colnames(cim_out$G) <- c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3)
    rownames(cim_out$G) <- c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3)
    res[[t]]$cim_waves = cim_out$G
    res[[t]]$cim_waves_time = proc.time()-start_time+resSkel$time_skel;
    
    # perform IP discovery
    resIP = IP_discovery_waves3(suffStat, resSkel$skel, GCM_wrap, alpha=0.01,
                                p=ncol(suffStat$data), waves)
    
    # run PC
    start_time <- proc.time()
    pc_out = pc_pre(suffStat, GCM_wrap, alpha=0.01, p=ncol(suffStat$data), 
                    verbose=FALSE, skeleton_pre=resSkel,waves=waves,
                    prior_know=prior_know)
    pc_out$G = pc_out$maag[order(resort_p),order(resort_p)]
    colnames(pc_out$G) <- c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3)
    rownames(pc_out$G) <- c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3)
    res[[t]]$pc = pc_out$G
    res[[t]]$pc_time = proc.time()-start_time+resSkel$time_skel;
    
    # run FCI
    start_time <- proc.time()
    fci_out = fci_pre(suffStat, GCM_wrap, alpha=0.01, p=ncol(suffStat$data),
                      verbose=FALSE, skeleton_pre=resIP)
    fci_out$G = fci_out$maag[order(resort_p),order(resort_p)]
    colnames(fci_out$G) <- c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3)
    rownames(fci_out$G) <- c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3)
    res[[t]]$fci = fci_out$G
    res[[t]]$fci_time = proc.time()-start_time+resIP$time_pdsep+resSkel$time_skel;
    
    # run RFCI
    start_time <- proc.time()
    rfci_out = rfci_pre(suffStat, GCM_wrap, alpha=0.01, p=ncol(suffStat$data),
                        verbose=FALSE, skeleton_pre=resSkel)
    rfci_out$G = rfci_out$maag[order(resort_p),order(resort_p)]
    colnames(rfci_out$G) <- c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3)
    rownames(rfci_out$G) <- c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3)
    res[[t]]$rfci = rfci_out$G
    res[[t]]$rfci_time = proc.time()-start_time+resSkel$time_skel;
    
    # run CCI
    start_time <- proc.time()
    cci_out = cci_pre(suffStat, GCM_wrap, alpha=0.01, p=ncol(suffStat$data),
                      verbose=FALSE, skeleton_pre=resIP)
    cci_out$G = cci_out$maag[order(resort_p),order(resort_p)]
    colnames(cci_out$G) <- c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3)
    rownames(cci_out$G) <- c(mixDAG$waves_L$w1,mixDAG$waves_L$w2,mixDAG$waves_L$w3)
    res[[t]]$cci = cci_out$G
    res[[t]]$cci_time = proc.time()-start_time+resIP$time_pdsep+resSkel$time_skel;
    
    res[[t]][[ss]]$resIP = resIP
  }
  
  # save outputs
  save(res, res_graphs, file = "synthetic_results.RData")
}


### RESULTS

metrics_fo = matrix(0,50,5) # fallout
metrics_OA = matrix(0,50,5) # overall score
metrics_sens = matrix(0,50,5) # sensitivity
metrics_time = matrix(0,50,5) # time
algs_index = c(2,4,6,8,10)
for (t in 1:50){
  load("synthetic_results.RData")
  for (a in algs_index){
    
    # compute fallout
    cum_res2 = res_graphs[[t]]$graph[res_graphs[[t]]$waves_L$w2,res_graphs[[t]]$waves_L$w2]
    ans_res2 = res[[t]][[a]][(length(res_graphs[[t]]$waves_L$w1)+1):(length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2)),
                             (length(res_graphs[[t]]$waves_L$w1)+1):(length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2))]
    err2 = length(intersect(which(t(cum_res2)!=0),which(t(ans_res2)==3)))
    
    cum_res3 = res_graphs[[t]]$graph[res_graphs[[t]]$waves_L$w3,res_graphs[[t]]$waves_L$w3]
    ans_res3 = res[[t]][[a]][(length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2)+1):(length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2)+length(res_graphs[[t]]$waves_L$w3)),
                             (length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2)+1):(length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2)+length(res_graphs[[t]]$waves_L$w3))]
    err3 = length(intersect(which(t(cum_res3)!=0),which(t(ans_res3)==3)))
    
    cres2 = sum(t(cum_res2)!=0)
    cres3 = sum(t(cum_res3)!=0)
    if ((cres2+cres3) > 0){
      errF = (err2+err3)/(cres2+cres3)
    } else{
      errF = 0;
    }
    metrics_fo[t,which(algs_index==a)]=errF # record fallout
  
   
    # compute sensitivity
    cum_res2 = res_graphs[[t]]$graph[res_graphs[[t]]$waves_L$w2,res_graphs[[t]]$waves_L$w2]
    ans_res2 = res[[t]][[a]][(length(res_graphs[[t]]$waves_L$w1)+1):(length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2)),
                             (length(res_graphs[[t]]$waves_L$w1)+1):(length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2))]
    err2 = length(intersect(which(cum_res2!=0),which(t(ans_res2)==3)))
    
    cum_res3 = res_graphs[[t]]$graph[res_graphs[[t]]$waves_L$w3,res_graphs[[t]]$waves_L$w3]
    ans_res3 = res[[t]][[a]][(length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2)+1):(length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2)+length(res_graphs[[t]]$waves_L$w3)),
                             (length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2)+1):(length(res_graphs[[t]]$waves_L$w1)+length(res_graphs[[t]]$waves_L$w2)+length(res_graphs[[t]]$waves_L$w3))]
    err3 = length(intersect(which(cum_res3!=0),which(t(ans_res3)==3)))
    
    cres2 = sum(cum_res2!=0)
    cres3 = sum(cum_res3!=0)
    if ((cres2+cres3) > 0){
      err = (err2+err3)/(cres2+cres3)
    } else{
      err = 0;
    }
    metrics_sens[t,which(algs_index==a)]=err # record sensitivty
    metrics_OA[t,which(algs_index==a)]=sqrt( metrics_fo[t,which(algs_index==a)]^2 + (1-metrics_sens[t,which(algs_index==a)])^2 )  # record overall score
    metrics_time[t,which(algs_index==a)]=res[[t]][[a+1]][3] # record time
  }
}

print(colMeans(metrics_sens))
print(colMeans(metrics_fo))
print(colMeans(metrics_OA))

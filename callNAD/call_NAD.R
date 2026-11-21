setwd('/lustre/user/user/')
library(depmixS4)

run_callNAD_genome = function(cell = char()){
  require(data.table)
  ### call NAD with depmixS4 HMM with gaussion log distribution
  call_NAD_HMM2 = function(mat,quantiles = c(0.05,0.95),
                           input = 'log2R'){
    require(depmixS4)
    require(cowplot)
    require(dplyr)
    require(ggplot2)
    ### quantile for outlier normalization
    mat = as.data.frame(mat)
    qupper = quantile(mat[,input],max(quantiles))
    qlower = quantile(mat[,input],min(quantiles))
    mat$log2Rq = mat[,input]
    
    mat$log2Rq[mat$log2Rq>qupper] = qupper
    mat$log2Rq[mat$log2Rq<qlower] = qlower
    
    ###get log2Ratio of NS/WGS seq from bedgraph(deepTools)
    ###train HMM based on gaussian(log) distribution (log normal distribution)
    set.seed(923)
    options(scipen = 4)
    
    mod = depmix(response = log2Rq~1,data = mat,nstates = 2,family = gaussian(link = 'log'))
    fm = fit(mod)
    summary(fm,which = 'prior')
    
    k = posterior(fm,type='viterbi')
    states = posterior(fm,type='global')
    posteriors = data.frame(posterior(fm,type='smoothing'))
    posteriors = cbind(posteriors,states)
    posteriors$p = apply(posteriors,1,function(x){
      p = x[x[3]]
      return(p)
    })
    
    mat$state = posteriors$state
    mat$p =  posteriors$p
    st = unique(mat$state)
    
    t = wilcox.test(x=mat$log2Rq[mat$state==st[1]],y=mat$log2Rq[mat$state==st[2]],alternative = 'greater')
    if(t$p.value<0.05){
      mat$state_new[mat$state==st[1]]='NAD'
      mat$state_new[mat$state==st[2]]='iNAD'
    }else{
      mat$state_new[mat$state==st[1]]='iNAD'
      mat$state_new[mat$state==st[2]]='NAD'
    }
    
    
    cols = c('red','grey50');names(cols) = c('NAD','iNAD')
    ratio = round(table(mat$state_new)['NAD']/sum(table(mat$state_new)),3)
    p = ggplot(mat, aes(x=log2R, fill=state_new)) +
      geom_histogram(position="identity", alpha=0.5, bins=300) + 
      scale_fill_manual(values = cols) +
      labs(title = paste0("Distribution of log2R by Hidden States\nRatio~",ratio), x = "log2R", y = "Frequency") +
      theme_cowplot()
    print(p)
    # return(mat)
    
    ### filter out NAD regions which mean HMM significance < threshold
    tmp = mat
    to = c(which(diff(tmp$state)!=0),dim(tmp)[1])
    from = c(1,which(diff(tmp$state)!=0)+1)
    tmp$frag = apply(data.frame(from,to,paste0('frag',1:length(from))),1, function(x){
      n = as.numeric(x[2])-as.numeric(x[1])+1
      x = rep(as.character(x[3]),n)
      return(x)
    }) %>% 
      unlist() %>% 
      factor(.,levels = paste0('frag',1:length(from)))
    # tmp = tmp[tmp$state_new=='NAD',]
    tmp = tmp %>% dplyr::group_by(frag) %>% 
      dplyr::mutate(mean_pval = median(p),
                    median_logR = median(log2R),
                    pos_cov = length(log2R[log2R>0])/length(log2R))
    mat = merge(mat,tmp[,c('chr','start','end','mean_pval','median_logR','pos_cov')],by = c('chr','start','end'),all.x=T, all.y=F)
    return(mat)
  }
  call_NAD_HMM3 = function(mat,quantiles = c(0.05,0.95),
                           input = 'log2R'){
    require(depmixS4)
    require(cowplot)
    require(dplyr)
    require(ggplot2)
    ### quantile for outlier normalization
    mat = as.data.frame(mat)
    qupper = quantile(mat[,input],max(quantiles))
    qlower = quantile(mat[,input],min(quantiles))
    mat$log2Rq = mat[,input]
    
    mat_pos = mat[mat[,input]>0,]
    hist(mat_pos[,input],breaks = 500)
    set.seed(923)
    options(scipen = 4)
    mod = depmix(response = log2Rq~1,data = mat,nstates = 3,family = gaussian(link = 'log'))
    fm = fit(mod)
    summary(fm,which = 'prior')
    k = posterior(fm,type='viterbi')
    states = posterior(fm,type='global')
    posteriors = data.frame(posterior(fm,type='smoothing'))
    posteriors = cbind(posteriors,states)
    posteriors$p = apply(posteriors,1,function(x){
      p = x[x[4]]
      return(p)
    })
    
    mat$state = posteriors$state
    mat$p =  posteriors$p
    mat$state = as.numeric(mat$state)
    st = unique(mat$state)
    median = c(median(mat$log2R[mat$state==st[1]]),
               median(mat$log2R[mat$state==st[2]]),
               median(mat$log2R[mat$state==st[3]]))
    names(median) = st
    median = sort(median)
    mat$state_new = ifelse(mat$state==names(median)[1],'iNAD',
                           ifelse(mat$state==names(median)[3],'NAD','rNAD'))
    
    
    cols = c('red','blue','grey50');names(cols) = c('NAD','rNAD','iNAD')
    ratio = round(table(mat$state_new)[c('NAD','rNAD')]/sum(table(mat$state_new)),3)
    p = ggplot(mat, aes(x=log2R, fill=state,group = state)) +
      geom_histogram(position="identity", alpha=0.5, bins=300) + 
      # scale_fill_manual(values = cols) +
      labs(title = paste0("Distribution of log2R by Hidden States\n Ratio~",ratio), x = "log2R", y = "Frequency") +
      theme_cowplot()
    print(p)
    
    ### filter out NAD regions which mean HMM significance < threshold
    tmp = mat
    to = c(which(diff(tmp$state)!=0),dim(tmp)[1])
    from = c(1,which(diff(tmp$state)!=0)+1)
    tmp$frag = apply(data.frame(from,to,paste0('frag',1:length(from))),1, function(x){
      n = as.numeric(x[2])-as.numeric(x[1])+1
      x = rep(as.character(x[3]),n)
      return(x)
    }) %>% 
      unlist() %>% 
      factor(.,levels = paste0('frag',1:length(from)))
    tmp = tmp %>% dplyr::group_by(frag) %>% 
      dplyr::mutate(mean_pval = median(p),
                    median_logR = median(log2R),
                    pos_cov = length(log2R[log2R>0])/length(log2R))
    mat = merge(mat,tmp[,c('chr','start','end','mean_pval','median_logR','pos_cov')],by = c('chr','start','end'),all.x=T, all.y=F)
    return(mat)
  }
  data.path = list.files('bedgraph',pattern = paste0(cell,'_rmBL'),full.names = T)
  dir.create('output',recursive = T)
  
  data = fread(data.path,header = F,col.names = c('chr','start','end','log2R'))
  data = data[data$chr%in%paste0('chr',c(1:22,'X')),]
  
  ### Hanning smoothing in 21kb window
  window_size <- 21  
  data$smoothed <- zoo::rollapply(data$log2R, width = window_size, FUN = function(x) {
    weights <- 0.5 * (1 - cos(2 * pi * seq_along(x) / (length(x) + 1)))
    sum(x * weights) / sum(weights)
  }, fill = NA, align = "center")
  data$smoothed[is.na(data$smoothed)]=0
  
  bdg = data[,c('chr','start','end','log2R')]
  write.table(bdg,paste0('output/',cell,'_NAD.bdg'),quote=F,sep = '\t',row.names = F,col.names = F)
  
  
  nad = call_NAD_HMM2(data,quantiles = c(0.05,0.95),input = 'smoothed')
  raw_NAD = nad[nad$state_new=='NAD',c('chr','start','end')]
  iNAD = nad[nad$state_new=='iNAD'&nad$mean_pval>=0.95,c('chr','start','end')]
  ### double check of identified NAD, remove low signal with low significance regions.
  NAD_score = quantile(nad$median_logR[nad$state_new=='NAD'],0.25,na.rm=T)
  filt_idx = which(nad$state_new=='NAD' & (nad$mean_pval>=0.95|(nad$median_logR>NAD_score&nad$pos_cov>=0.8)))
  filter_NAD = nad[filt_idx,c('chr','start','end')]
  
  write.table(raw_NAD,paste0('output/',cell,'_NAD.bed'),quote=F,sep = '\t',row.names = F,col.names = F)
  write.table(filter_NAD,paste0('output/',cell,'_NAD_filter.bed'),quote=F,sep = '\t',row.names = F,col.names = F)
  write.table(iNAD,paste0('output/',cell,'_iNAD.bed'),quote=F,sep = '\t',row.names = F,col.names = F)
}

sample = read_lines('sample.txt')
cl = makeCluster(length(sample))
parLapply(cl,sample,run_callNAD_genome)
stopCluster(cl)

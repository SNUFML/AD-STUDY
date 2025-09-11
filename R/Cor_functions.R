#normality test
NT<-function(data){
  outputs<-c()
  for(i in 1:ncol(data)){
    output<-c()
    st<-shapiro.test(data[,i])
    output<-c(st$statistic,st$p.value)
    outputs<-rbind(outputs,output)
  }
  colnames(outputs)<-c("shapiro test W","shapiro test p value")
  rownames(outputs)<-colnames(data)
  write.csv(outputs,paste("normality.test_",substitute(data),".csv",sep=""))
  
  
  pdf(paste("normality.test_",substitute(data),".pdf",sep=""))
  par(mfrow = c(3, 2))
  for(j in 1:ncol(data)){
    hist(data[,j],xlab  = colnames(data)[j],main=paste("Histogram of",colnames(data)[j],sep=" "))
    qqnorm(data[,j],main=paste("Normal Q-Q plot of",colnames(data)[j],sep=" ")) 
  }
  dev.off()
}
###part1

SC<-function(data1,data2){
  SCCmatrix<-matrix(nrow=ncol(data1),ncol=ncol(data2))
  pmatrix<-matrix(nrow=ncol(data1),ncol=ncol(data2))
  for(i in 1:ncol(data1)){
    for(j in 1:ncol(data2)){
      output<-c()
      scctest<-cor.test(data1[,i],data2[,j],method="spearman",exact=FALSE)
      pmatrix[i,j]<-scctest$p.value
      SCCmatrix[i,j]<-scctest$estimate
    }
  }
  colnames(SCCmatrix)<-colnames(data2)
  rownames(SCCmatrix)<-colnames(data1)
  colnames(pmatrix)<-colnames(data2)
  rownames(pmatrix)<-colnames(data1)
  return(list(SCCmatrix,pmatrix))
}

###part2
OR<-function(form,meta.d,clinic.d,outcome.d){
  outputs1<-c()#form1 OR and P
  for (i in 1:ncol(meta.d)){
    output1<-c()
    tempdata<-data.frame(meta.d[,i]/sd(meta.d[,i]),clinic.d,outcome.d)
    colnames(tempdata)<-c("meta",colnames(clinic.d),colnames(outcome.d))
    #model1
    cm1 <- try(clogit(form, data = tempdata),silent = T)
    if ('try-error' %in% class(cm1)) {
      output1<-rbind(outputs1,rep(0,4))
      next
    }else{
      CI<-exp(confint(cm1))[1,]
      OR<-exp(coef(cm1))[1]
      p<-as.data.frame(summary(cm1)$coefficient)$`Pr(>|z|)`[1]
      output1<-cbind(CI[1],CI[2],OR,p)
      outputs1<-rbind(outputs1,output1)
    }
  }
  fdr<-p.adjust(outputs1[,4],"fdr")
  outputs<-data.frame(outputs1,fdr)
  colnames(outputs)<-c("2.5% CI","97.5% CI","OR","p","fdr")
  rownames(outputs)<-colnames(meta.d)
  return(outputs)
}

#set training dataset and validation set
split_dataset<-function(n.seed,proportion.subdataset1){
  id<-table(outcome.d$match)
  set.seed(n.seed)
  train_row<-sample(1:length(id),length(id)*proportion.subdataset1)
  train_ID<-names(id)[train_row]
  #training dataset
  t_rows<-c()
  for(i in 1:length(train_ID)){
    t_row<-c()
    t_row<-which(outcome.d$match==train_ID[i])
    t_rows<-c(t_rows,t_row)
  }
  t_meta.d<-meta.d[t_rows,]
  t_clinic.d<-clinic.d[t_rows,]
  t_outcome.d<-outcome.d[t_rows,]
  #validation dataset
  v_meta.d<-meta.d[-t_rows,]
  v_clinic.d<-clinic.d[-t_rows,]
  v_outcome.d<-outcome.d[-t_rows,]
  output<-list(t_meta.d,t_clinic.d,t_outcome.d,v_meta.d,v_clinic.d,v_outcome.d)
  return(output)
}

#select validation dataset randomly
random_dataset<-function(n.seed,percentage){
  id<-table(outcome.d$match)
  set.seed(n.seed)
  select_row<-sample(1:length(id),length(id)*percentage)
  select_ID<-names(id)[select_row]
  #select dataset
  s_rows<-c()
  for(i in 1:length(select_ID)){
    s_row<-c()
    s_row<-which(outcome.d$match==select_ID[i])
    s_rows<-c(s_rows,s_row)
  }
  s_meta.d<-meta.d[s_rows,]
  s_clinic.d<-clinic.d[s_rows,]
  s_outcome.d<-outcome.d[s_rows,]
  output<-list(s_meta.d,s_clinic.d,s_outcome.d)
  return(output)
}
random_OR<-function(times.random,percentage){
  ORs<-c()
  ps<-c()
  fdrs<-c()
  for(i in 1:times.random){
    random.d<-random_dataset(i,percentage)
    r_ORresult<-OR(form1,random.d[[1]],random.d[[2]],random.d[[3]])
    OR<-r_ORresult$OR
    p<-r_ORresult$p
    fdr<-r_ORresult$fdr
    ORs<-cbind(ORs,OR)
    ps<-cbind(ps,p)
    fdrs<-cbind(fdrs,fdr)
  }
  ORs_g1_p<-c()
  ORs_l1_p<-c()
  ORs_g1_fdr<-c()
  ORs_l1_fdr<-c()
  for (i in 1:nrow(ORs)){
    OR_g1_p<-length(which(ORs[i,]>1 & ps[i,]<0.05))
    OR_l1_p<-length(which(ORs[i,]<1 & ps[i,]<0.05))
    OR_g1_fdr<-length(which(ORs[i,]>1 & fdrs[i,]<0.05))
    OR_l1_fdr<-length(which(ORs[i,]<1 & fdrs[i,]<0.05))
    ORs_g1_p<-c(ORs_g1_p,OR_g1_p)
    ORs_l1_p<-c(ORs_l1_p,OR_l1_p)
    ORs_g1_fdr<-c(ORs_g1_fdr,OR_g1_fdr)
    ORs_l1_fdr<-c(ORs_l1_fdr,OR_l1_fdr)
  }
  output<-data.frame(ORs_g1_p,ORs_g1_fdr,ORs_l1_p,ORs_l1_fdr)
  colnames(output)<-c("OR>1,p<0.05","OR>1,fdr<0.05","OR<1,p<0.05","OR<1,fdr<0.05")
  rownames(output)<-colnames(meta.d)
  return(output)
}

###Part3

#mediation function
CM<-function(form.a,form.b,form.c){
  outputs<-c()
  for(i in 1:ncol(meta.d)){
    for(j in 1:ncol(clinic.d)){
      output<-c()
      tempdata<-data.frame(outcome.d$case,(meta.d[,i]-mean(meta.d[,i]))/sd(meta.d[,i]),clinic.d[,j])
      colnames(tempdata)<-c("case","meta","clinic")
      cm.a<-lm(form.a,data=tempdata) ##ls as mediation meta-ls-t2dm
      cm.b<-glm(form.b,data=tempdata)##t2dm as outcome
      cm.c<-glm(form.c,data=tempdata)
      set.seed(001)
      a<-c()
      b<-c()
      c<-c()
      cm<-mediate(cm.a,cm.b,treat="meta",mediator="clinic",boot=T)
      a<-c(summary(cm.a)$coefficients[2,1],summary(cm.a)$coefficients[2,4])
      b<-c(summary(cm.b)$coefficients[3,1],summary(cm.b)$coefficients[3,4])
      c<-c(summary(cm.c)$coefficients[2,1],summary(cm.c)$coefficients[2,4])
      
      output<-c(colnames(meta.d)[i],colnames(clinic.d)[j],
                a,
                b,
                c,
                cm$d0,cm$d0.ci,cm$d0.p,
                cm$z0,cm$z0.ci,cm$z0.p,
                cm$tau.coef,cm$tau.ci,cm$tau.p,
                cm$n0,cm$n0.ci,cm$n0.p)
      outputs<-rbind(outputs,output)
    }
  }
  colnames(outputs)<-c("metabolites","clinical parameters",
                       "beta a","p a",
                       "beta b","p b",
                       "beta c","p c",
                       "acme","2.5%CI","97.5%CI","p acme",
                       "ade","2.5%CI","97.5%CI","p ade",
                       "te","2.5%CI","97.5%CI","p te",
                       "prop","2.5%CI","97.5%CI","p prop")
  return (outputs)
}

###Part4

VP<-function(form){
  corrs<-c()
  vars<-c()
  for(i in 1:ncol(meta.d)){
    tempdata<-data.frame(meta.d[,i],clinic.d)
    tempdata$sex<-factor(tempdata$sex)
    colnames(tempdata)<-c("meta",colnames(clinic.d))
    fit<-lm(form,tempdata)
    var<-as.matrix(calcVarPart(fit)*100)
    corr<-as.matrix(coef(fit))
    vars<-cbind(vars,var)
    corrs<-cbind(corrs,corr)
  }
  colnames(vars)<-colnames(meta.d)
  colnames(corrs)<-colnames(meta.d)
  vars<-vars[-nrow(vars),]
  corrs<-corrs[-1,]
  corr_type<-array(,c(nrow(corrs),ncol(corrs)))
  for(k in 1:nrow(corrs)){
    for(z in 1:ncol(corrs)){
      if(corrs[k,z]>=0 ){
        corr_type[k,z]<-"pos"
      }else {
        corr_type[k,z]<-"neg"
      }
    }
  }
  colnames(corr_type)<-colnames(meta.d)
  outputs<-list(vars,corr_type)
  return(outputs)
}

#plot
VPplot<-function(var,corr_type){
  plot_data1=melt(var[,1,drop=F])
  plot_data2=melt(corr_type[,1,drop=F])
  plot_data2$value=factor(plot_data2$value)
  
  ggplot(data=plot_data1,mapping=aes(x=value,y=Var1,pattern = plot_data2$value ))+
    geom_bar_pattern(stat="identity",color = "black",pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6)+
    scale_pattern_manual(values = c(neg = "stripe", pos = "none"))+
    labs(x="VPplot",y="")
  
}


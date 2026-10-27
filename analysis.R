prepare_data <- function(path) {
  d <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  needed <- c('center','report_group','sex','bmi','age','htn','dm','ckd',
              'procedure_flag','moderate_count','severe_count','additional_fx_yes',
              'med_90d','clinic_90d','followup_months')
  if (!all(needed %in% names(d))) stop('Missing required columns. See data_dictionary.csv.')
  d <- d[,needed]
  if (anyNA(d)) stop('Missing data: this pipeline reproduces the complete-data analysis, not imputation.')
  for (v in c('report_group','additional_fx_yes','med_90d','clinic_90d'))
    if (!all(d[[v]] %in% 0:1)) stop(paste('Invalid binary coding:',v))
  stopifnot(all(d$procedure_flag %in% 0:3),all(d$center %in% c('Center 1','Center 2')),
            all(d$sex %in% c('Female','Male')))
  for(v in c('age','bmi','moderate_count','severe_count','followup_months'))
    stopifnot(is.numeric(d[[v]]),all(is.finite(d[[v]])),all(d[[v]]>=0))
  for(v in c('moderate_count','severe_count')) stopifnot(all(d[[v]]==floor(d[[v]])))
  # Preserve original factor order for the bootstrap strata and model matrices.
  d$center <- factor(d$center,levels=c('Center 2','Center 1'))
  d$sex <- factor(d$sex,levels=c('Female','Male'))
  for(v in c('htn','dm','ckd')) {
    stopifnot(all(d[[v]] %in% c('No','Yes')))
    d[[v]] <- factor(d[[v]],levels=c('No','Yes'))
  }
  d$procedure <- factor(ifelse(d$procedure_flag %in% c(1,3),'Yes','No'),levels=c('No','Yes'))
  d$clinic_or_med_90d <- as.integer(d$med_90d==1 | d$clinic_90d==1)
  d
}

model_formula <- function(y,cohort='Pooled',adjustment='Reduced') {
  vars <- c('report_group','age','sex','bmi','moderate_count','severe_count',
            if(cohort=='Pooled')'center',
            if(adjustment=='Full')c('htn','dm','ckd'),
            if(adjustment!='Restricted')'procedure')
  reformulate(vars,y)
}

bind_rows <- function(x) {z<-do.call(rbind,x);rownames(z)<-NULL;z}
write_result <- function(x,name,out) write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE)

descriptives <- function(d,out) {
  groups<-list(All=seq_len(nrow(d)),Reported=which(d$report_group==1),Unreported=which(d$report_group==0))
  rows<-list()
  for(g in names(groups)) {
    z<-d[groups[[g]],,drop=FALSE]
    for(v in c('age','bmi','moderate_count','severe_count','followup_months'))
      rows[[paste(g,v)]]<-data.frame(group=g,variable=v,level='',n=nrow(z),mean=mean(z[[v]]),sd=sd(z[[v]]),count=NA,percent=NA)
    for(v in c('sex','htn','dm','ckd','center','procedure')) for(l in levels(z[[v]])) {
      count<-sum(z[[v]]==l)
      rows[[paste(g,v,l)]]<-data.frame(group=g,variable=v,level=l,n=nrow(z),mean=NA,sd=NA,count=count,percent=100*count/nrow(z))
    }
  }
  write_result(bind_rows(rows),'table1',out)
}

firth_models <- function(d,out) {
  results<-diagnostics<-cells<-counts<-list()
  for(adj in c('Reduced','Full','Restricted')) for(co in if(adj=='Restricted')'Pooled' else c('Pooled','Center 1','Center 2')) {
    z<-if(co=='Pooled')d else d[d$center==co,,drop=FALSE]
    if(adj=='Restricted') z<-z[z$procedure=='No',,drop=FALSE]
    for(y in c('additional_fx_yes','clinic_or_med_90d')) {
      fm<-model_formula(y,co,adj); key<-paste(adj,co,y)
      for(v in intersect(all.vars(fm),c('report_group','sex','center','procedure','htn','dm','ckd'))) {
        tab<-table(droplevels(factor(z[[v]])),factor(z[[y]],levels=0:1))
        cells[[paste(key,v)]]<-data.frame(model=key,variable=v,min_cell=min(tab),zero_cell=any(tab==0))
      }
      j<-which(colnames(model.matrix(fm,z))=='report_group'); warns<-character()
      f<-withCallingHandlers(logistf::logistf(fm,data=z,pl=TRUE,plconf=j),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
      stopifnot(all(is.finite(c(f$coefficients[j],f$ci.lower[j],f$ci.upper[j],f$prob[j]))))
      results[[key]]<-data.frame(cohort=co,outcome=y,adjustment_set=adj,n=nrow(z),events=sum(z[[y]]),odds_ratio=exp(f$coefficients[j]),ci_low=exp(f$ci.lower[j]),ci_high=exp(f$ci.upper[j]),p_value=f$prob[j])
      k<-length(f$coefficients)-1; epv<-min(sum(z[[y]]),nrow(z)-sum(z[[y]]))/k
      diagnostics[[key]]<-data.frame(model=key,parameters=k,EPV=epv,EPV_below_10=epv<10,warning=paste(warns,collapse=' | '),convergence=paste(f$conv,collapse=' | '),profile_convergence=paste(f$pl.conv,collapse=' | '),formula=paste(deparse(fm),collapse=' '))
      for(g in 0:1) counts[[paste(key,g)]]<-data.frame(cohort=co,outcome=y,adjustment_set=adj,reported=g,n=sum(z$report_group==g),events=sum(z[[y]][z$report_group==g]))
    }
  }
  r<-bind_rows(results)
  write_result(r,'firth_estimates',out)
  write_result(bind_rows(diagnostics),'firth_diagnostics',out)
  write_result(bind_rows(cells),'sparse_cell_diagnostics',out)
  write_result(bind_rows(counts),'outcome_counts',out)
  r
}

absolute_risks <- function(d,out,B=2000L) {
  # Set once, before the two outcomes, to reproduce the original bootstrap stream.
  RNGkind('Mersenne-Twister','Inversion','Rejection'); set.seed(42)
  risks<-diagnostics<-list()
  standardize<-function(m,z) {
    z0<-z1<-z;z0$report_group<-0;z1$report_group<-1
    p0<-mean(predict(m,newdata=z0,type='response'));p1<-mean(predict(m,newdata=z1,type='response'))
    c(unreported=p0,reported=p1,risk_difference=p1-p0,risk_ratio=p1/p0)
  }
  for(y in c('additional_fx_yes','clinic_or_med_90d')) {
    fm<-model_formula(y); z0<-model.frame(fm,d)
    fit<-glm(fm,data=z0,family=binomial());point<-standardize(fit,z0)
    strata<-split(seq_len(nrow(z0)),interaction(z0$center,z0$report_group,drop=TRUE))
    draws<-matrix(NA_real_,B,5,dimnames=list(NULL,c('OR',names(point))))
    flags<-data.frame(iteration=seq_len(B),converged=FALSE,full_rank=FALSE,finite=FALSE,extreme_probability=FALSE,warning='')
    for(b in seq_len(B)) {
      idx<-unlist(lapply(strata,function(i)sample(i,length(i),replace=TRUE)),use.names=FALSE)
      z<-z0[idx,,drop=FALSE];warns<-character()
      g<-tryCatch(withCallingHandlers(glm(fm,data=z,family=binomial(),control=glm.control(maxit=100)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
      flags$warning[b]<-paste(warns,collapse=' | ')
      if(inherits(g,'error')) {flags$warning[b]<-conditionMessage(g);next}
      flags$converged[b]<-g$converged;flags$full_rank[b]<-g$rank==length(coef(g))
      flags$finite[b]<-all(is.finite(coef(g)))
      flags$extreme_probability[b]<-any(fitted(g)<1e-8 | fitted(g)>1-1e-8)
      if(flags$converged[b] && flags$full_rank[b] && flags$finite[b]) draws[b,]<-c(exp(coef(g)['report_group']),standardize(g,z))
    }
    valid<-complete.cases(draws) & apply(draws,1,function(x)all(is.finite(x)))
    if(sum(valid)<=.95*B) stop('More than 5% bootstrap replicates invalid. Inspect diagnostics.')
    qs<-apply(draws[valid,,drop=FALSE],2,quantile,probs=c(.025,.975))
    risks[[y]]<-data.frame(outcome=y,estimand=names(point),estimate=unname(point),ci_low=qs[1,names(point)],ci_high=qs[2,names(point)])
    diagnostics[[y]]<-data.frame(outcome=y,attempted=B,valid=sum(valid),failed=sum(!valid),warnings=sum(nzchar(flags$warning)),extreme_probability=sum(flags$extreme_probability))
    write_result(cbind(iteration=seq_len(B),draws),paste0(y,'_bootstrap_draws'),out)
    write_result(flags,paste0(y,'_bootstrap_diagnostics'),out)
    cat('Bootstrap completed:',y,sum(valid),'/',B,'\n')
  }
  write_result(bind_rows(risks),'adjusted_risks',out)
  write_result(bind_rows(diagnostics),'bootstrap_summary',out)
}

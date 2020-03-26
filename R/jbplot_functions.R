jabba.libs <- function(){
  # Install required packages if missing
  list.of.packages <- c("gplots", "coda","rjags","R2jags","fitdistrplus","reshape","mvtnorm","scales",'snpar')
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  # Load Packages
  library(gplots);library(coda);library(rjags);library(R2jags);library("fitdistrplus");library(reshape);library('snpar')
  library(mvtnorm);library(scales)
  
}


#---------------------------------------------------------------------------------
# JABBA plotting functions
#---------------------------------------------------------------------------------- 


#----------------
# Total Landings
#----------------
jbplot_catch <- function(jabba,output.dir=getwd()){
  cat(paste0("\n","><> jbplot_catch()  <><","\n"))
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Landings_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  cord.x <- c(jabba$yr,rev(jabba$yr))
  y<-rep(0,length(jabba$yr))
  plot(jabba$yr,(jabba$catch),type="l",ylim=c(0,max(jabba$settings$TC,na.rm=T)),lty=1,lwd=1.3,xlab="Year",ylab=paste0("Catch ",jabba$settings$catch.metric),main="")
  polygon(cord.x,c(jabba$catch,rev(y)),col="gray",border=1,lty=1)
  dev.off()
}

#------------------------------
# Catch estimated with CV
#------------------------------
jbplot_catcherror <- function(jabba,output.dir=getwd()){  
  if(jabba$settings$add.catch.CV==TRUE){
    cat(paste0("\n","><> jbplot_catcherror()  <><","\n"))  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.7,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Catch.fit_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  # estimated Catch
  predC = jabba$est.catch 
  years = jabba$yr
  cord.x <- c(jabba$yr,rev(jabba$yr))
  cord.y<-c(predC[,3],rev(predC[,4]))
  plot(years,(jabba$catch),type="n",ylim=c(0,max(predC,na.rm=T)),lty=1,lwd=1.3,xlab="Year",ylab=paste0("Catch ",jabba$settings$catch.metric),main="")
  polygon(cord.x,cord.y,col="gray",border=0,lty=1)
  lines(years,predC[,2],lwd=2,col=4)
  points(years,(jabba$catch),pch=21,bg=0,cex=1.5)
  legend("topright",c("Observed","Predicted"),pch=c(21,-1),bg=0,lwd=c(-1,2),col=c(1,4),bty="n")
  dev.off() 
  } else {
    cat(paste0("\n","><> jbplot_catcherror() only available if add.catch.CV=TRUE <><","\n"))
}
}

#------------------------------
# Plot Posteriors
#------------------------------

jbplot_ppdist <- function(jabba, output.dir=getwd()){  
  cat(paste0("\n","><> jbplot_ppist() - prior and posterior distributions  <><","\n"))
  out =   jabba$pars_posterior
  node_id = names(out)
  #informative priors
  Prs = as.matrix(cbind(jabba$settings$K.pr,jabba$settings$r.pr,c(0,0),jabba$settings$psi.pr))
  
  #Posteriors
  Par = list(mfrow=c(round(length(node_id)/3+0.33,0),3),mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Posteriors_",jabba$assessment,"_",jabba$scenario,".png"),width  = 8, height = 2.5*round(length(node_id)/3,0), 
      res = 200, units = "in")
  par(Par)
  
  for(i in 1:length(node_id))
  {
    
    post.par = as.numeric(unlist(out[paste(node_id[i])]))
    
    if(i==1){
      
      rpr =  rlnorm(10000,log(Prs[1,i]),Prs[2,i]) 
      pdf = stats::density(post.par,adjust=2)  
      prior = dlnorm(sort(rpr),log(Prs[1,i]),Prs[2,i])   
      plot(pdf,type="l",ylim=range(prior,pdf$y),xlim=range(c(pdf$x,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab="K",ylab="",xaxs="i",yaxs="i",main="")
      
      polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
      polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
      legend('right',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
      PPVM = round(mean(post.par)/mean(rpr),3)
      legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")  
      
    }  
    if(i==2){
      rpr = rlnorm(10000,log(Prs[1,i]),Prs[2,i]) 
      pdf = stats::density(post.par,adjust=2) 
      prior = dlnorm(sort(rpr),log(Prs[1,i]),Prs[2,i])   
      plot(pdf$x,pdf$y,type="l",ylim=range(prior,pdf$y),xlim=range(c(post.par,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")
      polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
      polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
      PPVM = round(mean(post.par)/mean(rpr),3)
      legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
      
      
    }
    
    if(i==3){
      if(jabba$settings$model.id<4){
        plot(1,1,type="n",xlim=range(0.5,2.5),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")  
        abline(v=jabba$pars["m",1],lwd=2)}
      if(jabba$settings$model.id==4){
        mpr = rlnorm(10000,log(jabba$jagsdata$mu.m),jabba$jagsdata$m.CV) 
        pdf = stats::density(post.par,adjust=2) 
        prior = dlnorm(sort(mpr),log(m),shape.CV)   
        plot(pdf$x,pdf$y,type="l",ylim=range(prior,pdf$y),xlim=range(c(post.par,quantile(rpr,c(0.0001,0.95)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")
        
        polygon(c(sort(mpr),rev(sort(mpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
        polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
        PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
        PPVM = round(mean(post.par)/mean(rpr),3)
        legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
        
      }
    }
    
    
    if(i==4){
      if(jabba$settings$psi.dist=="beta"){
        parm = fitdist(post.par[post.par<1 & post.par>0.01], "beta")$estimate
        rpr = rbeta(10000,(Prs[1,4]),Prs[2,4]) 
        pdf = stats::density(post.par,adjust=2)  
        prior = dbeta(sort(rpr),psi.pr[1],psi.pr[2])   
        PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
        PPVM = round(mean(post.par)/mean(rpr),3)
        legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
        
      } else {
        rpr = rlnorm(10000,log(Prs[1,4]),Prs[2,4]) 
        pdf = stats::density(post.par,adjust=2)  
        prior = dlnorm(sort(rpr),log(Prs[1,4]),Prs[2,4])}
      plot(pdf,type="l",ylim=range(quantile(c(prior,pdf$y,c(0,0.95)))),xlim=range(c(0.5,post.par,pdf$x,quantile(rpr,c(0.001,0.999)))),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
      polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.4,1))
      polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
      PPVM = round(mean(post.par)/mean(rpr),3)
      legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
      
      #legend('topright',c("Prior","Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
    }        
    
    if(i>4){
      if(jabba$settings$sigma.proc!=TRUE & i==length(node_id)) {
        plot(1,1,type="n",xlim=range(0,0.15^2),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i")  
        abline(v=sigma.proc^2,lwd=2)} else {
          
          pdf = stats::density(post.par,adjust=2)  
          plot(pdf,type="l",xlim=range(0,post.par),yaxt="n",xlab=paste(node_id[i]),ylab="",xaxs="i",yaxs="i",main="")
          if(i==length(node_id)& jabba$settings$igamma[1]>0.9){
            rpr = 1/rgamma(10000,jabba$settings$igamma[1],jabba$settings$igamma[2])
            prior = stats::density(rpr,adjust=2)
            polygon(c(prior$x,rev(prior$x)),c(prior$y,rep(0,length(prior$y))),col=gray(0.4,1))
            PPVR = round((sd(post.par)/mean(post.par))^2/(sd(rpr)/mean(rpr))^2,3)  
            PPVM = round(mean(post.par)/mean(rpr),3)
            legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=1,bty="n")
            
          }
          
          polygon(c(pdf$x,rev(pdf$x)),c(pdf$y,rep(0,length(pdf$y))),col=gray(0.7,0.7))
          #legend('topright',c("Posterior"),pch=22,pt.cex=1.5,pt.bg = c(grey(0.8,0.6)),bty="n")
        } }         
    
  }
  mtext(paste("Density"), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
  dev.off()   
} # End of ppdist plot


#-----------------------------
# MCMC chains of posteriors
#-----------------------------

jbplot_mcmc <- function(jabba, output.dir=getwd()){
  cat(paste0("\n","><> jbplot_mcmc() - mcmc chains  <><","\n"))
  out =   jabba$pars_posterior
  node_id = names(out)
  
  Par = list(mfrow=c(round(length(node_id)/3+0.33,0),3),mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/MCMC_",jabba$assessment,"_",jabba$scenario,".png"), width = 8, height = 2.5*round(length(node_id)/3,0), 
      res = 200, units = "in")
  par(Par)
  for(i in 1:length(node_id)){
    
    post.par = as.numeric(unlist(out[paste(node_id[i])]))
    plot(out[,i],xlab=paste(node_id[i]),ylab="",type="l",col=4)
    lines(rep(mean(out[,i]),length(out[,i])),col=2,lwd=2)   
  }
  dev.off()
}

#---------------------------------------------------------------------------
# Plot CPUE with expectected mean CIs and Posterior Predictive Distributions
#----------------------------------------------------------------------------

jbplot_cpue <- function(jabba, output.dir=getwd()){  
  if(jabba$settings$CatchOnly==FALSE){
    cat(paste0("\n","><> jbplot_cpue() - fits to CPUE <><","\n"))
    
    N = jabba$settings$N 
    years = jabba$yr
    n.indices = jabba$settings$nI
    series = 1:jabba$settings$nI
    CPUE = jabba$settings$I 
    indices = unique(jabba$diags$name)
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    #CPUE FITS
    Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
    png(file = paste0(output.dir,"/Fits_",jabba$assessment,"_",jabba$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0), 
        res = 200, units = "in")
    par(Par)
    
    
    for(i in 1:n.indices){
      
      # set observed vs predicted CPUE
      Yr = jabba$yr
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1
      
      fit =  t(jabba$cpue.ppd[,c(2,1,3),i])   
      fit.hat = t(jabba$cpue.hat[,c(2,1,3),i])
      mufit = mean(fit[2,])
      fit = fit/mufit
      fit.hat = fit.hat/mufit
      
      cpue.i = CPUE[is.na(CPUE[,i])==F,i]
      yr.i = Yr[is.na(CPUE[,i])==F]
      se.i = sqrt(jabba$settings$SE2[is.na(CPUE[,i])==F,(i)])
      
      ylim = c(min(fit*0.9,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit*1.05,exp(log(cpue.i)+1.96*se.i)/mufit))
      
      cord.x <- c(Yr,rev(Yr))
      cord.y <- c(fit[1,yr],rev(fit[3,yr]))
      cord.yhat <- c(fit.hat[1,yr],rev(fit.hat[3,yr]))
      # Plot Observed vs predicted CPUE
      plot(years,CPUE[,i],ylab="",xlab="",ylim=ylim,xlim=range(jabba$yr),type='n',xaxt="n",yaxt="n")
      axis(1,labels=TRUE,cex=0.8)
      axis(2,labels=TRUE,cex=0.8)
      polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
      polygon(cord.x,cord.yhat,col=grey(0.3,0.5),border=grey(0.3,0.5),lty=2)
      
      lines(Yr,fit[2,yr],lwd=2,col=1)
      if(jabba$settings$SE.I  ==TRUE | max(jabba$settings$SE2)>0.01){ plotCI(yr.i,cpue.i/mufit,ui=exp(log(cpue.i)+1.96*se.i)/mufit,li=exp(log(cpue.i)-1.96*se.i)/mufit,add=T,gap=0,pch=21,xaxt="n",yaxt="n")}else{
        points(yr.i,cpue.i/mufit,pch=21,xaxt="n",yaxt="n",bg="white")}
      
      legend('top',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.9)
      }
    
    mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
    mtext(paste("Normalized Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
    dev.off()
  } else {
    cat(paste0("\n","><> jbplot_cpue() not available CatchOnly=TRUE <><","\n"))
  }
} # End of CPUE plotting function

#-------------------------
# Plot logfits CPUE
#-------------------------

jbplot_logfits <- function(jabba, output.dir= getwd()){  
  if(jabba$settings$CatchOnly==FALSE){
    cat(paste0("\n","><> jbplot_logfits()  <><","\n"))
    
    N = jabba$settings$N 
    years= jabba$yr
    n.indices = jabba$settings$nI
    series = 1:jabba$settings$nI
    CPUE = jabba$settings$I 
    indices = unique(jabba$diags$name)
    check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
    cpue.yrs = years[check.yrs>0]
    
    #log CPUE FITS
    Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
    png(file = paste0(output.dir,"/logFits_",jabba$assessment,"_",jabba$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0), 
        res = 200, units = "in")
    par(Par)
    
    
    for(i in 1:n.indices){
      
      Yr = jabba$yr
      Yr = min(Yr):max(Yr)
      yr = Yr-min(years)+1
      
      fit = t(jabba$cpue.hat[,c(2,1,3),i])
      mufit = mean(fit[2,])
      fit = fit/mufit
      cpue.i = CPUE[is.na(CPUE[,i])==F,i]
      yr.i = Yr[is.na(CPUE[,i])==F]
      se.i = sqrt(jabba$settings$SE2[is.na(CPUE[,i])==F,(i)])
      
      ylim = log(c(min(fit[,yr[is.na(CPUE[,i])==F]]*0.8,exp(log(cpue.i)-1.96*se.i)/mufit), max(fit[,yr[is.na(CPUE[,i])==F]]*1.3,exp(log(cpue.i)+1.96*se.i)/mufit)))
      
      # Plot Observed vs predicted CPUE
      plot(years,CPUE[,i],ylab="",xlab="",ylim=ylim,xlim=range(yr.i),type='n',xaxt="n",yaxt="n")
      axis(1,labels=TRUE,cex=0.8)
      axis(2,labels=TRUE,cex=0.8)
      #polygon(cord.x,cord.y,col=grey(0.5,0.5),border=0,lty=2)
      
      
      lines(Yr,log(fit[2,yr]),lwd=2,col=4)
      if(jabba$settings$SE.I ==TRUE | max(jabba$settings$SE2)>0.01){ plotCI(yr.i,log(cpue.i/mufit),ui=log(exp(log(cpue.i)+1.96*se.i)/mufit),li=log(exp(log(cpue.i)-1.96*se.i)/mufit),add=T,gap=0,pch=21,xaxt="n",yaxt="n")}else{
        points(yr.i,log(cpue.i/mufit),pch=21,xaxt="n",yaxt="n",bg="white")}
      legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
    }
    
    mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
    mtext(paste("Log Index"), side=2, outer=TRUE, at=0.5,line=1,cex=1)
    dev.off()
  }  else {
    cat(paste0("\n","><> jbplot_logfit() not available CatchOnly=TRUE <><","\n"))
  }
} # End of logfit


#-------------------------------------------------
# JABBA Residual Plot
#-------------------------------------------------
jbplot_residuals <- function(jabba, output.dir=getwd()){
  if(jabba$settings$CatchOnly==FALSE){
  cat(paste0("\n","><> jbplot_residuals() - JABBA residual plot  <><","\n"))
  years = jabba$yr  
  check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
  cpue.yrs = years[check.yrs>0]
  Resids = jabba$residuals
  Yr = jabba$yr
  n.years = length(Yr) 
  n.indices = jabba$settings$nI
  indices = unique(jabba$diags$name)
  series = 1:jabba$settings$nI
  
  # JABBA-residual plot
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Residuals_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  
  
  plot(Yr,Yr,type = "n",ylim=ifelse(rep(max(Resids,na.rm = T),2)>0.9,range(1.2*Resids,na.rm = T),range(c(-1.3,1.2))),xlim=range(cpue.yrs),ylab="log residuals",xlab="Year")
  boxplot(Resids,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
  abline(h=0,lty=2)
  positions=runif(n.indices,-0.2,0.2)
  
  for(i in 1:n.indices){
    for(t in 1:n.years){
      lines(rep((Yr+positions[i])[t],2),c(0,Resids[i,t]),col=jabba$settings$cols[i])}
    points(Yr+positions[i],Resids[i,],col=1,pch=21,bg=jabba$settings$cols[i])}
  mean.res = apply(Resids[,as.numeric(colnames(Resids))%in%cpue.yrs],2,mean,na.rm =TRUE)
  smooth.res = predict(loess(mean.res~cpue.yrs),data.frame(cpue.yrs=cpue.yrs))
  lines(cpue.yrs,smooth.res,lwd=2)
  # get degree of freedom
  Nobs =length(as.numeric(Resids)[is.na(as.numeric(Resids))==FALSE])
  
  RMSE = round(jabba$stats[5,2],1)
  
  legend('topright',c(paste0("RMSE = ",RMSE,"%")),bty="n")
  legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,pt.cex=1.1,cex=0.75,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba$settings$col[series],1),lwd=c(rep(-1,n.indices),2))
  dev.off()
  }  else {
    cat(paste0("\n","><> jbplot_residuals() not available CatchOnly=TRUE <><","\n"))
  }
} # End of functions

#---------------------------------------
# Plot Stadardized Residuals
#--------------------------------------

jbplot_stdresiduals <- function(jabba, output.dir=getwd()){
  if(jabba$settings$CatchOnly==FALSE){
    cat(paste0("\n","><> jbplot_staresiduals() - standardized residuals  <><","\n"))
  years = jabba$yr  
  check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
  cpue.yrs = years[check.yrs>0]
  Resids = jabba$residuals
  Yr = jabba$yr
  n.years = length(Yr) 
  n.indices = jabba$settings$nI
  indices = unique(jabba$diags$name)
  series = 1:jabba$settings$nI
  StResid = jabba$std.residuals
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/StandardizedResids_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  # Standardized Residuals
  plot(Yr,Yr,type = "n",ylim=c(min(-1,-1.2*max(abs(StResid),na.rm = T)),max(1,1.2*max(abs(StResid),na.rm = T))),xlim=range(cpue.yrs),ylab="Standardized residuals",xlab="Year")
  boxplot(StResid,add=TRUE,at=c(Yr),xaxt="n",col=grey(0.8,0.5),notch=FALSE,outline = FALSE)
  abline(h=0,lty=2)
  positions=runif(n.indices,-0.2,0.2)
  
  for(i in 1:n.indices){
    for(t in 1:n.years){
      lines(rep((Yr+positions[i])[t],2),c(0,StResid[i,t]),col=jabba$settings$cols[i])}
    points(Yr+positions[i],StResid[i,],col=1,pch=21,bg=jabba$settings$cols[i])}
  mean.res = apply(StResid[,as.numeric(colnames(StResid))%in%cpue.yrs],2,mean,na.rm =TRUE)
  smooth.res = predict(loess(mean.res~cpue.yrs),data.frame(cpue.yrs=cpue.yrs))
  lines(cpue.yrs,smooth.res,lwd=2)
  SDNR = round(sqrt(sum(StResid^2,na.rm =TRUE)/(jabba$stats[1,2]-1)),2)
  Crit.value = (qchisq(.95, df=(jabba$stats[1,2]-1))/(jabba$stats[1,2]-1))^0.5
  legend('topright',c(paste0("SDNR = ",SDNR,"(",round(Crit.value,2),")")),bty="n")
  legend('bottomright',c(paste(indices),"Loess"),bty="n",col=1,cex=0.75,pt.cex=1.1,pch=c(rep(21,n.indices),-1),pt.bg=c(jabba$settings$cols[series],1),lwd=c(rep(-1,n.indices),2))
  dev.off()
  }  else {
    cat(paste0("\n","><> jbplot_stdresiduals() not available CatchOnly=TRUE <><","\n"))
  }
  
} # end of function



#-------------------------------------------------
# Function to do runs.test and 3 x sigma limits  
#------------------------------------------------
runs.sig3 <- function(x,type=NULL) {
  if(is.null(type)) type="resid"
  if(type=="resid"){mu = 0}else{mu = mean(x, na.rm = TRUE)} 
  # Average moving range
  mr  <- abs(diff(x - mu))
  amr <- mean(mr, na.rm = TRUE)
  # Upper limit for moving ranges
  ulmr <- 3.267 * amr
  # Remove moving ranges greater than ulmr and recalculate amr, Nelson 1982
  mr  <- mr[mr < ulmr]
  amr <- mean(mr, na.rm = TRUE)
  # Calculate standard deviation, Montgomery, 6.33
  stdev <- amr / 1.128
  # Calculate control limits
  lcl <- mu - 3 * stdev
  ucl <- mu + 3 * stdev
  if(nlevels(factor(sign(x)))>1){ 
    runstest = snpar::runs.test(x) 
    pvalue = round(runstest$p.value,3)} else {
      pvalue = 0.001  
    }
  
  return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
}

#----------------------------------------------------
# Runs test plots
#---------------------------------------------------


jbplot_runstest <- function(jabba, output.dir=getwd()){
  if(jabba$settings$CatchOnly==FALSE){
    cat(paste0("\n","><> jbplot_runstest()   <><","\n"))

    
  years = jabba$yr  
 check.yrs = abs(apply(jabba$residuals,2,sum,na.rm=TRUE))
 cpue.yrs = years[check.yrs>0]
 Resids = jabba$residuals
 n.years = length(years) 
 n.indices = jabba$settings$nI
 indices = unique(jabba$diags$name)
 series = 1:jabba$settings$nI


Par = list(mfrow=c(round(n.indices/2+0.01,0),ifelse(n.indices==1,1,2)),mai=c(0.35,0.15,0,.15),omi = c(0.2,0.25,0.2,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Residual_RunsTests_",jabba$assessment,"_",jabba$scenario,".png"), width = 7, height = ifelse(n.indices==1,5,ifelse(n.indices==2,3.,2.5))*round(n.indices/2+0.01,0), 
    res = 200, units = "in")
par(Par)


for(i in 1:n.indices){
  
  resid = (Resids[i,is.na(Resids[i,])==F])  
  res.yr = years[is.na(Resids[i,])==F]
  runstest = runs.sig3(x=as.numeric(resid),type="resid")
  # CPUE Residuals with runs test
  plot(res.yr,rep(0,length(res.yr)),type="n",ylim=c(min(-1,runstest$sig3lim[1]*1.25),max(1,runstest$sig3lim[2]*1.25)),lty=1,lwd=1.3,xlab="Year",ylab=expression(log(cpue[obs])-log(cpue[pred])))
  abline(h=0,lty=2)
  lims = runstest$sig3lim
  cols =  c(rgb(1,0,0,0.5),rgb(0,1,0,0.5))[ifelse(runstest$p.runs<0.05,1,2)]
  rect(min(years-1),lims[1],max(years+1),lims[2],col=cols,border=cols) # only show runs if RMSE >= 0.1
  for(j in 1:length(resid)){
    lines(c(res.yr[j],res.yr[j]),c(0,resid[j]))  
  }
  points(res.yr,resid,pch=21,bg=ifelse(resid < lims[1] | resid > lims[2],2,"white"),cex=1)
  legend('topright',paste(indices[i]),bty="n",y.intersp = -0.2,cex=0.8)
}  


mtext(paste("Year"), side=1, outer=TRUE, at=0.5,line=1,cex=1)
mtext(expression(log(cpue[obs])-log(cpue[pred])), side=2, outer=TRUE, at=0.5,line=1,cex=1)
dev.off()
  }  else {
    cat(paste0("\n","><> jbplot_stdresiduals() not available CatchOnly=TRUE <><","\n"))
  }
  
} # end of function




#------------------------------
# Plot process error deviation
#------------------------------

jbplot_procdev <- function(jabba, output.dir=getwd()){ 
  cat(paste0("\n","><> jbplot_procdev() - Process error diviations on log(biomass)  <><","\n"))
  
  years=jabba$yr
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/ProcDev_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  ylim = c(min(-0.22,jabba$timeseries[,,"procB"]),max(0.22,jabba$timeseries[,,"procB"]))#range(proc.dev)*1.1
  cord.x <- c(years,rev(years))
  cord.y <- c(jabba$timeseries[,2,"procB"],rev(jabba$timeseries[,3,"procB"]))
  # Process Error
  plot(years,jabba$timeseries[,1,"procB"],ylab="Process Error Deviates",xlab="Year",ylim=ylim,type="n")
  polygon(cord.x,cord.y,col='grey',border=0,lty=2)
  lines(years,jabba$timeseries[,1,"procB"],lwd=2)
  lines(years,rep(0,length(years)),lty=5)
  dev.off()
} # end of plot function


#-----------------------------------------------------------
# <><<><<><<><<>< JABBA Management Plots ><>><>><>><>><>><>
#-----------------------------------------------------------

#-----------------------
# Plot trajectores: B,F,Bmsy,FFmsy,BB0 
#-----------------------

jbplot_trj <-  function(jabba, type = c("B","F","BBmsy","FFmsy","BB0") ,output.dir=getwd()){ 
  
  for(i in 1:length(type)){  
    
    Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.5, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
    png(file = paste0(output.dir,"/",type[i],"_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 3.5, 
        res = 200, units = "in")
    par(Par)
    cat(paste0("\n","><> jbplot_trj() - ", type[i]," trajectory  <><","\n"))
    
    ylabs = c(paste("Biomass",jabba$settings$catch.metric),ifelse(jabba$settings$harvest.label=="Fmsy","Fishing mortality F","Harvest rate H"),expression(B/B[MSY]),ifelse(jabba$settings$harvest.label=="Fmsy",expression(F/F[MSY]),expression(H/h[MSY])),expression(B/B[0])) 
    trj = jabba$timeseries[,,paste(type[i])] 
    years = jabba$yr
    ylim = c(0, max(trj[,3]))
    cord.x <- c(years,rev(years))
    cord.y <- c(trj[,2],rev(trj[,3]))
    plot(years,trj[,1],ylab=ylabs[i],xlab="Year",ylim=ylim,type="n")
    polygon(cord.x,cord.y,col='grey',border=0,lty=2)
    lines(years,trj[,1],lwd=2,col=1)
    if(i==1) lines(years,rep(jabba$refpts$bmsy[1],length(years)),lty=5)
    if(i==2) lines(years,rep(jabba$refpts$fmsy[1],length(years)),lty=5)
    if(i>2 & i <5) lines(years,rep(1,length(years)),lty=5)
    if(i==1) text((max(years)-min(years))/30+years[1],jabba$refpts$bmsy[1]*1.11,expression(paste(B[MSY])))
    if(i==2) text((max(years)-min(years))/30+years[1],jabba$refpts$fmsy[1]*1.11,ifelse(jabba$settings$harvest.label=="Fmsy",expression(F[MSY]),expression(H[MSY])))
    if(i==5){
      lines(years,rep(jabba$refpts$bmsy[1]/jabba$refpts$k[1],length(years)),lty=5)
      text((max(years)-min(years))/30+years[1],jabba$refpts$bmsy[1]/jabba$refpts$k[1]*1.11,expression(paste(B[MSY])))
    }
    dev.off()
  }} # end of plot function



#-----------------------------------------
# Produce JABBA SP-phase plot
#-----------------------------------------
jbplot_spphase <-  function(jabba ,output.dir=getwd()){ 
  cat(paste0("\n","><> jbplot_spphase() - JABBA Surplus Production Phase Plot  <><","\n"))
  
  # extract pars  
  m = jabba$pfun$m[1]
  Bit = jabba$pfun$SB_i
  Cmsy = jabba$pfunc$Cmsy
  B = jabba$timeseries[,"mu","B"]
  Hmsy.sp = jabba$pfunc$Hmsy[1] 
  SB0.sp =jabba$pfunc$SB0[1]
  SP = jabba$pfunc$SP
  Bmsy.sp = jabba$estimates["SBmsy",1]
  MSY.sp = jabba$estimates["MSY",1]
  N = jabba$settings$N
  years = jabba$yr
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/SPphase_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  
  green.x = c(max(Bit,B),max(Bit,B),Bmsy.sp,Bmsy.sp,max(Bit))
  green.y = c(Bmsy.sp,0,0,max(SP),max(Cmsy))
  red.x = c(0,0,Bmsy.sp,Bmsy.sp,0)
  red.y = c(SB0.sp,0,max(SP),SB0.sp,SB0.sp)
  plot(Bit,SP,type = "n",ylim=c(0,max(c(max(jabba$catch,na.rm=T)*1.05,max(MSY.sp*1.1)))),xlim=c(0,max(Bit,B)),ylab="Surplus Production",xlab="Spawning Biomass",xaxs="i",yaxs="i")
  rect(0,0,SB0.sp*1.1,SB0.sp*1.1,col="green",border=0)
  rect(0,0,SB0.sp,SB0.sp,col="yellow",border=0)
  rect(0,max(SP),SB0.sp,SB0.sp,col="orange",border=0)
  polygon(green.x,green.y,border = 0,col="green")
  polygon(red.x,red.y,border = 0,col="red")
  
  ry.sp = Bit[Bit<=Bmsy.sp]
  for(i in 1:length(ry.sp)){
    
    lines(rep(Bit[i],2),c(Cmsy[i],SP[i]),col=ifelse(i %% 2== 0,"yellow","red"),lty=3)  
    #i = i+1
  }
  
  polygon(c(-10000,10^7,10^7,-10000),c(rep(MSY.sp[1],2),rep(MSY.sp[3],2)),border = FALSE,col=rgb(0,0,1,0.4))
  lines(Bit,SP,col=4,lwd=2)
  lines(B,jabba$catch,lty=1,lwd=1)
  points(B,jabba$catch,cex=0.8,pch=16)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(B[sel.yr],jabba$catch[sel.yr],col= 1,pch=c(22,21,24),bg="white",cex=1.7)
  abline(h=max(SP),col=4,lty=5)
  sel.years =years[sel.yr]
  lines(rep(Bmsy.sp,2),c(-1000,max(SP)),lty=2,col=4)
  
  legend('topright', 
         c(expression(B[MSY]),"MSY","SP","Catch",paste(sel.years)), 
         lty=c(2,5,1,1,1,1,1),pch=c(-1,-1,-1,16,22,21,24),pt.bg=c(0,0,0,0,rep("white",3)), 
         col=c(4,4,4,rep(1,4)),lwd=c(1,1,2,1,1,1),cex=0.8,pt.cex=c(-1,-1,-1,0.5,rep(1.3,3)),bty="n")
  
  dev.off()
} #end of plotting function


#------------------------------------------------------
# status plot (kobe, biplot)
#------------------------------------------------------

jbplot_kobe <-  function(jabba ,output.dir=getwd()){ 
  cat(paste0("\n","><> jbplot_kobe() - Stock Status Plot  <><","\n"))
  
  mu.f = jabba$timeseries[,,"FFmsy"]   
  mu.b = jabba$timeseries[,,"BBmsy"]
  f = jabba$kobe$harvest
  b = jabba$kobe$stock
  years=jabba$yr
  N = length(years)
  # fit kernel function
  kernelF <- gplots::ci2d(b,f,nbins=151,factor=1.5,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,xlab= ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])))
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(output.dir,"/Kobe_",jabba$assessment,"_",jabba$scenario,".png"), width = 5, height = 4.5, 
      res = 200, units = "in")
  par(Par)
  
  #Create plot
  plot(1000,1000,type="b", xlim=c(0,max(1/(jabba$refpts$bmsy/jabba$refpts$k)[1],mu.b[,1]) +0.05), ylim=c(0,max(mu.f[,1],quantile(f,0.85),2.)),lty=3,ylab=ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab=expression(paste(B/B[MSY])),xaxs="i",yaxs="i")
  c1 <- c(-1,100)
  c2 <- c(1,1)
  
  # extract interval information from ci2d object
  # and fill areas using the polygon function
  zb2 = c(0,1)
  zf2  = c(1,100)
  zb1 = c(1,100)
  zf1  = c(0,1)
  polygon(c(zb1,rev(zb1)),c(0,0,1,1),col="green",border=0)
  polygon(c(zb2,rev(zb2)),c(0,0,1,1),col="yellow",border=0)
  polygon(c(1,100,100,1),c(1,1,100,100),col="orange",border=0)
  polygon(c(0,1,1,0),c(1,1,100,100),col="red",border=0)
  
  polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
  polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
  polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
  
  points(mu.b[,1],mu.f[,1],pch=16,cex=1)
  lines(c1,c2,lty=3,lwd=0.7)
  lines(c2,c1,lty=3,lwd=0.7)
  lines(mu.b[,1],mu.f[,1], lty=1,lwd=1.)
  sel.yr = c(1,round(quantile(1:N,0.7),0),N)
  points(mu.b[sel.yr,1],mu.f[sel.yr,1],col=
           1,pch=c(22,21,24),bg="white",cex=1.9)
  
  # Get Propability
  Pr.green = sum(ifelse(b>1 & f<1,1,0))/length(b)*100
  Pr.red = sum(ifelse(b<1 & f>1,1,0))/length(b)*100
  Pr.yellow = sum(ifelse(b<1 & f<1,1,0))/length(b)*100
  Pr.orange = sum(ifelse(b>1 & f>1,1,0))/length(b)*100
    
  
  
  sel.years = c(years[sel.yr])
  ## Add legend
  
    legend('topright', 
           c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.orange,Pr.green),1),"%")), 
           lty=c(1,1,1,rep(-1,8)),pch=c(22,21,24,rep(22,8)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","orange","green"), 
           col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.2,4)),bty="n")  
    
  dev.off()
} # End of Kobe plot 


#---------------------------------------------------------
# Produce 'post-modern' biplot (see Quinn and Collie 2005)
#---------------------------------------------------------

jbplot_biplot <-  function(jabba ,output.dir=getwd()){ 
cat(paste0("\n","><> jbplot_biplot() - Stock Status Plot  <><","\n"))
mu.f = jabba$timeseries[,,"FFmsy"]   
mu.b = jabba$timeseries[,,"BBmsy"]
f = jabba$kobe$harvest
b = jabba$kobe$stock
years=jabba$yr
N = length(years)

# fit kernel function
kernelF <- gplots::ci2d(f,b,nbins=201,factor=2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,ylab= ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),xlab=expression(paste(B/B[MSY])))


Par = list(mfrow=c(1,1),mai=c(0.2,0.15,0,.15),omi = c(0.3,0.25,0.2,0) + 0.1, mgp =c(3,1,0), tck = -0.02,cex=0.8)
png(file = paste0(output.dir,"/Biplot_",jabba$assessment,"_",jabba$cenario,".png"), width = 5, height = 4.5, 
    res = 200, units = "in")
par(Par)

#Create plot
plot(1000,1000,type="b", ylim=c(0,max(1/(jabba$refpts$bmsy/jabba$refpts$k)[1],mu.b[,1]) +0.05), xlim=c(0,max(mu.f[,1],quantile(f,0.85),2.)),lty=3,xlab=ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])),xaxs="i",yaxs="i")

# and fill areas using the polygon function
fint = seq(0.001,100,0.01)
# read ftarget,bthreshold
ftarget<-0.8
bthreshold<-0.2

#Zone X
xb=bthreshold+(1.0-bthreshold)/ftarget*fint
xf =  ifelse(xb>1,0.8,fint)
polygon(c(0,0,xf),c(max(xb),bthreshold,xb),col="green")
zb = bthreshold+(1.0-bthreshold)*fint
zf  = ifelse(zb>1,1,fint) 
polygon(c(zf,rep(max(fint),2),rep(0,2)),c(zb,max(zb),0,0,bthreshold),col="red")

polygon(c(xf,rev(zf)),c(xb,rev(zb)),col="yellow")

c1 <- c(-1,100)
c2 <- c(1,1)

# extract interval information from ci2d object
# and fill areas using the polygon function
polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
points(mu.f[2,],mu.b[2,],pch=16,cex=1)

lines(c1,c2,lty=3,lwd=0.7)
lines(c2,c1,lty=3,lwd=0.7)
lines(mu.f[,1],mu.b[,1], lty=1,lwd=1.)
sel.yr = c(1,round(quantile(1:N,0.7),0),N)
points(mu.f[sel.yr,1],mu.b[sel.yr,1],col=
         1,pch=c(22,21,24),bg="white",cex=1.9)


sel.years = years[sel.yr]
## Add legend
legend('topright', 
       c(paste(sel.years),"50% C.I.","80% C.I.","95% C.I."), 
       lty=c(1,1,1,-1,-1,-1),pch=c(22,21,24,22,22,22),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4"), 
       col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,4),1.7,1.7,1.7),bty="n")

Zone  = NULL
Status = NULL
X  = 0.15
Y = 0
Z = -0.15

for(i  in 1:length(f))
{
  if(b[i]>1.0){
    if(f[i]<ftarget){
      Zone[i]<-X
    } else if (f[i]>1.0){
      Zone[i]<-Z
    } else {
      Zone[i]<-Y
    }
  } else {
    if(b[i]>bthreshold+(1.0-bthreshold)/ftarget*f[i]){
      Zone[i]<-X
    } else if(b[i]<bthreshold+(1.0-bthreshold)*f[i]){
      Zone[i]<-Z
    } else {
      Zone[i]<-Y
    }
    
    
  }}

perGreen = round(length(Zone[Zone==0.15])/length(Zone)*100,1) 
perYellow = round(length(Zone[Zone==0])/length(Zone)*100,1) 
perRed = round(length(Zone[Zone==-0.15])/length(Zone)*100,1)

mtext(expression(paste(B/B[MSY])), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
mtext(ifelse(jabba$settings$harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))), side=1, outer=TRUE, at=0.5,line=1,cex=0.9)

text(0.65,2.4,paste0(perGreen,"%"))
text(0.9,2.4,paste0(perYellow,"%"))
text(1.2,2.4,paste0(perRed,"%"))

dev.off()

} # End of biplot function

#-------------------------
# wrapper plot function
#------------------------
jbplots = function(jabba,output.dir = getwd(),statusplot ="kobe"){
jbplot_catch(jabba,output.dir=output.dir) # catch.metric
jbplot_catcherror(jabba,output.dir=output.dir) # posteriors
jbplot_cpue(jabba,output.dir=output.dir) # check years
jbplot_logfits(jabba,output.dir=output.dir) # check n.indices
jbplot_mcmc(jabba,output.dir=output.dir)
jbplot_ppdist(jabba,output.dir=output.dir) # check m
jbplot_procdev(jabba,output.dir=output.dir)
jbplot_trj(jabba,output.dir=output.dir) 
jbplot_spphase(jabba,output.dir=output.dir) # check TC
jbplot_residuals(jabba,output.dir=output.dir) # check years
jbplot_stdresiduals(jabba,output.dir=output.dir)
jbplot_runstest(jabba,output.dir=output.dir)
if(statusplot =="kobe"){
jbplot_kobe(jabba,output.dir=output.dir)} else {
jbplot_biplot(jabba,output.dir=output.dir)}
}


#-----------------------
# Plot B/Bmsy and H/Hmsy
#-----------------------

#mu.f = jabba$timeseries[,,"FFmsy"]   
#mu.b = jabba$timeseries[,,"BBmsy"]
#f = jabba$kobe$harvest
#b = jabba$kobe$stock
#years=jabba$yr

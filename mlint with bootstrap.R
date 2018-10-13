########  ml.int em icin bootstrap yapiyrm


library(MASS)
library(survival)
library(joineR)
library(nlme)
library(parallel)
library(statmod)
library(dplyr)


simjointml<- function(m, noh, beta, beta2, sig.v, sig.hv, sig.err, gamma, ntms=8, theta0=-3,censlam=exp(-3), censoring=T, truncation=T){
  trunctime=max(ntms)
  n<- m*noh
  vi<-rnorm(n,0,sig.v)#vi<-mvrnorm(n ,mu=c(0,0),Sigma=var.v)
  vih<-rnorm(noh, 0, sig.hv)
  ctsx<-rnorm(n)
  grp<-rbinom(n,1,0.5)# binary covariate for treatment
  err<-rnorm(n*ntms, 0, sig.err)
  t<-c(0:(ntms-1)) #time covariate
  y<- rep(0, n*ntms)
  X<- cbind(1,rep(t,n),rep(grp, each=ntms),rep(ctsx,each=ntms))#create design matrix
  Xbeta<- X%*%beta#fixed effects
  b<- cbind(rep(vi,each=ntms))#b<- cbind(rep(vi[,1],each=10), rep(vi[,2],each=10))
  h1<-cbind(rep(vih, each=ntms*m))
  Y<- Xbeta+(b+h1)+err #generate Y
  ##### generate data for coxph
  X2 <- cbind(ctsx,grp)
  Xbeta2 <-X2%*%beta2
  svih<-rep(vih,each=m)
  uu<-runif(n)
  cens<-rep(1,n)
  survtime<--log(uu)/exp(theta0+Xbeta2+gamma[1]*vi+gamma[2]*svih)
  censtime <- -log(runif(n))/censlam
  censtime <- pmin(censtime,trunctime)
  ii<-censtime<survtime
  survtime[ii]<-censtime[ii]
  cens[ii]<-0
  
  
  idl<-rep(1:n,each=ntms)
  hosidl<-rep(1:noh, each=ntms*m)
  obstime<-rep(0:(ntms-1),n)
  lsurvtime<-rep(survtime, each=ntms)
  l.cens<-rep(cens, each=ntms)
  ctsxl<-rep(ctsx,each=ntms)
  grpl<-rep(grp,each=ntms)
  ul<-rep(vi, each=ntms)
  uhl<-rep(vih, each=m*ntms)
  ###put all the data into a dataframe
  
  obstimecopy<-rep(0:(ntms-1),n)
  Y<- Y[obstime<=lsurvtime]
  idl<- idl[obstime<=lsurvtime]
  hosidl<-hosidl[obstime<=lsurvtime]
  l.cens<-l.cens[obstime<=lsurvtime]
  ctsxl<-ctsxl[obstime<=lsurvtime]
  grpl<- grpl[obstime<=lsurvtime]
  ul<-ul[obstime<=lsurvtime]
  uhl<-uhl[obstime<=lsurvtime]
  obsertime<-obstimecopy[obstime<=lsurvtime]
  lsurvtime<-lsurvtime[obstime<=lsurvtime]
  cat(100*sum(cens)/n,"% experienced event\n")
  event.rate <- 100*sum(cens)/n
  fulldata<-data.frame(y=Y,IDL=idl ,hospitalid=hosidl, 
                       survtime=lsurvtime, obntime=obsertime, ctsxl= ctsxl, 
                       group=grpl, censoring=l.cens, u0=ul, u1=uhl)
  longdat<- fulldata
  # survdat<-  longdat[!duplicated(longdat$IDL),] 
  list(longdat=longdat, eventrate=event.rate)
}


"sort.dat"<-function(longdat,survdat){
  # Sorting survival & longitudinal data by survival time
  sort.long=matrix(0,dim(longdat)[1],dim(longdat)[2])
  index=longdat[,2]
  n.obs=diff(match(unique(index),index))
  n.obs[length(n.obs)+1]=length(index)-sum(n.obs)
  index=rep(survdat[,4],n.obs)
  sort.long=longdat[order(index),]
  sort.surv=matrix(0,dim(survdat)[1],dim(survdat)[2])
  sort.surv[,c(1:3,5:(dim(survdat)[2]))]=as.matrix(survdat[order(survdat[,4]),c(1:3,5:(dim(survdat)[2]))])
  sort.surv[,4]=sort(survdat[,4])
  list(long.s=data.frame(sort.long),surv.s=data.frame(sort.surv))
}
"longst"<-function(longdat){
  model.lme<- lme(y~obntime+group+ctsxl, random=list(hospitalid=~1, IDL=~1), control=lmeControl(returnObject=TRUE), data=longdat)
  sigma.z=model.lme$sigma^2
  v <- VarCorr(model.lme)
  sig.v <- as.numeric(v[4,1])
  sig.hv <- as.numeric(v[2,1])
  sig <-  matrix(0,2,2)
  sig[1,1]=sig.v
  sig[2,2]=sig.hv
  sig[1,2]=0
  sig[2,1]=0
  log.like=model.lme$logLik
  b1 <- (model.lme$coefficients$fixed)
  list(b1=data.frame(b1),sigma.z=data.frame(sigma.z),sigma.u=data.frame(sig),log.like=data.frame(log.like))
}

"survst" <- function(survdat) {
  n <- length(survdat[, 2])
  s <- survdat[,4]
  cen <- survdat[,8]
  if (cen[1]==0){ cen[1] = 1 }
  X2 <- as.matrix(survdat[,6:7])
  surv.start <- coxph(Surv(s, cen) ~ X2)
  surv.start.f <- survfit(surv.start)
  sf <- surv.start.f$time[surv.start.f$n.event != 0]
  nf <- length(sf)
  nev <- surv.start.f$n.event[surv.start.f$n.event != 0]
  haz <- coxph.detail(surv.start)$hazard	  
  rs <- rep(1:nf, c(diff(match(sf, s)), n + 1 - match(sf, s)[nf]))
  b2 <- coef(surv.start)
  names(b2) <- NULL
  noties <- (nf == sum(cen))
  
  ll <- surv.start$loglik[2] - sum(cen)
 # ll2 <- surv.start$loglik[2]
  
  # if (noties) {
  #   ll <- surv.start$loglik[2] - sum(cen)
  # } else {
  #   ll <- tiedsurvlike(survdat=survdat, b2=b2, haz=haz, sf=sf, rs=rs)
  # }
  list(b2 = b2, haz = haz, rs = rs, sf = sf, nev = nev, 
       log.like = ll)
}
tiedsurvlike <- function(survdat, b2, haz, sf, rs) {
  s <- survdat[,4]
  cen <- survdat[,8]
  if (cen[1] == 0) {
    cen[1] <-  1
  }
  X2 <- as.matrix(survdat[,6:7]) 
  bX <- X2 %*% b2
  f <- match(s, sf) * cen
  f[is.na(f)] <- 0
  chaz <- cumsum(haz)
  l1 <- log(haz[f]) + bX[cen == 1]
  l2 <- chaz[rs] * exp(bX)
  ll <- sum(l1) - sum(l2)
  ll
}

getD <- function(q, arg) {
  D <- matrix(0, q, length(arg))
  for (i in 1:q) {
    D[i, ] <- arg^(i - 1)
  }
  D
}

"em.alg"<-function(longdat,survdat,paraests,tol=0.001){
  index=longdat[,2]
  Y=longdat[,1]
  tt=longdat[,5]
  age<-longdat[,6]
  group<-longdat[,7]
  #X1=as.matrix(longdat[,4:dim(longdat)[2]])
  X879<- data.frame(rep(1, dim(longdat)[1]), tt, group, age)
  X1<- as.matrix(X879)
  n=length(survdat[,2])
  s=survdat[,4]
  cen=survdat[,8]
  p1=3 ##dim(X1)[2] (it is 3(age-group and time) only for this case)
  p2=2    #dim(survdat)[2]-3 (it is 2(age and group) only for this case)
  X2=0
  b1 <- paraests$b1[, 1]
  if(p2>0){X2=as.matrix(survdat[,6:7])}
  #corr=paraests$corr[,1]
  sig <- paraests$sigma.u
  haz <- paraests$haz
  sf <- paraests$sf
  id <- paraests$rs
  nev <- paraests$nev
  sigma.z=paraests$sigma.z[,1]
  n.obs=diff(match(unique(index),index))
  n.obs[length(n.obs)+1]=length(index)-sum(n.obs)
  N=sum(n.obs)
  maxn=max(n.obs)
  gpt=3
  ran=2
  lat=1
  b2 <- c(paraests$b2, rep(0, lat))
  ab=c(-1.224744870,0,1.224744870)
  w=c(0.295408975,1.1816359,0.295408975)
  gammat=matrix(0,gpt^2,ran)
  gammat[,1]=rep(ab,each=gpt)
  gammat[,2]=rep(ab,gpt)
  wvec=as.vector(w%*%t(w))
  EU=matrix(0,n,2)
  EUU=matrix(0,n,3)
  EexpU=matrix(0,n,length(haz))
  EU0expU=matrix(0,n,length(haz))
  EU1expU=matrix(0,n,length(haz))
  EU0U0expU=matrix(0,n,length(haz))
  EU0U1expU=matrix(0,n,length(haz))
  EU1U1expU=matrix(0,n,length(haz))
  W11=diag(maxn)
  W3=matrix(0,maxn,2)
  cvar=matrix(0,ran,ran)
  cvarch=matrix(0,ran,ran)
  major.it=20
  minor.it=10
  iter=0
  s1<-rep(1, length(s))
  sf1<-rep(1, length(sf))
  Ds1<-getD(ran, s1)
  Dst1<-t(Ds1)
  Dsf1<-getD(ran, sf1)
  cnn <- c(0, cumsum(n.obs))
  
  for (it in 1:major.it){
    for (it.2 in 1:minor.it){
      iter=iter+1
      W22=sig
      b2temp=0
      b2x=matrix(0,n,1)
      if(p2>0){b2temp=c(b2[1:p2])}
      if(p2>0){b2x=X2%*%b2temp}
      rlong=Y-X1%*%b1
      count=1
      nf<- length(sf)
      eb2x<-exp(b2x)
      for (i in 1:n){
        W21=matrix(0,2,n.obs[i])
        #rvec <- rlong[(cnn[i] + 1):cnn[i + 1]]
        
        rvec=rlong[count:(count+n.obs[i]-1)]
        W11=(sig[1,1]+2*sig[1,2]+sig[2,2])
        W21[1,1:n.obs[i]]=sig[1,1]+sig[1,2]
        W21[2,1:n.obs[i]]=sig[1,2]+sig[2,2]
        W11=W11+(sigma.z*diag(n.obs[i]))
        count=count+n.obs[i]
        W3=solve(W11,t(W21))
        cvar=W22-W21%*%W3
        cvar=cvar*2
        cvarch=chol(cvar)
        cvarch=(cvarch)
        # cm=t(W3)%*%rvec
        cm <- rvec%*%W3
        cmmat<-matrix(0,gpt^2,ran)
        cmmat[,1]<-rep(cm[1],gpt^2)
        cmmat[,2]<-rep(cm[2],gpt^2)
        #newumat=cvarch%*%t(gammat)+t(cmmat)
        newumat=gammat %*% cvarch+(cmmat)
        tnewumat <- t(newumat)
        egDUs<-exp((newumat)%*%(Dst1[i,]*b2[3]))
        
        eb2x<-exp(b2x)
        
        egDUsf<-exp((newumat)%*%(Dsf1[,1:id[i]]*b2[3]))
        ess<-exp(-(eb2x[i,]*egDUsf)%*%haz[1:id[i]])
        fvec<-egDUs*ess*wvec
        
        # fvec=exp(cen[i]*b2[p2+lat]*(newumat[1,]+newumat[2,]))
        #  ssvec=exp(b2[p2+lat]*(newumat[1,]+newumat[2,]))%*%haz[1:id[i]]
        #  ssvec=ssvec*exp(b2x[i,])
        #  fvec=fvec*wvec*exp(-ssvec)
        den=sum(fvec)
        EU[i,1:2]=t(tnewumat%*%fvec/den)
        EUU[i,1:2]=t((tnewumat^2)%*%fvec/den)
        EUU[i,3]=(tnewumat[1,]*tnewumat[2,])%*%fvec/den
        const=egDUsf[,1:id[i]]
        EexpU[i,1:id[i]]=t(fvec)%*%const/den
        EU0expU[i,1:id[i]]=t(fvec)%*%(tnewumat[1,]*const)/den
        EU1expU[i,1:id[i]]=t(fvec)%*%(tnewumat[2,]*const)/den
        EU0U0expU[i,1:id[i]]=t(fvec)%*%(tnewumat[1,]^2*const)/den
        EU0U1expU[i,1:id[i]]=t(fvec)%*%(tnewumat[1,]*tnewumat[2,]*const)/den
        EU1U1expU[i,1:id[i]]=t(fvec)%*%(tnewumat[2,]^2*const)/den
      }
      # M-STEP
      parac<-data.frame(c(b1,b2,sigma.z,sig))
      # UPDATE FOR BASELINE HAZARD
      SF=matrix(sf,n,length(sf),TRUE)
      S=matrix(rep(s,each=length(sf)),n,length(sf),TRUE)
      CEN=(S==SF)*cen
      CEN=colSums(CEN)
      RISK=(S>=SF)
      RISKe=colSums(RISK*EexpU*exp(b2x[,1]))
      haz=CEN/RISKe
      EUmat=matrix(0,N,2) 
      EUUmat=matrix(0,N,3)
      EUmat[,1]=rep(EU[,1],n.obs)
      EUmat[,2]=rep(EU[,2],n.obs)
      EUUmat[,1]=rep(EUU[,1],n.obs)
      EUUmat[,2]=rep(EUU[,2],n.obs)
      EUUmat[,3]=rep(EUU[,3],n.obs)
      EUsh=matrix(0,n,6)
      EUsh[,1]=rowSums(t(t(EexpU)*haz))
      EUsh[,2]=rowSums(t(t(EU0expU)*haz))
      EUsh[,3]=rowSums(t(t(EU1expU)*haz))
      EUsh[,4]=rowSums(t(t(EU0U0expU)*haz))
      EUsh[,5]=rowSums(t(t(EU1U1expU)*haz))
      EUsh[,6]=rowSums(t(t(EU0U1expU)*haz))
      # UPDATE BETA_1
      
      sum=EUmat[,1]+EUmat[,2]
      Ystar=Y-sum
      XTX=t(X1)%*%X1 
      XTY=t(X1)%*%Ystar
      b1=solve(XTX,XTY)
      # UPDATE NOISE VARIANCE
      bx=X1%*%b1
      r=Y-bx
      sum2=r^2-2*r*(sum)+EUUmat[,1]+(EUUmat[,2])
      sum2=sum2+2*EUUmat[,3]
      sigma.z=sum(sum2)/N
      # RANDOM EFFECTS COVARIANCE MATRIX
      sig[1,1]=sum(EUU[,1])/n
      sig[2,2]=sum(EUU[,2])/n
      sig[1,2]=sum(EUU[,3])/n
      sig[2,1]=sig[1,2]
      corr=sig[1,2]/sqrt(sig[1,1]*sig[2,2])
      # UPDATE BETA_2
      fd<-vector("numeric",p2+lat)
      sd<-matrix(0,p2+lat,p2+lat)
      eb2x=exp(b2x)
      if(p2>0){
        fd[1:p2]=c(colSums((cen*X2)-(X2*eb2x[,1]*EUsh[,1])))
        sd[1:p2,p2+lat]=-colSums((X2*eb2x[,1]*(EUsh[,2]+EUsh[,3])))
        sd=sd+t(sd)
        for (i in 1:p2){
          for (j in 1:p2){
            sd[i,j]=(-sum(X2[,i]*X2[,j]*eb2x*EUsh[,1]))}}}
      fd[p2+lat]=sum(cen*(EU[,1]+EU[,2]))-sum( eb2x*(EUsh[,2]+EUsh[,3]))
      sd[p2+lat,p2+lat]=(-sum(eb2x*(EUsh[,4]+2*EUsh[,6]+EUsh[,5])))
      b2=b2-solve(sd,fd)
    }
    para<-data.frame(c(b1,b2,sigma.z,sig))
    dd=abs(parac-para)
    if(max(dd)<tol){break}
  }
  if(it>major.it){print("Not converged")}
  list(b1=data.frame(b1),b2=data.frame(b2),sigma.z=data.frame(sigma.z),sigma.u=data.frame(sig),corr=data.frame(corr),haz=data.frame(haz))
}

jlike <- function(longdat, survdat, ran, likeests, lgpt=10) {
  lat <- dim(likeests$b2)[1]-2
  id=longdat[,2]
  Y=longdat[,1]
  tt=longdat[,5]
  age<-longdat[,6]
  group<-longdat[,7]
  #X1=as.matrix(longdat[,4:dim(longdat)[2]])
  X879<- data.frame(rep(1, dim(longdat)[1]), tt, group, age)
  X1<- as.matrix(X879)
  n=length(survdat[,2])
  s=survdat[,4]
  cen=survdat[,8]
  p1=3 ##dim(X1)[2] (it is 3(age-group and time) only for this case)
  p2=2    #dim(survdat)[2]-3 (it is 2(age and group) only for this case)
  X2=0
  ran <- 2
  lgpt <- 10
  if(p2>0){X2=as.matrix(survdat[,6:7])}
  nn <- diff(match(unique(id), id))
  nn[length(nn) + 1] <- length(id) - sum(nn)
  
  b1 <- likeests$b1[, 1]
  sigu <- likeests$sigma.u
  sigu <- as.matrix(sigu)
  vare <- (likeests$sigma.z[,1])
  b2 <- likeests$b2[, 1]
  haz <- likeests$haz
  haz <- as.matrix(haz)
  sf <- likeests$sf
  sf1<-rep(1, length(sf))
  s1<-rep(1, length(s))
  tt1<-rep(1, length(tt))
  rs <- likeests$rs
  N <- sum(nn)
  g <- gauss.quad.prob(lgpt, "normal", sigma = sqrt(0.5))
  ab <- g$nodes
  w <- g$weights * sqrt(pi)
  gmat <- matrix(0, lgpt^ran, ran)
  gmat[, 1] <- rep(ab, each = lgpt^(ran - 1))
  gmat[, 2] <- rep(ab, lgpt)
  w <- as.vector(w %x% w)
  l1 <- 0
  l2 <- 0
  r <- Y - X1 %*% b1
  Dtt <- getD(ran, tt1)
  Dtt2 <- t(Dtt)
  Ds <- getD(ran, s1)
  Dst <- t(Ds)
  Dsf <- getD(ran, sf1)
  cnn <- c(0, cumsum(nn))
  b2x <- X2 %*% b2[1:p2]
  if (p2 == 0) {
    b2x <- matrix(0, n, 1)
  }
  varei <- vare * diag(max(nn))
  cov <- sigu %*% Dtt
  tcov <- Dtt2 %*% sigu
  for (i in 1:n) {
    rv <- r[(cnn[i] + 1):cnn[i + 1]]
    ttv <- Dtt2[(cnn[i] + 1):cnn[i + 1], ]
    W21 <- cov[, (cnn[i] + 1):cnn[i + 1]]
    W12 <- tcov[(cnn[i] + 1):cnn[i + 1], ]
    W11 <- ttv %*% W21 + varei[1:nn[i], 1:nn[i]]
    if (nn[i] == 1) {
      W3 <- W12/as.vector(W11)
      cvch <- chol((sigu - tcrossprod(W21, W3)) *  2)
      cm <- matrix(W3 * rv, lgpt^ran, ran, TRUE)
    }    else {
      W3 <- solve(W11, W12)
      cvch <- chol((sigu - W21 %*% W3) * 2)
      cm <- matrix(rv %*% W3, lgpt^ran, ran, TRUE)
    }
    newu <- gmat %*% cvch + cm
    DUs <- newu %*% Ds[, i]
    DUsf <- newu %*% Dsf[, 1:rs[i]]
    ss <- exp(b2x[i, ]) * exp(newu %*% (Dsf[, 1:rs[i]] * b2[(p2 + 1):(p2 + lat)])) %*% haz[1:rs[i]]
    den <- sum(exp(cen[i] * (newu %*% (Dst[i, ] * b2[(p2 + 1):(p2 + lat)]) + b2x[i, ])) * (haz[rs[i]]^cen[i]) * w * exp(-ss))
    l2 <- l2 + 0
    if (den > 0) {
      l2 <- l2 + log(den)
    }
    l1 <- l1 - nn[i] * 0.5 * log(2 * pi) - 0.5 * log(det(W11)) - 
      0.5 * sum(rv * solve(W11, rv))
  }
  ll <- l1 + l2 - 0.5 * ran * n * log(pi)
  list(log.like = ll, longlog.like = l1, survlog.like = l2)
}

sampdata <- function(longdat){
  l.data<-longdat
  s.data<- l.data[!duplicated(l.data$IDL),] 
  n <- dim(s.data)[1]
  idSamp <- sample(s.data$IDL, n, replace = TRUE)
  sidSamp <- sort(idSamp)
  survdat22 <-s.data[sidSamp, ]
  survdat22$IDL <- 1:n
  mid <- match(sidSamp, l.data$IDL)
  rmid <- rep(mid, table(l.data$IDL)[sidSamp])
  nseq <- sequence(table(l.data$IDL)[sidSamp])
  longdat <- l.data[rmid + nseq - 1, ]
  nobs <- table(l.data$IDL)[sidSamp]
  longdat$IDL <- rep(1:n, nobs)
  row.names(longdat) <- NULL
  row.names(survdat22) <- NULL
  list(longitudinal = longdat, survival = survdat22, sidSamp=sidSamp)
}


bootstrap3.ml <- function(longdat, nboot){
  myfun2 <- function(i){
   sampledata <- sampdata(longdat)
   samp.long <- sampledata$longitudinal
   samp.survv <- sampledata$survival
   sort=sort.dat(samp.long,samp.survv)
   samp.long=(sort$long.s)
   samp.survv=(sort$surv.s)
   ldaests <-  longst(samp.long)
   survests <- survst(survdat=samp.survv)
   paraests=c(ldaests,survests)
   jointests=em.alg(samp.long,samp.survv,paraests)
 #  likeests <- c(jointests, list(rs = survests$rs, sf = survests$sf))
#   likelihoods <- jlike(samp.long, samp.survv, likeests,ran=2, lgpt=10)
  }
 # myfun <- function(i){x <- jointScore.sep(sampdata(l.data)$longitudinal)$U1}
  mycores <- detectCores() - 1
  test <<- mclapply(1:nboot, myfun2, mc.cores = mycores, mc.preschedule = FALSE)
  b.mat <- matrix(NA, ncol = 13, nrow=nboot)
  for (i in 1: nboot){
    b.mat[i,1:4] <- test[[i]]$b1$b1
    b.mat[i,5:7] <- test[[i]]$b2$b2
    b.mat[i, 8] <- test[[i]]$sigma.z$sigma.z
    b.mat[i,9] <- sqrt(test[[i]]$sigma.u[1,1])
    b.mat[i, 10:11] <- rep(test[[i]]$sigma.u[1,2], 2)
    b.mat[i, 12] <- sqrt(test[[i]]$sigma.u[2,2])
    b.mat[i, 13] <- test[[i]]$corr$corr
  }
  colnames(b.mat) <- c("long.int", "long.time", "long.cts", "long.group", "surv.cts", "surv.group", "assoc",
                                 "sigmaz", "sigmau0", "sigmau0uh", "sigmau0uh", "sigmauh", "corr")  #  print(TEST)
  vb <- apply(b.mat,2, var)
  vmean <- apply(b.mat, 2, mean)
  list(b.mat=b.mat, vb=vb, vmean =vmean)
}

options(scipen = 999) # get rid of e-values

# Start the clock!
ptm <- proc.time()

#resjoint.mlint<-function(nsim, m=25, noh=20, beta=c(7, -1, 0.5, 0.1), beta2=c(0.6,-0.6), var.v= matrix(c(0.3,0.1,0.1,0.2),2,2), var.hv=matrix(c(0.2,0.001,0.001,0.1),2,2), vare=0.1, gamma=c(0.4,0.7, 0.2,0.1), ntms=8, theta0=-2,theta1=0.2,censlam=exp(-3), censoring=T,  truncation=T){
  

resjoint.mlint<-function(nsim, m=25, noh=20,nboot=10,  beta=c(7, -1, 0.5, 0.1), beta2=c(0.6,-0.6), sig.v=0.4, sig.hv=0.12, sig.err=0.1, gamma=c(0.4,0.4), ntms=8, theta0=-2,censlam=exp(-3), censoring=T,  truncation=T){
#  result <- matrix(NA, nrow=nsim, ncol=13) #create a matrix to hold outcome
  listofdfs <- list()
  ev.rate <- vector("numeric", length = nsim)
  for(i in 1:nsim) {
    mldata<-simjointml(m=m, noh=noh, beta=beta, beta2=beta2, sig.v=sig.v, sig.hv=sig.hv, sig.err=sig.err, gamma=gamma, ntms=ntms, theta0=theta0, censlam=censlam, censoring=censoring, truncation=truncation)
  #  mldataintslope<-simjointmlintslope(m=m, noh=noh, beta=beta, beta2=beta2, var.v= var.v, var.hv=var.hv, vare=vare, gamma=gamma, ntms=ntms, theta0=theta0, theta1=theta1,censlam=censlam, censoring=censoring,  truncation=truncation)
    longdat <- mldata$longdat
    # survdat<- longdat[!duplicated(longdat$IDL),]
    # sort=sort.dat(longdat,survdat)
    # longdat=(sort$long.s)
    # survdat=(sort$surv.s)
    # ldaests=longst(longdat)
    # survests=survst(survdat)
    # paraests=c(ldaests,survests)
    # jointests=em.alg(longdat,survdat,paraests)
    # likeests <- c(jointests, list(rs = survests$rs, sf = survests$sf))
    # likelihoods <- jlike(longdat, survdat, likeests,ran=2, lgpt=10)
    # likelihoods.sep <- c(ldaests$log.like, survests$log.like)
    # l.est<-t(jointests[[1]])
    # s.est<-t(jointests[[2]])
    # sigmaz.est<-t(jointests[[3]])
    # sigmau.est<-t(jointests[[4]])
    # corr.est<-t(jointests[[5]])
    # result[i,]<-c(l.est, s.est,sigmaz.est,sigmau.est, corr.est)
        tryCatch({
        listofdfs[[i]] <- bootstrap3.ml(longdat, nboot = nboot)
        print(i)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    ev.rate[i] <- mldata$eventrate
  } 
  fitse.mean <- t(sapply(listofdfs, function(m) m$vmean))
  fitse.sd <- t(sapply(listofdfs, function(m) m$vb))
  # fitse.mean.all <- apply(fitse.mean, 2, mean)
  # fitse.sd.all <- apply(fitse.sd, 2, mean)
  # seandmean <- rbind(fitse.mean.all, fitse.sd.all)
  # rownames(seandmean) <- c("Estimate", "SE")
  # seandmean <- t(seandmean)
  list(listofdfs=listofdfs, fitse.mean=fitse.mean, fitse.sd=fitse.sd, event.rate=ev.rate)
}

finalresult70er <- resjoint.mlint(nboot=100, nsim=100, m=25, noh=20, gamma = rep(0.3,2), theta0 = -1.1)
finalresult20er <- resjoint.mlint(nboot=100, nsim=100, m=25, noh=20, gamma = rep(0.3,2), theta0 = -3)
stoptime<- proc.time() - ptm # Stop the clock
time.mlintboot <- stoptime[3]/60
print(time.mlintboot)

save.image("mlint bootres.RData")

#after loading the results, 
 gamma <- c(0.3,0.3)
 ort70 <- apply(finalresult70er$fitse.mean, 2, mean)
 sd70er <- apply(finalresult70er$fitse.sd, 2, mean)
 ort20 <- apply(finalresult20er$fitse.mean, 2, mean)
 sd20er <- apply(finalresult20er$fitse.sd, 2, mean)
 true.values <- c(beta, beta2, gamma[1], sig.err^2, sig.v,0,0,sig.hv, 0)
 tablo <- data.frame(true.values, ort70, sd70er, ort20, sd20er)
 
 mean(finalresult70er$event.rate)
 mean(finalresult20er$event.rate)
# m <- 25
# noh <- 20
# nsim <- 100
# beta=c(7, -1, 0.5, 0.1)
# beta2=c(0.6,-0.6)
# sig.v=0.4
# sig.hv=0
# vare=0.1
# gamma=0.3
# sigu <- c(sig.v, sig.hv)
# 
# n <- m*noh
# 
# 
# t05=qt(0.975,n-1)
# 
# res.mean <- finalresult$fitse.mean[,c(1,2,3,4,5,6,7,8,9,12, 10, 11, 13)]
# #res.sd <- finalresult$fitse.sd[,c(1,2,3,4,5,6,7,8,9,12, 10, 11, 13)]
# res.sd <- apply(finalresult$fitse.sd[,c(1,2,3,4,5,6,7,8,9,12, 10, 11, 13)] , 2, sqrt)
# 
# 
# joint.mean <-apply(res.mean, 2, mean)
# joint.sd <-apply(res.mean, 2, sd)
# coverage.b1 <- colSums((res.mean[,1:4]- t05*res.sd[,1:4]/sqrt(n)<=beta)&(res.mean[,1:4]+t05*res.sd[,1:4]/sqrt(n)>=beta))/nsim
# coverage.b2 <- colSums((res.mean[,5:6]- t05*res.sd[,5:6]/sqrt(n)<=beta2)&(res.mean[,5:6]+t05*res.sd[,5:6]/sqrt(n)>=beta2))/nsim
# coverage.assoc <- sum((res.mean[,7]- t05*res.sd[,7]/sqrt(n)<=gamma)&(res.mean[,7]+t05*res.sd[,7]/sqrt(n)>=gamma))/nsim
# coverage.sigu<- colSums((res.mean[,9:10]- t05*res.sd[,9:10]/sqrt(n)<=sigu)&(res.mean[,9:10]+t05*res.sd[,9:10]/sqrt(n)>=sigu))/nsim
# coverage.vare <- sum((res.mean[,8]- t05*res.sd[,8]/sqrt(n)<=vare)&(res.mean[,8]+t05*res.sd[,8]/sqrt(n)>=vare))/nsim
# biasb1 <- joint.mean[1:4]-beta#bias for long coefficients
# biasb2 <- joint.mean[5:6]-beta2#bias for surv coefficients
# bias.assoc <- joint.mean[7]-gamma#bias for gamma
# bias.sigu <- joint.mean[9:10]-sigu #bias for sigu
# bias.vare <- joint.mean[8]-vare
# mseb1 <- (joint.sd[1:4])^2+biasb1^2
# mseb2 <- (joint.sd[5:6])^2+biasb2^2
# mse.assoc <- (joint.sd[7])^2+bias.assoc^2
# mse.sigu <- (joint.sd[9:10])^2+bias.sigu^2
# mse.vare <- (joint.sd[8])^2+bias.vare^2
# true.values <- c(beta, beta2, gamma, sigu, vare)
# bias <- c(biasb1, biasb2, bias.assoc, bias.sigu, bias.vare)
# mse <- c(mseb1, mseb2, mse.assoc, mse.sigu, mse.vare)
# coverage <- c(coverage.b1, coverage.b2, coverage.assoc, coverage.sigu, coverage.vare)
# lastres <- data.frame(true.values, bias, mse, coverage)
# 
# 
# true.values2 <- c(beta, beta2, gamma, (vare)^2, c(0.4, 0,0, 0.12),NA )
# Estimates <- apply(t(sapply(finalresult$listofdfs, function(m) m$vmean)), 2, mean)
# SE <- apply(t(sapply(finalresult$listofdfs, function(m) m$vb)), 2, mean)
# ress <- t(rbind(true.values2, Estimates , SE))
# 

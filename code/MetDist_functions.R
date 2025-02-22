
#RAC function takes in a list of activity units and returns the list ranked
RAC <- function(x = ""){
  x = as.vector(x)
  x.ab = x[x >0]
  rac = x.ab[order(x.ab, decreasing = TRUE)]
  ranks <- as.vector(seq(1,length(rac)))
  rac <- cbind(ranks, rac)
  return(as.data.frame(rac))
}

Pareto <- function(x = "", pareto_para = ""){
  pareto <- fitdist(x[,1], "pareto", start = list(shape = 0.1, scale = 10))
  pareto_para <- rbind(pareto_para, c(pareto$estimate[1], pareto$estimate[2]))
  return(pareto_para)
}

Gamma <- function(x = "", gamma_para = ""){
  gamma <- fitdist(x[,1], "gamma", start = list(shape = 0.1, scale = 10))
  gamma_para <- rbind(gamma_para, c(gamma$estimate[1], gamma$estimate[2]))
  return(gamma_para)
}

Lognorm <- function(x = "", lnorm_para = ""){
  lnorm <- fitdist(x[,1], "lnorm")
  lnorm_para <- rbind(lnorm_para, c(lnorm$estimate[1], lnorm$estimate[2]))
  return(lnorm_para)
}

process <- function(x = "", channel = "", scale = ""){
  if(scale == TRUE){
    fcs <- flowCore::read.FCS(file = x, transformation = "scale")
  }
  else{
    fcs <- flowCore::read.FCS(file = x)
  }
  flow <- flowCore::exprs(fcs)
  RSG_H <- as.data.frame(flow[,channel])
  RSG_H <- unlist(RSG_H)
  RSG_ab = RSG_H[RSG_H >0]
  return(as.data.frame(RSG_ab))
}

BP_fxn <- function(CPMs, Samples, set){
  ##extract whole experiment info from top of csv
  
  #date experiment was run
  date <- as.Date(as.character(CPMs[1,2]), "%m/%d/%Y")
  
  
  #date the standard was produced
  date_std <- as.Date(as.character(CPMs[2,2]), "%m/%d/%Y")
  
  #DPM of the standard at date of production
  DPM_std <- as.double(as.character(CPMs[3,2]))
  
  #DPM of the standard based on scintillation counter on experiment date
  DPM_curr <- as.double(as.character(CPMs[4,2]))
  
  #half life of tritium - 12.346 years
  half_life <- as.double(as.character(CPMs[5,2]))
  
  #Mols of leucine in each tube based on hot leucine stock concentration
  M_Leu <- as.double(as.character(CPMs[6,2]))
  
  #CPMs of the voucher on experiment date
  Voucher <- as.double(as.character(CPMs[7,2]))
  
  ##remove whole experiment info from top of dataframe
  CPMs <- CPMs[-c(1:9),]
  colnames(CPMs) <- c("Sample", "CPM", "Kill")
  
  
  ##calculate time from the experiment date to the standard production date
  t <- as.numeric(date - date_std)/365
  
  ##calculate the expected DPMs of the standard based on t
  DPM_exp <- (DPM_std)*exp((-0.693/half_life)*t)
  
  ##calculate scintillation efficiency as DPM ratio
  efficiency <- DPM_curr/DPM_exp
  
  #divide CPMs into kill and sample dataframes
  Kills <- subset(CPMs, Kill == "T")
  CPMs <- subset(CPMs, Kill == "F")
  
  #convert CPMs to DPMs, DPMs = CPMs/efficiency
  CPMs$CPM <- as.numeric(as.character(CPMs$CPM))
  CPMs$DPM <- CPMs$CPM / efficiency
  
  Kills$CPM <- as.numeric(as.character(Kills$CPM))
  Kills$DPM <- Kills$CPM / efficiency
  
  #average DPMs for each sample and add to Tau
  for(x in unique(CPMs$Sample)){
    Samples[Samples$Sample == x, "DPM"] <- as.numeric(mean(CPMs[CPMs$Sample == x, "DPM"]))
  }
  
  
  #for each sample, subtract the corresponding kill DPM
  for (x in unique(Samples$Sample)){
    if(Samples[Samples$Sample == x, "Set"] == set){
      Samples[Samples$Sample == x, "DPMKills"] <- Samples[Samples$Sample ==x, "DPM"] - (as.numeric(Kills[Kills$Sample == x, "CPM"])/efficiency)
    }
  }
  
  #Determine Mols Leucine based on MLeu_sample = MLeu * DPM/voucher
  Samples$molLeu <- Samples$DPMKills * (1/2.22e12)*(1/161)*(1/1000)
  Samples$molLeuLhr <- Samples$molLeu * 1 * (1/0.0015)
  #Convert molLeu/L/hr to gC/L/hr
  Samples$gCLhr <- Samples$molLeuLhr * 131.2 * (1/0.073)* 0.86 * 2
  #Conversion of g to mol to umol to get umolC/hr
  Samples$uMChr <- Samples$gCLhr * 0.083333 * 100000
  Samples$log_uMChr <- log(Samples$uMChr, 10)
  
  return(Samples)
}

Met_fxn <- function(met_samples, Samples){
  for(x in unique(met_samples$Sample)){
    Samples$Shape[Samples$Sample == x] <- met_samples$shape[met_samples$Sample == x]
    Samples$Scale[Samples$Sample == x] <- met_samples$scale[met_samples$Sample == x]
  }
  
  return(Samples)
}

CDist <- function(x = ""){
  x <- as.vector(x)
  sum <- sum(x)
  rank <- 1
  total <- length(x)
  cdist <- as.vector((x[1]/sum)*100)
  Per <- as.vector((rank/total) * 100)
  act <- as.vector(x[1])
  for(num in x){
    rank <- rank + 1
    x <- x[-1]
    current <- cdist[length(cdist)] + ((x[1]/sum)*100)
    act <- c(act, x[1])
    Per <- c(Per, ((rank/total) *100))
    cdist <- c(cdist, current)
  }
  ranked <- cbind(cdist, Per, act)
  return(as.data.frame(ranked))
  
}

scientific_10 <- function(x) {
  x <- 10^x
  parse(text=gsub("e[+]?", " %*% 10^", scientific_format()(x)))
}


Evar<-function(x)   {
  v<-as.numeric(ncol(x)-1)
  for(k in 1:nrow(x)) {
    y<-x[k,2:ncol(x)]
    a<-as.numeric(ncol(x)-1); b<-a; d<-a
    for(i in 1:(ncol(x)-1)) {
      a[i]<-if(y[i]!=0) log(y[i]) else 0 }
    S<-sum(y>0)
    for(j in 1:(ncol(x)-1))    {
      b[j]<-a[[j]]/S}
    c<-sum(b)
    for(i in 1:(ncol(x)-1)) {
      d[i]<-if(y[i]!=0) ((a[[i]]-c)^2/S) else 0 }
    f<-sum(d)
    v[k]<-(1-2/pi*atan(f))   }
  v }



SimpE.cal <- function(x = ""){
  S <- rowSums(x > 0) * 1
  x = as.data.frame(x)
  D <- diversity(x, index = "inv")
  E <- (D)/S
  return(E)
}
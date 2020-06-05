### Delhaye, G., Hardy, O., Séleck, M., Ilunga wa Ilunga, E., Mahy, G., Meerts, P. (2019) 
### Plant community assembly along a natural metal gradient in central Africa: 
### functional and phylogenetic approach. Journal of Vegetation Science. 
## R code 

require(picante)
require(FD)
require(Hmisc)
require(AICcmodavg)
library(ade4)
library(vegan)
library(ape)
library(pez)
library(geiger)

###### Data 
## Traits data. Here we use the raw data although in the paper the log10 transformed data were
## used for the CWM (but not for the other indices as it can cause mathematical problems)
## This does not change the relationships between the indices and the soil Cu gradient
traits <- read.table("Traits.csv", sep = ";", dec = ".", header = T, row.names = 1)
df.fun.s <- as.data.frame(traits[,c(2,3,4,5,6,7,8,9,13,14)])
colnames(df.fun.s) <- c("VH","Prof", "LA",'LDMC',"LT","SLA","N","SM","Cu","Co")
head(df.fun.s)

## Community data (abundances)
abon <- read.table("Abundance.csv", sep = ";", dec = ".", header = T, row.names = 1)
comm.fun <- as.data.frame(abon)
head(comm.fun)

## Soil data. Here we use only the log 10 of the soil Cu content 
Soil <- read.table("Soil.csv", sep = ";", dec = ".", header = T, row.names = 1)
Cu<-log10(Soil$Cu) 

#######################################
### CWM and sesFRic
###################################
null.fd <- function(comm,trait,df,nreps) {
  
  # trait vector:
  df <- as.data.frame(df)
  df[,trait]=as.numeric( df[,trait])
  t=na.omit(df2vec(df,trait))
  t=t[names(t)%in%names(comm)]
  
  # Community matrix reduced to species for which the trait is known
  comm <- comm[,names(t)]
  comm <- comm/rowSums(comm) 
  
  # Number of community samples
  n_sites <- nrow(comm)
  
  ################ Results vector
  obs.cwms=NULL
  obs.frics = NULL
  mean.null.fric = NULL
  sd.null.fric = NULL
  
  # Looping on the n sites
  for (a in 1:n_sites) {
    com.obs=comm[a,comm[a,]>0] 
    t.obs=t[names(com.obs)]
    
    # Weighted mean	
    obs.cwms[a] = wtd.mean(t.obs,com.obs)   
    
    # Observed FRic
    fric.obs = max(t.obs)-min(t.obs)
    null_fric_array=NULL
    
    # looping on the nreps null samplings
    for (k in 1:nreps) {
      # Sample from the whole pool of species
      s.com=randomizeMatrix(comm[a,], null.model = "richness", iterations = 1) 
      s.com=s.com[1,s.com[1,]>0]
      t.obs=t[names(s.com)]
      null_fric_array[k] = max(t.obs)-min(t.obs)
    }
    
    # Null values+observed value
    null_fric_array=c(null_fric_array,fric.obs)
    obs.frics[a]=fric.obs
    mean.null.fric[a]= mean(null_fric_array)     
    sd.null.fric[a]= sd(null_fric_array)
  }
  
  # results data frame
  results=data.frame(round(obs.cwms,3),round(obs.frics,3),
                     round(mean.null.fric,3),round(sd.null.fric,3))
  names(results)=c("obs.cwm","obs.FRic","mean.null.FRic","sd.null.FRic")
  row.names(results)=row.names(comm)
  return(results)
}

#####################################################
##### Results CWM + FRic
height <- null.fd(comm.fun, "VH", df.fun.s, nreps=99)
leafarea <- null.fd(comm.fun, "LA", df.fun.s, nreps=99)
sla <- null.fd(comm.fun, "SLA", df.fun.s, nreps=99)
leafthick<- null.fd(comm.fun, "LT", df.fun.s, nreps=99)
ldmc <- null.fd(comm.fun, "LDMC", df.fun.s, nreps=99)
lnc<- null.fd(comm.fun, "N", df.fun.s, nreps=99)
lcoc<- null.fd(comm.fun, "Co", df.fun.s, nreps=99)
lcuc<- null.fd(comm.fun, "Cu", df.fun.s, nreps=99)
seedmass<- null.fd(comm.fun, "SM", df.fun.s, nreps=99)
rootdepth<- null.fd(comm.fun, "Prof", df.fun.s, nreps=99)

############################
## sesFRic
height.ses<-(height$obs.FRic-height$mean.null.FRic)/height$sd.null.FRic
leafarea.ses<-(leafarea$obs.FRic-leafarea$mean.null.FRic)/leafarea$sd.null.FRic 
sla.ses<-(sla$obs.FRic-sla$mean.null.FRic)/sla$sd.null.FRic 
leafthick.ses<-(leafthick$obs.FRic-leafthick$mean.null.FRic)/leafthick$sd.null.FRic
ldmc.ses<-(ldmc$obs.FRic-ldmc$mean.null.FRic)/ldmc$sd.null.FRic
lnc.ses<-(lnc$obs.FRic-lnc$mean.null.FRic)/lnc$sd.null.FRic
seedmass.ses<-(seedmass$obs.FRic-seedmass$mean.null.FRic)/seedmass$sd.null.FRic
rootdepth.ses<-(rootdepth$obs.FRic-rootdepth$mean.null.FRic)/rootdepth$sd.null.FRic 
lcuc.ses<-(lcuc$obs.FRic-lcuc$mean.null.FRic)/lcuc$sd.null.FRic 
lcoc.ses<-(lcoc$obs.FRic-lcoc$mean.null.FRic)/lcoc$sd.null.FRic 

####### Spearman correlations CWM
(cor.height.cwm<-rcorr(Cu,height$obs.cwm, type="spearman"))
(cor.leafarea.cwm<-rcorr(Cu,leafarea$obs.cwm, type="spearman"))
(cor.sla.cwm<-rcorr(Cu,sla$obs.cwm, type="spearman"))
(cor.leafthick.cwm<-rcorr(Cu,leafthick$obs.cwm, type="spearman"))
(cor.ldmc.cwm<-rcorr(Cu,ldmc$obs.cwm, type="spearman"))
(cor.lnc.cwm <-rcorr(Cu,lnc$obs.cwm, type="spearman"))
(cor.rootdepth.cwm<-rcorr(Cu,rootdepth$obs.cwm, type="spearman"))
(cor.seedmass.cwm<-rcorr(Cu,seedmass$obs.cwm, type="spearman"))
(cor.lcuc.cwm<-rcorr(Cu,lcuc$obs.cwm, type="spearman"))
(cor.lcoc.cwm<-rcorr(Cu,lcoc$obs.cwm, type="spearman"))

## Wilcoxon test FRic
wilcox.test(height.ses, mu=0, alternative = "two.sided")
wilcox.test(leafarea.ses, mu=0, alternative = "two.sided")
wilcox.test(sla.ses, mu=0, alternative = "two.sided")
wilcox.test(leafthick.ses, mu=0, alternative = "two.sided")
wilcox.test(ldmc.ses, mu=0, alternative = "two.sided")
wilcox.test(lnc.ses, mu=0, alternative = "two.sided")
wilcox.test(seedmass.ses, mu=0, alternative = "two.sided")
wilcox.test(rootdepth.ses, mu=0, alternative = "two.sided")
wilcox.test(lcuc.ses, mu=0, alternative = "two.sided")
wilcox.test(lcoc.ses, mu=0, alternative = "two.sided")

####### Spearman correlations FRic
(cor.height.ses<-rcorr(Cu,height.ses, type="spearman"))
plot(Cu,height.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.leafarea.ses<-rcorr(Cu,leafarea.ses, type="spearman"))
plot(Cu,leafarea.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.sla.ses<-rcorr(Cu,sla.ses, type="spearman"))
plot(Cu,sla.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.leafthick.ses<-rcorr(Cu,leafthick.ses, type="spearman"))
plot(Cu,leafthick.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.ldmc.ses<-rcorr(Cu,ldmc.ses, type="spearman"))
plot(Cu,ldmc.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.lnc.ses <-rcorr(Cu,lnc.ses, type="spearman"))
plot(Cu,lnc.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.rootdepth.ses<-rcorr(Cu,rootdepth.ses, type="spearman"))
plot(Cu,rootdepth.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.seedmass.ses<-rcorr(Cu,seedmass.ses, type="spearman"))
plot(Cu,seedmass.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.lcuc.ses<-rcorr(Cu,lcuc.ses, type="spearman"))
plot(Cu,lcuc.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.lcoc.ses<-rcorr(Cu,lcoc.ses, type="spearman"))
plot(Cu,lcoc.ses)
abline(h=0, col = "red", lwd = 3, lty =2)

################################################
######## FDis 
##########################
null.FDis <- function(comm,trait,df,nreps) {
  
  # trait vector:
  df <- as.data.frame(df)
  df[,trait]=as.numeric( df[,trait])
  t=na.omit(df2vec(df,trait))
  t=t[names(t)%in%names(comm)]
  
  # Community matrix reduced to species for which the trait is known
  comm <- comm[,names(t)]
  comm <- comm/rowSums(comm)  # transform abundance data in relative abundance data
  
  # Number of community samples
  n_sites <- nrow(comm)
  
  #Results vectors
  obs.fdiss=NULL
  mean.null.fdis=NULL
  sd.null.fdis=NULL
  
  # Looping on the n sites
  for (a in 1:n_sites) {
    com.obs=comm[a,comm[a,]>0]
    t.obs=t[names(com.obs)]
    
    # FD obs
    fd.obs <- dbFD(t.obs, com.obs, corr ="none", messages=F)
    fdis.obs = fd.obs$FDis
    null_fdis_array=NULL
    
    # looping on the nreps null samplings
    for (k in 1:nreps)
    {
      # Shuffle abundance keeping species and relative cover constant so richness does not change
      #Sample from the local community
      s.com=randomizeMatrix(com.obs, null.model = "richness", iterations = 1) 
      s.com=s.com[1,s.com[1,]>0]
      t.obs=t[names(s.com)]
      
      # FD for null communities
      fd.null <- dbFD(t.obs, s.com, corr = "none", messages=F)
      null_fdis_array[k] = fd.null$FDis
    }
    
    # Distributions nulles + valeurs obs:
    null_fdis_array=c(null_fdis_array,fdis.obs)
    obs.fdiss[a]=fdis.obs
    mean.null.fdis[a]= mean(null_fdis_array)     
    sd.null.fdis[a]= sd(null_fdis_array)
  }
  
  # results data frame
  results=data.frame(round(obs.fdiss,3),round(mean.null.fdis,3),round(sd.null.fdis,3)
  )
  names(results)=c("obs.FDis","mean.null.FDis","sd.null.FDis"
  )
  row.names(results)=row.names(comm)
  return(results)
}

###################################################"
##### Results FDis
height.FDis <- null.FDis(comm.fun, "VH", df.fun.s, nreps=99)
leafarea.FDis <- null.FDis(comm.fun, "LA", df.fun.s, nreps=99)
sla.FDis <- null.FDis(comm.fun, "SLA", df.fun.s, nreps=99)
leafthick.FDis<- null.FDis(comm.fun, "LT", df.fun.s, nreps=99)
ldmc.FDis <- null.FDis(comm.fun, "LDMC", df.fun.s, nreps=99)
lnc.FDis<- null.FDis(comm.fun, "N", df.fun.s, nreps=99)
lcoc.FDis<- null.FDis(comm.fun, "Co", df.fun.s, nreps=99)
lcuc.FDis<- null.FDis(comm.fun, "Cu", df.fun.s, nreps=99)
seedmass.FDis<- null.FDis(comm.fun, "SM", df.fun.s, nreps=99)
rootdepth.FDis<- null.FDis(comm.fun, "Prof", df.fun.s, nreps=99)

## sesFDis
height.FDis.ses<-(height.FDis$obs.FDis-height.FDis$mean.null.FDis)/height.FDis$sd.null.FDis
leafarea.FDis.ses<-(leafarea.FDis$obs.FDis-leafarea.FDis$mean.null.FDis)/leafarea.FDis$sd.null.FDis
sla.FDis.ses<-(sla.FDis$obs.FDis-sla.FDis$mean.null.FDis)/sla.FDis$sd.null.FDis
leafthick.FDis.ses<-(leafthick.FDis$obs.FDis-leafthick.FDis$mean.null.FDis)/leafthick.FDis$sd.null.FDis
ldmc.FDis.ses<-(ldmc.FDis$obs.FDis-ldmc.FDis$mean.null.FDis)/ldmc.FDis$sd.null.FDis
lnc.FDis.ses<-(lnc.FDis$obs.FDis-lnc.FDis$mean.null.FDis)/lnc.FDis$sd.null.FDis
seedmass.FDis.ses<-(seedmass.FDis$obs.FDis-seedmass.FDis$mean.null.FDis)/seedmass.FDis$sd.null.FDis
rootdepth.FDis.ses<-(rootdepth.FDis$obs.FDis-rootdepth.FDis$mean.null.FDis)/rootdepth.FDis$sd.null.FDis
lcuc.FDis.ses<-(lcuc.FDis$obs.FDis-lcuc.FDis$mean.null.FDis)/lcuc.FDis$sd.null.FDis 
lcoc.FDis.ses<-(lcoc.FDis$obs.FDis-lcoc.FDis$mean.null.FDis)/lcoc.FDis$sd.null.FDis

###### Test Wilcox FDis
wilcox.test(height.FDis.ses, mu=0, alternative = "two.sided")
wilcox.test(leafarea.FDis.ses, mu=0, alternative = "two.sided")
wilcox.test(sla.FDis.ses, mu=0, alternative = "two.sided")
wilcox.test(leafthick.FDis.ses, mu=0, alternative = "two.sided")
wilcox.test(ldmc.FDis.ses, mu=0, alternative = "two.sided")
wilcox.test(lnc.FDis.ses, mu=0, alternative = "two.sided")
wilcox.test(seedmass.FDis.ses, mu=0, alternative = "two.sided")
wilcox.test(rootdepth.FDis.ses, mu=0, alternative = "two.sided")
wilcox.test(lcuc.FDis.ses, mu=0, alternative = "two.sided")
wilcox.test(lcoc.FDis.ses, mu=0, alternative = "two.sided")

####### Spearman correlations FDis
(cor.height.FDis.ses<-rcorr(Cu,height.FDis.ses, type="spearman"))
plot(Cu,height.FDis.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.leafarea.FDis.ses<-rcorr(Cu,leafarea.FDis.ses, type="spearman"))
plot(Cu,leafarea.FDis.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.sla.FDis.ses<-rcorr(Cu,sla.FDis.ses, type="spearman"))
plot(Cu,sla.FDis.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.leafthick.FDis.ses<-rcorr(Cu,leafthick.FDis.ses, type="spearman"))
plot(Cu,leafthick.FDis.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.ldmc.FDis.ses<-rcorr(Cu,ldmc.FDis.ses, type="spearman"))
plot(Cu,ldmc.FDis.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.lnc.FDis.ses <-rcorr(Cu,lnc.FDis.ses, type="spearman"))
plot(Cu,lnc.FDis.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.rootdepth.FDis.ses<-rcorr(Cu,rootdepth.FDis.ses, type="spearman"))
plot(Cu,rootdepth.FDis.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.seedmass.FDis.ses<-rcorr(Cu,seedmass.FDis.ses, type="spearman"))
plot(Cu,seedmass.FDis.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.lcuc.FDis.ses<-rcorr(Cu,lcuc.FDis.ses, type="spearman"))
plot(Cu,lcuc.FDis.ses)
abline(h=0, col = "red", lwd = 3, lty =2)
(cor.lcoc.FDis.ses<-rcorr(Cu,lcoc.FDis.ses, type="spearman"))
plot(Cu,lcoc.FDis.ses)
abline(h=0, col = "red", lwd = 3, lty =2)

###########################################
#### Taxonomic changes and composition turn-over

#### Taxonomic richness
rich.tax <- rowSums(abon != 0)
cor.test(Cu, rich.tax, method="spearman")
plot(Cu, rich.tax, pch=1,cex=2, bty="l", ylim=c(18,36), ylab = "Species richness", xlab = expression('Soil Cu content (log10 mg.kg-1)'))

#### Composition turn over
abon.dist<-as.matrix(1-vegdist(abon, type="Euclidean", upper=T))
cu.dist <- as.matrix(vegdist(as.matrix(Soil$Cu), type="Euclidean",upper=T)) #Using raw Cu data to compute distance matrix
mantel(cu.dist, abon.dist)
abon.dist<-1-vegdist(abon, type="sorensen", upper=F)
cu.dist <- vegdist(as.matrix(Cu), type="Euclidean")
plot(cu.dist, abon.dist, pch=1,cex=2, ylim=c(0,0.85), bty="l", ylab = "Community Sørensen similarity", xlab = expression('Euclidean Cu dissimilarity'))

################################################
### Phylogenetic signal
#################################################
# Data
#Same abundances but with species name identical as phylogenetic tree
abon.p <- read.table("Abundance.phylo.csv", sep = ";", dec = ".", header = T, row.names = 1)
abon.p <- t(abon.p)

#Same traits but with species name identical as phylogenetic tree
traits <- read.table("Traits.phylo.csv", sep = ";", dec = ".", header = T, row.names = 1)
traits <- traits[,-c(1,10,11,12,15,16,17,18)]

## Phylogenetic tree
phylo <- read.tree("supertree.tre")
phylo <- multi2di(phylo, random = TRUE)#avoid polytomies

### ses PD
sespd <- ses.pd(abon.p, phylo, include.root = FALSE, runs =999)
mod.pd<-lmorigin(sespd$pd.obs.z~Cu, origin = FALSE, nperm =999)
(r.sq <- summary(mod.pd$reg)$r.squared)
(p.perm <- mod.pd$p.perm.t.2tail[2])
plot(Cu,sespd$pd.obs.z, xlab="Cu concentration (log10)", ylab = "z-value (PD)")
abline(h=0, col = "red", lwd = 3, lty =2)

############################################
#### Blomberg's K
tab.res <- as.data.frame(NULL)

for (f in 1:ncol(traits)){ #for all traits
  #select trait and remove na
   X <- traits[,f]
  names(X)<- rownames(traits)
  X <- na.omit(X)
  #compute the K and associated p
  tab.ind<- phylosignal(X[phylo$tip.label], phylo, reps = 999, checkdata = TRUE)
  tab.res<-rbind(tab.res,tab.ind)
  }
  rownames(tab.res) <- colnames(traits)
tab.res

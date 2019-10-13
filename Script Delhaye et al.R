### Delhaye, G., Hardy, O., Séleck, M., Ilunga wa Ilunga, E., Mahy, G., Meerts, P. (2019) 
### Plant community assembly along a natural metal gradient in central Africa: 
### functional and phylogenetic approach. Journal of Vegetation Science. 

## R code 

#######################################
### CWM and sesFRic
###################################

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
traits <- read.table("Traits.csv", sep = ";", dec = ".", header = T, row.names = 1)
raits$VH=sqrt(traits$HV)
traits$Prof=log(traits$Prof)
traits$LA=log(traits$LA)
traits$LDMC=sqrt(traits$LDMC)
traits$LT=sqrt(traits$LT)
traits$SLA=log(traits$SLA)
traits$N=log(traits$N)
traits$SM=log(traits$SM)
traits$Cu=log(traits$Cu)
traits$Co=log(traits$Co)
traits.s <- scale(traits)
colMeans(traits.s, na.rm =T)  
apply(traits.s, 2, sd ,na.rm =T)
df.fun.s <- as.data.frame(traits.s)

abon <- read.table("Abundance.csv", sep = ";", dec = ".", header = T, row.names = 1)
comm.fun <- as.data.frame(abon)
Cu <- read.table("Cu.csv", sep = ";", dec = ".", header = T, row.names = 1)
Cu<-log10(Cu)

#### Function FRic
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
    obs.cwm = wtd.mean(t.obs,com.obs)   
    
    # Observed FRic
    fric.obs = max(t.obs)-min(t.obs)
    null_fric_array=NULL
   
    # looping on the nreps null samplings
    for (k in 1:nreps)
    {
      s.com=randomizeMatrix(comm[1,], null.model = "richness", iterations = 1) 
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
  results=data.frame(round(obs.cwms,3),round(obs.frics,3),round(mean.null.fric,3),round(sd.null.fric,3)
                    )
    names(results)=c("obs.cwm","obs.FRic","mean.null.FRic","sd.null.FRic"
                  )
  row.names(results)=row.names(comm)
  return(results)
}

##### Results CWM + FRic
height <- null.fd(comm.fun, "HV", df.fun.s, nreps=9999)
leafarea <- null.fd(comm.fun, "LA", df.fun.s, nreps=9999)
sla <- null.fd(comm.fun, "SLA", df.fun.s, nreps=9999)
leafthick<- null.fd(comm.fun, "LT", df.fun.s, nreps=9999)
ldmc <- null.fd(comm.fun, "LDMC", df.fun.s, nreps=9999)
lnc<- null.fd(comm.fun, "N", df.fun.s, nreps=9999)
lcoc<- null.fd(comm.fun, "Copl", df.fun.s, nreps=9999)
lcuc<- null.fd(comm.fun, "Cupl", df.fun.s, nreps=9999)
seedmass<- null.fd(comm.fun, "SM", df.fun.s, nreps=9999)
rootdepth<- null.fd(comm.fun, "Prof", df.fun.s, nreps=9999)

## CWM
height.cwm<-(height$obs.cwm)
leafarea.cwm<-(leafarea$obs.cwm)
sla.cwm<-(sla$obs.cwm) 
leafthick.cwm<-(leafthick$obs.cwm)
ldmc.cwm<-(ldmc$obs.cwm)
lnc.cwm<-(lnc$obs.cwm)
seedmass.cwm<-(seedmass$obs.cwm)
rootdepth.cwm<-(rootdepth$obs.cwm)
lcuc.cwm<-(lcuc$obs.cwm)
lcoc.cwm<-(lcoc$obs.cwm)

####### Spearman correlations CWM
cor.height.cwm<-rcorr(Cu,height.cwm, type="spearman")
cor.leafarea.cwm<-rcorr(Cu,leafarea.cwm, type="spearman")
cor.sla.cwm<-rcorr(Cu,sla.cwm, type="spearman")
cor.leafthick.cwm<-rcorr(Cu,leafthick.cwm, type="spearman")
cor.ldmc.cwm<-rcorr(Cu,ldmc.cwm, type="spearman")
cor.lnc.cwm <-rcorr(Cu,lnc.cwm, type="spearman")
cor.rootdepth.cwm<-rcorr(Cu,rootdepth.cwm, type="spearman")
cor.seedmass.cwm<-rcorr(Cu,seedmass.cwm, type="spearman")
cor.lcuc.cwm<-rcorr(Cu,lcuc.cwm, type="spearman")
cor.lcoc.cwm<-rcorr(Cu,lcoc.cwm, type="spearman")

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

##### Test Wilcox FRic
wilcox.test(height.ses, zero, paired=TRUE)
wilcox.test(leafarea.ses, zero, paired=TRUE)
wilcox.test(sla.ses, zero, paired=TRUE)
wilcox.test(leafthick.ses, zero, paired=TRUE)
wilcox.test(ldmc.ses, zero, paired=TRUE)
wilcox.test(lnc.ses, zero, paired=TRUE)
wilcox.test(seedmass.ses, zero, paired=FALSE)
wilcox.test(rootdepth.ses, zero, paired=FALSE)
wilcox.test(lcuc.ses, zero, paired=FALSE)
wilcox.test(lcoc.ses, zero, paired=FALSE)

####### Spearman correlations FRic
cor.height.ses<-rcorr(Cu,height.ses, type="spearman")
cor.leafarea.ses<-rcorr(Cu,leafarea.ses, type="spearman")
cor.sla.ses<-rcorr(Cu,sla.ses, type="spearman")
cor.leafthick.ses<-rcorr(Cu,leafthick.ses, type="spearman")
cor.ldmc.ses<-rcorr(Cu,ldmc.ses, type="spearman")
cor.lnc.ses <-rcorr(Cu,lnc.ses, type="spearman")
cor.rootdepth.ses<-rcorr(Cu,rootdepth.ses, type="spearman")
cor.seedmass.ses<-rcorr(Cu,seedmass.ses, type="spearman")
cor.lcuc.ses<-rcorr(Cu,lcuc.ses, type="spearman")
cor.lcoc.ses<-rcorr(Cu,lcoc.ses, type="spearman")

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
height.FDis <- null.FDis(comm.fun, "HV", df.fun.s, nreps=9999)
leafarea.FDis <- null.FDis(comm.fun, "LA", df.fun.s, nreps=9999)
sla.FDis <- null.FDis(comm.fun, "SLA", df.fun.s, nreps=9999)
leafthick.FDis<- null.FDis(comm.fun, "LT", df.fun.s, nreps=9999)
ldmc.FDis <- null.FDis(comm.fun, "LDMC", df.fun.s, nreps=9999)
lnc.FDis<- null.FDis(comm.fun, "N", df.fun.s, nreps=9999)
lcoc.FDis<- null.FDis(comm.fun, "Copl", df.fun.s, nreps=9999)
lcuc.FDis<- null.FDis(comm.fun, "Cupl", df.fun.s, nreps=9999)
seedmass.FDis<- null.FDis(comm.fun, "SM", df.fun.s, nreps=9999)
rootdepth.FDis<- null.FDis(comm.fun, "Prof", df.fun.s, nreps=9999)

## sesFDis
height.FDis.ses<-(height.FDis$obs.FRic-height.FDis$mean.null.FRic)/height.FDis$sd.null.FRic
leafarea.FDis.ses<-(leafarea.FDis$obs.FRic-leafarea.FDis$mean.null.FRic)/leafarea.FDis$sd.null.FRic 
sla.FDis.ses<-(sla.FDis$obs.FRic-sla.FDis$mean.null.FRic)/sla.FDis$sd.null.FRic 
leafthick.FDis.ses<-(leafthick.FDis$obs.FRic-leafthick.FDis$mean.null.FRic)/leafthick.FDis$sd.null.FRic
ldmc.FDis.ses<-(ldmc.FDis$obs.FRic-ldmc.FDis$mean.null.FRic)/ldmc.FDis$sd.null.FRic
lnc.FDis.ses<-(lnc.FDis$obs.FRic-lnc.FDis$mean.null.FRic)/lnc.FDis$sd.null.FRic
seedmass.FDis.ses<-(seedmass.FDis$obs.FRic-seedmass.FDis$mean.null.FRic)/seedmass.FDis$sd.null.FRic
rootdepth.FDis.ses<-(rootdepth.FDis$obs.FRic-rootdepth.FDis$mean.null.FRic)/rootdepth.FDis$sd.null.FRic 
lcuc.FDis.ses<-(lcuc.FDis$obs.FRic-lcuc.FDis$mean.null.FRic)/lcuc.FDis$sd.null.FRic 
lcoc.FDis.ses<-(lcoc.FDis$obs.FRic-lcoc.FDis$mean.null.FRic)/lcoc.FDis$sd.null.FRic 

###### Test Wilcox FDis
wilcox.test(height.FDis.ses, zero, paired=TRUE)
wilcox.test(leafarea.FDis.ses, zero, paired=TRUE)
wilcox.test(sla.FDis.ses, zero, paired=TRUE)
wilcox.test(leafthick.FDis.ses, zero, paired=TRUE)
wilcox.test(ldmc.FDis.ses, zero, paired=TRUE)
wilcox.test(lnc.FDis.ses, zero, paired=TRUE)
wilcox.test(seedmass.FDis.ses, zero, paired=FALSE)
wilcox.test(rootdepth.FDis.ses, zero, paired=FALSE)
wilcox.test(lcuc.FDis.ses, zero, paired=FALSE)
wilcox.test(lcoc.FDis.ses, zero, paired=FALSE)

####### Spearman correlations FRic
cor.height.FDis.ses<-rcorr(Cu,height.FDis.ses, type="spearman")
cor.leafarea.FDis.ses<-rcorr(Cu,leafarea.FDis.ses, type="spearman")
cor.sla.FDis.ses<-rcorr(Cu,sla.FDis.ses, type="spearman")
cor.leafthick.FDis.ses<-rcorr(Cu,leafthick.FDis.ses, type="spearman")
cor.ldmc.FDis.ses<-rcorr(Cu,ldmc.FDis.ses, type="spearman")
cor.lnc.FDis.ses <-rcorr(Cu,lnc.FDis.ses, type="spearman")
cor.rootdepth.FDis.ses<-rcorr(Cu,rootdepth.FDis.ses, type="spearman")
cor.seedmass.FDis.ses<-rcorr(Cu,seedmass.FDis.ses, type="spearman")
cor.lcuc.FDis.ses<-rcorr(Cu,lcuc.FDis.ses, type="spearman")
cor.lcoc.FDis.ses<-rcorr(Cu,lcoc.FDis.ses, type="spearman")

################################################
### Phylogenetic signal
#################################################

abon.p <- t(comm.fun)
soil <- read.table("SolFung.csv", dec = ",", header = T, row.names = 1, sep = ";")
soil$Nbsp <- as.numeric(soil$Nbsp)
cu.p <- Cu
phylo <- read.tree("supertree.tre")
phylo <- multi2di(phylo, random = TRUE)

### ses PD
coph <- cophenetic(phylo)
sespd <- ses.pd(abon.p, phylo, include.root = FALSE, runs =9999)
mod.pd<-lmorigin(sespd$pd.obs.z~log10(cu.p$Cu), origin = FALSE, nperm =9999)
r.sq <- summary(mod.pd$reg)$r.squared
p.perm <- mod.pd$p.perm.t.2tail[2]
plot(log10(cu.p$Cu),sespd$pd.obs.z, xlab="Cu concentration (log10)", ylab = "z-value (PD)")
abline(lm(sespd$pd.obs.z~log10(cu.p$Cu)))
abline(h=0, col = "red", lwd = 3, lty =2)
legend("top", legend =c(paste("R² =", round(r.sq,2)),
                        paste("p value =", round(p.perm,4))),
       xpd=NA, inset = -0.1)

#### Computes Blomberg's K+ p
## Sort species to match tips labels
phylo$tip.label
traits.n<- cbind(traits, rownames(traits))
traits.n <- traits[match(phylo$tip.label, traits.n[,1]),] ## vary the parameter from first to last
traits.n<-traits.n[,-1]
colnames(traits.n)

## result table
tableau.res <- as.data.frame(NULL)

## choose species for which trait is present
X <- traits.n[,1]
names(X)<- rownames(traits.n)
X <- na.omit(X)

## Tree for which all species are known
list <-X
phyl.temp=add_species_to_genera(names(list),btree)
match.temp=match.phylo.comm(btree,t(list))
phylo.temp=match.temp$phy

## Compute index 
ordered.X <- X[phylo.temp$tip.label] #Sort the values of X following tips order
phylo.temp<-multi2di(phylo.temp)
tab.ind<-phylosignal(ordered.X, phylo.temp, reps = 9999, checkdata = TRUE)
tableau.res<-rbind(tableau.res,tab.ind)

rownames(tableau.res) <- colnames(traits.n)
write.table(tableau.res,"Phylo signal.txt")

###########################################
#### Taxonomic changes and composition turn-over

#### Taxonomic richness
rich.tax <- rowSums(abon != 0)
cor.test(cu, rich.tax, method="spearman")
plot(cu, rich.tax, pch=1,cex=2, bty="l", ylim=c(18,36), ylab = "Species richness", xlab = expression('Soil Cu content (log10 mg.kg-1)'), main=expression(paste(rho, " = -0.89***")))
legend("topleft", legend = c("a)"),cex = 1.3, text.font=2, box.lty=0,inset=c(-0.3,-0.3) )

#### Composition turn over
abon.dist<-as.matrix(1-vegdist(abon, type="Euclidean", upper=T))
cu.dist <- as.matrix(vegdist(as.matrix(Cu$Cu), type="Euclidean",upper=T))
mantel(cu.dist, abon.dist)
abon.dist<-1-vegdist(abon, type="sorensen", upper=F)
cu.dist <- vegdist(as.matrix(Cu$Cu), type="Euclidean")
plot(cu.dist, abon.dist, pch=1,cex=2, ylim=c(0,0.85), bty="l", ylab = "Community Sørensen similarity", xlab = expression('Euclidean Cu dissimilarity'), main=expression(paste(r, " = -0.82***")))
legend("topleft", legend = c("b)"),cex = 1.3, text.font=2, box.lty=0,inset=c(-0.3,-0.3) )


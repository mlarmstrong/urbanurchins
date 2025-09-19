setwd("~/Documents/IndividualProjects/Maddie/Ch3_popgen/")

library(rrBLUP)
library(SNPRelate)
library(LEA)
library(vroom)
library(ggplot)
library(viridis)

##Read in metadata
inds <- read.table("7.urbanurchfull.fam")
meta <- read.csv("urbanurchins_metadata_thinned.csv")
meta.order <- meta[match(inds$V1,meta$Sample_ID),]

#####Imputing genotypes - only need to do this once
#snpgdsBED2GDS("7.urbanurchfull.new.bed","7.urbanurchfull.new.fam","7.urbanurchfull.new.bim","7.urbanurchfull.new.gds",cvt.chr="char")
#genofile <- snpgdsOpen("7.urbanurchfull.new.gds")
#g <- read.gdsn(index.gdsn(genofile, "genotype"))
#g[g==3] <- 9
#write.table(t(g),"7.urbanurchfull.new.geno",sep="",row.names = F,col.names = F,quote = F)
#project.snmf = snmf("7.urbanurchfull.new.geno",K=1,
#                    entropy=T,repetitions=1,
#                    project="new")
#impute(project.snmf,"7.urbanurchfull.new.geno",method='mode')
#imp.gens <- vroom("7.urbanurchfull.new.lfmm_imputed.lfmm",delim=" ",col_names=1:ncol(g))
#snpgdsCreateGeno("7.urbanurchfull.new_imputed.gds",
#                 genmat=as.matrix(t(imp.gens)),
#                 sample.id=meta.order$Sample_ID,
#                 snp.id=1:ncol(g),
#                 snp.position=sites$Pos,
#                 snp.chromosome=sites$Chr)


##Read in site info
sites <- read.table("7.urbanurchfull.new.bim", header = FALSE,sep = "\t")[,c(1,4)]
names(sites) <- c("Chr","Pos")
sites$Snp <- paste(sites$Chr,"_",sites$Pos,sep="")

##write environment file
phen <- as.numeric(as.factor(meta.order$Dev))
#write.table(phen,file="Dev.env",row.names=F,col.names=F,quote=F)

##Read imputed genotypes
imp.g <- read.table("7.urbanurchfull.new.lfmm_imputed.lfmm",sep=" ")

###Run LFMM on all samples
lfmm <- lfmm2(input=as.matrix(imp.g),
              env=as.matrix(phen),
              K=2,
              effect.sizes=T)
pv <- lfmm2.test(object=lfmm,
                 input=as.matrix(imp.g),
                 env=as.matrix(phen),
                 linear=T,genomic.control=T)
plot(-log10(pv$pvalues))
q.vals <- p.adjust(pv$pvalues,method="fdr")
length(which(q.vals<0.1)) #no outliers

lfmm.frame <- data.frame(sites,pvalue=pv$pvalues,qvalue=q.vals)
write.table(lfmm.frame,"ubanurchin.lfmm.results.txt",quote=F,row.names=F)

#####################
###Polygenic scores
n=100
res.frame <- data.frame(LA.n=NA,
                        LA.u=NA,
                        SD.n=NA,
                        SD.u=NA,
                        VIC.n=NA,
                        VIC.u=NA,
                        p=NA)

##subset training set currently 60% (110 individuals)
for (i in 1:n) {
print(i)
train <- as.matrix(sample(1:nrow(meta.order),110))
test <- setdiff(1:nrow(meta.order),train)
train.g <- imp.g[train,]
train.e <- phen[train]
valid.g <- imp.g[test,]
valid.e <- phen[test]

##Run LFMM
lfmm <- lfmm2(input=as.matrix(train.g),
              env=as.matrix(train.e),
              K=2,
              effect.sizes=T)


## Effect sizes
b.values <- lfmm@B[,1]

#Prediction
pred <- scale(valid.g,scale=F) %*% matrix(b.values)

frame <- data.frame(Dev=meta.order$Dev[test],
                    Reg=meta.order$region[test],
                    Score=pred)

agg <- aggregate(frame$Score,list(frame$Dev,frame$Reg),mean,na.rm=T)
res.frame[i,1:6] <- agg$x
mod <- anova(lm(Score~Dev+Reg, data=frame))
res.frame[i,7] <- mod$`Pr(>F)`[1]
}

gath <- gather(res.frame,"metric","value")
write.csv(gath,"Polygenic_results.csv")

#How many bootstraps are significant?
length(which(res.frame$p<0.05))

#Make boxplot
gath.format <- gath %>% filter(metric!="p") %>% 
  separate(metric,into=c("Region","Dev"),sep="\\.") %>%
  mutate(Dev=case_when(Dev=="n"~"nonurban",
                       Dev=="u"~"urban"))
ggplot(gath.format,aes(x=Region,y=value,fill=Dev)) + 
  geom_boxplot() + theme_bw() +
  ylab("Polygenic Score") + scale_fill_viridis(discrete=T,begin=0.3,end=0.7)


##################
###Null for polygenic scores using randomized environment 
##################
n=100
res.frame <- data.frame(LA.n=NA,
                        LA.u=NA,
                        SD.n=NA,
                        SD.u=NA,
                        VIC.n=NA,
                        VIC.u=NA,
                        p=NA)

##subset training set currently 60% (110 individuals)
for (i in 1:n) {
  print(i)
  train <- as.matrix(sample(1:nrow(meta.order),110))
  test <- setdiff(1:nrow(meta.order),train)
  train2 <- sample(1:nrow(meta.order),110) #for randomizations
  test2 <- setdiff(1:nrow(meta.order),train2) #for randomizations
  train.g <- imp.g[train,]
  train.e <- phen[train2] #for randomizations
  valid.g <- imp.g[test,]
  valid.e <- phen[test2]
  
  ##Run LFMM
  lfmm <- lfmm2(input=as.matrix(train.g),
                env=as.matrix(train.e),
                K=2,
                effect.sizes=T)
  
  
  ## Effect sizes
  b.values <- lfmm@B[,1]
  
  ## Prediction
  pred <- scale(valid.g,scale=F) %*% matrix(b.values)
  
  frame <- data.frame(Dev=meta.order$Dev[test],
                      Reg=meta.order$region[test],
                      Score=pred)
  
  agg <- aggregate(frame$Score,list(frame$Dev,frame$Reg),mean,na.rm=T)
  res.frame[i,1:6] <- agg$x
  mod <- anova(lm(Score~Dev+Reg, data=frame))
  res.frame[i,7] <- mod$`Pr(>F)`[1]
}

#How many bootstraps are significant?
length(which(res.frame$p<0.05))
gath <- gather(res.frame,"metric","value")
write.csv(gath,"Polygenic_randomizations.csv")

#Make boxplot
gath.format <- gath %>% filter(metric!="p") %>% 
  separate(metric,into=c("Region","Dev"),sep="\\.") %>%
  mutate(Dev=case_when(Dev=="n"~"nonurban",
                       Dev=="u"~"urban"))
ggplot(gath.format,aes(x=Region,y=value,fill=Dev)) + 
  geom_boxplot() + theme_bw() +
  ylab("Polygenic Score") + scale_fill_viridis(discrete=T,begin=0.3,end=0.7)




##################################
####Playing with rrBLUP - don't use this
################################

##impute
impute = A.mat(g,impute.method="mean",return.imputed=T)
g.imp <- impute$imputed

##train
train=as.matrix(sample(1:nrow(meta.order),110))
test <- setdiff(1:nrow(meta.order),train)

pheno_train=meta.order[train,]
m_train=g.imp[train,]
pheno_valid=meta.order[test,]
m_valid=g.imp[test,]

answer <- mixed.solve(as.numeric(as.factor(pheno_train$Dev)),Z=m_train,K=NULL,
                      SE=FALSE,return.Hinv=FALSE)
e=as.matrix(answer$u)
pred_valid <- m_valid %*% e
pred <- (pred_valid[,1]) + answer$beta

boxplot(pred~pheno_valid$Dev+pheno_valid$region)
summary(lm(pred~pheno_valid$Dev*pheno_valid$region))




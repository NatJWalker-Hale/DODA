setwd("~/Analyses/2018/betalain_reconstruction/1pergenus/Caryophyllales_3524/round1")
library("phytools") # loads ape
library("phangorn")

tr <- read.tree("RAxML_bestTree.Caryophyllales_3524_outaln.cn_JFWT.rr.pl.tre.gennames.rot")
is.ultrametric(tr) # false
tr <- nnls.phylo(tr,cophenetic(tr),rooted = TRUE)
is.ultrametric(tr) # true
tr <- rotate(tr,649) # fixes Aizoaceae node

gen_pigment_df <- read.csv("../../../betalains_genuswise.csv",header=T)
intree <- gen_pigment_df[gen_pigment_df$Genus %in% tr$tip.label,]$Betalains
names(intree) <- as.character(gen_pigment_df[gen_pigment_df$Genus %in% tr$tip.label,]$Genus)
intree_not_df <- tr$tip.label[!(tr$tip.label %in% gen_pigment_df$Genus)]
intree_not_df_states <- c(rep(NA,length(intree_not_df)-82),rep(0,82)); names(intree_not_df_states) <- intree_not_df
pigment_states <- c(intree,intree_not_df_states)
pigment_states <- pigment_states[tr$tip.label]
pigment_states_df <- as.data.frame(cbind(names(pigment_states),as.integer(pigment_states)))

er_ace <- ace(pigment_states,tr,type="d",use.expm = T,use.eigen = F,model = "ER",marginal=F)
plot(tr,show.tip.label = F,type="f")
title(main="Marginal Reconstruction, Equal Rates")
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=er_ace$lik.anc,piecol=c("#0066ff","#ff3399"),cex=0.2)
tiplabels(pie =pigment_states[tr$tip.label],piecol=c("#ff3399","#0066ff"),cex=0.1)
AIC(er_ace)

ard_ace <- ace(pigment_states,tr,type="d",use.expm = T,use.eigen = F,model = "ARD",marginal=F)
plot(tr,show.tip.label = F,type="f")
title(main="Marginal Reconstruction, All Rates Different")
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=ard_ace$lik.anc,piecol=c("#0066ff","#ff3399"),cex=0.2)
tiplabels(pie =pigment_states[tr$tip.label],piecol=c("#ff3399","#0066ff"),cex=0.1)
AIC(ard_ace)

pigment_states_simmap <- matrix(data=NA,nrow=length(pigment_states),ncol= 2)
for (i in 1:length(pigment_states)) {
  if (is.na(pigment_states[i])) {
    pigment_states_simmap[i,] <- c(1,1)
  } else if (pigment_states[i] == 0) {
    pigment_states_simmap[i,] <- c(1,0) 
  } else if (pigment_states[i] == 1) {
    pigment_states_simmap[i,] <- c(0,1)
  }
}
rownames(pigment_states_simmap) <- names(pigment_states)
colnames(pigment_states_simmap) <- c("A","B")

#par(mfrow=c(1,2),mar=c(1,1,1,1))
#title(main="b)",adj=0,line=-6)
cols <- setNames(c("#0066ff","#ff3399"),c("A","B"))
er_simmap <- make.simmap(tr,pigment_states_simmap,model = "ER",nsim = 1000,pi=c(1.0,0.0))
# test nsim 
er_simmap_10k <- make.simmap(tr,pigment_states_simmap,model = "ER",nsim = 10000,pi=c(1.0,0.0))
desc_er_simmap <- describe.simmap(er_simmap)
desc_er_simmap_10k <- describe.simmap(er_simmap_10k)
###
#1000 trees with a mapped discrete character with states:
# A, B 
#
#trees have 11.394 changes between states on average
#
#changes are of the following types:
#  A,B   B,A
#x->y 4.973 6.421
#
#mean total time spent in each state is:
#  A            B    total
#raw  4825.9952056 8600.9249730 13426.92
#prop    0.3594268    0.6405732     1.00
###
plot(sample(er_simmap,1),colors = cols,type="fan",lwd=1,fsize=0.001)
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_er_simmap$ace,piecol=cols,cex=0.2)
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_er_simmap_10k$ace,piecol=cols,cex=0.2)
dens_er_simmap <- densityMap(er_simmap,plot=FALSE)
ard_simmap <- make.simmap(tr,pigment_states_simmap,model = "ARD",nsim = 1000,pi=c(1.0,0.0))
desc_ard_simmap <- describe.simmap(ard_simmap)
###
#1000 trees with a mapped discrete character with states:
#  A, B 
#
#trees have 8.62 changes between states on average
#
#changes are of the following types:
#  A,B   B,A
#x->y 6.774 1.846

#mean total time spent in each state is:
#  A            B    total
#raw  4891.8803879 8535.0397907 13426.92
#prop    0.3643338    0.6356662     1.00
###
plot(sample(er_simmap,1),colors = cols,type="fan",lwd=1,fsize=0.001)
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_ard_simmap$ace,piecol=cols,cex=0.2)

# polyMk stuff

pigment_states_polymk <- matrix(data=NA,nrow=length(pigment_states),ncol= 3)
for (i in 1:length(pigment_states)) {
  if (is.na(pigment_states[i])) {
    pigment_states_polymk[i,] <- c(1,1,1)
  } else if (pigment_states[i] == 0) {
    pigment_states_polymk[i,] <- c(1,0,0) 
  } else if (pigment_states[i] == 1) {
    pigment_states_polymk[i,] <- c(0,1,0)
  }
}
rownames(pigment_states_polymk) <- names(pigment_states)
colnames(pigment_states_polymk) <- c("A","B","A+B")

polymat <- matrix(c(0,0,1,0,0,2,3,4,0),ncol=3)
colnames(polymat) <- c("A","B","A+B")
rownames(polymat) <- c("A","B","A+B")

polymk_simmap <- make.simmap(tr,pigment_states_polymk,model = polymat,nsim = 1000,pi=c(1.0,0.0,0.0))
desc_polymk_simmap <- describe.simmap(polymk_simmap)
multicols <- setNames(c("#0066ff","purple","#ff3399"),c("A","A+B","B"))
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_polymk_simmap$ace,piecol=multicols,cex=0.2)

# plotting edges only 

edge_cols <- rep("lightgrey",nrow(tr$edge))
row_in_edge <- which(tr$edge[,2] %in% c(1:640))
#pigment_states[tr$edge[which(tr$edge[,2] %in% c(1:640)),2]]
#which(pigment_states[tr$edge[which(tr$edge[,2] %in% c(1:640)),2]] == 0)
edge_cols[row_in_edge[which(pigment_states[tr$edge[which(tr$edge[,2] %in% c(1:640)),2]] == 0)]] <- "#0066ff"
edge_cols[row_in_edge[which(pigment_states[tr$edge[which(tr$edge[,2] %in% c(1:640)),2]] == 1)]] <- "#ff3399"
plot.phylo(tr,show.tip.label = F,edge.color = edge_cols,type="f",y.lim = c(100,-100))

## for bea11 constraint - illustrative purposes

tr <- read.tree("RAxML_bestTree.Caryophyllales_3524_outaln.cn_bea11.rr.pl.tre.gennames.rot")
is.ultrametric(tr) # false
tr <- nnls.phylo(tr,cophenetic(tr),rooted = TRUE)
is.ultrametric(tr) # true
tr <- rotate(tr,650) # fixes Aizoaceae node

er_ace <- ace(pigment_states,tr,type="d",use.expm = T,use.eigen = F,model = "ER",marginal=F)
plot(tr,show.tip.label = F,type="f")
title(main="Marginal Reconstruction, Equal Rates")
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=er_ace$lik.anc,piecol=c("#0066ff","#ff3399"),cex=0.2)
tiplabels(pie =pigment_states[tr$tip.label],piecol=c("#ff3399","#0066ff"),cex=0.1)
AIC(er_ace) # 48.13557
 
ard_ace <- ace(pigment_states,tr,type="d",use.expm = T,use.eigen = F,model = "ARD",marginal=F)
#ard_ace_scale <- ace(pigment_states,tr,type="d",use.expm = T,use.eigen = F,model = "ARD")
# not used. Not a clue what is happening here, but it is sort of fucked. 
plot(tr,show.tip.label = F,type="f")
title(main="Marginal Reconstruction, All Rates Different")
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=ard_ace$lik.anc,piecol=c("#0066ff","#ff3399"),cex=0.2)
tiplabels(pie =pigment_states[tr$tip.label],piecol=c("#ff3399","#0066ff"),cex=0.1)
AIC(ard_ace) #49.85259

cols <- setNames(c("#0066ff","#ff3399"),c("A","B"))
er_simmap <- make.simmap(tr,pigment_states_simmap,model = "ER",nsim = 1000,pi=c(1.0,0.0))
er_simmap_5_5 <- make.simmap(tr,pigment_states_simmap,model = "ER",nsim = 1000,pi=c(0.5,0.5))
desc_er_simmap <- describe.simmap(er_simmap)
###
#1000 trees with a mapped discrete character with states:
#  A, B 
#
#trees have 10.397 changes between states on average
#
#changes are of the following types:
#  A,B   B,A
#x->y 5.495 4.902
#
#mean total time spent in each state is:
#  A            B    total
#raw  4680.3956263 5092.8182053 9773.214
#prop    0.4789004    0.5210996    1.000
###
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_er_simmap$ace,piecol=cols,cex=0.2)
desc_er_simmap_5_5 <- describe.simmap(er_simmap_5_5)
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_er_simmap_5_5$ace,piecol=cols,cex=0.2)
# prior sensitivity for ARD
ard_simmap_1_0 <- make.simmap(tr,pigment_states_simmap,model = "ARD",nsim = 1000,pi=c(1.0,0.0))
ard_simmap_9_1 <- make.simmap(tr,pigment_states_simmap,model = "ARD",nsim = 1000,pi=c(0.9,0.1))
ard_simmap_8_2 <- make.simmap(tr,pigment_states_simmap,model = "ARD",nsim = 1000,pi=c(0.8,0.2))
ard_simmap_7_3 <- make.simmap(tr,pigment_states_simmap,model = "ARD",nsim = 1000,pi=c(0.7,0.3))
ard_simmap_6_4 <- make.simmap(tr,pigment_states_simmap,model = "ARD",nsim = 1000,pi=c(0.6,0.4))
ard_simmap_5_5 <- make.simmap(tr,pigment_states_simmap,model = "ARD",nsim = 1000,pi=c(0.5,0.5))

desc_ard_simmap_1_0 <- describe.simmap(ard_simmap_1_0)
### 
#1000 trees with a mapped discrete character with states:
#  A, B 
#
#trees have 8.86 changes between states on average
#
#changes are of the following types:
#  A,B   B,A
#x->y 6.432 2.428
#
#mean total time spent in each state is:
#  A            B    total
#raw  4694.3459102 5078.8679215 9773.214
#prop    0.4803278    0.5196722    1.000
###
desc_ard_simmap_9_1 <- describe.simmap(ard_simmap_9_1)
desc_ard_simmap_8_2 <- describe.simmap(ard_simmap_8_2)
desc_ard_simmap_7_3 <- describe.simmap(ard_simmap_7_3)
desc_ard_simmap_6_4 <- describe.simmap(ard_simmap_6_4)
desc_ard_simmap_5_5 <- describe.simmap(ard_simmap_5_5)

nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_ard_simmap_1_0$ace,piecol=cols,cex=0.2)
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_ard_simmap_9_1$ace,piecol=cols,cex=0.2)
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_ard_simmap_8_2$ace,piecol=cols,cex=0.2)
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_ard_simmap_7_3$ace,piecol=cols,cex=0.2)
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_ard_simmap_6_4$ace,piecol=cols,cex=0.2)
nodelabels(node=1:tr$Nnode+Ntip(tr),pie=desc_ard_simmap_5_5$ace,piecol=cols,cex=0.2)

## let's investigate this a bit for rayDISC - prior sensitivity for ARD model 

## prep data for ray disc 

colnames(pigment_states_df) <- c("Species","Pigment")
pigment_states_df$Pigment <- as.character(pigment_states_df$Pigment)
pigment_states_df$Pigment[which(pigment_states_df$Pigment == 0)] <- "A"
pigment_states_df$Pigment[which(pigment_states_df$Pigment == 1)] <- "B"
pigment_states_df$Pigment[which(is.na(pigment_states_df$Pigment))] <- "A&B"
pigment_states_df$Pigment <- as.character(pigment_states_df$Pigment)

library(corHMM)
rd <- rayDISC(phy = tr,data = pigment_states_df,model = "ARD",root.p = c(1,0))
rd_5_5 <- rayDISC(phy = tr,data = pigment_states_df,model = "ARD",root.p = c(0.5,0.5))
rddef <- rayDISC(phy = tr,data = pigment_states_df,model = "ARD")


# R scripts for analysis of Arabidopsis genomes

library(dplyr)
library(seqinr)
library(GenomicRanges)
library(Biostrings)
library(segmentSeq)
library(Rsamtools)
library(GenomicAlignments)
library(parallel)
library(kmer)
library(ape)
library(phytools)
library(phangorn)
library(stats)
library(ade4)
library(ShortRead)
library(ips)
library(plot.matrix)
library(methimpute)

chrs <- c('Chr1','Chr2','Chr3','Chr4','Chr5')
chrss <- c('chr1','chr2','chr3','chr4','chr5')

colcen.seq <- read.fasta('ColCEN.fasta')
colcen.chrlen <- c(32537951,22117479,25741668,21576983,29477880,366924,154478)
colcen.cum <- c(0,cumsum(colcen.chrlen))
colcen.chrnames <- names(colcen.seq)
colcen.totlen <- sum(colcen.chrlen[1:5])

cen.rough.start <- c(14841110,3823792,13597188,4203902,11784131)
cen.rough.end <- c(17559778,6045243,15733925,6977949,14551809)

lerhifi.seq <- read.fasta(file="Ler-0_110x.ragtag_scaffolds.fa")
lerhifi.seq <- lerhifi.seq[1:5]
getLength(lerhifi.seq)
lerhifi.chrlen <- getLength(lerhifi.seq)[1:5]
lerhifi.cum <- c(0,cumsum(lerhifi.chrlen))
lerhifi.totlen <- sum(lerhifi.chrlen)

cos.f2 <- read.table('Table.CO_list.F2.v2.txt',header=T)
# 12410
cos.bc1.fem <- read.table(file='Table.CO_list.BC1.female.txt',header=T)
# 1009 9 
cos.bc1.male <- read.table(file='Table.CO_list.BC1.male.txt',header=T)
# 978 9
cos.all <- rbind(cos.f2,cos.bc1.male,cos.bc1.fem)
# 14397 9
cos.widths <- cos.all[,5]-cos.all[,4]
mean(cos.widths)
# 3992.111
cos.all <- cbind(cos.all,cos.widths)
cos.all <- cbind(cos.all,cos.all[,4]+(round(cos.widths/2)))
dim(cos.all)
# 14397 11
colnames(cos.all) <- c('Data_set','Sample_id','Chr_id','COs_start','COs_stop','Geno_start','Geno_stop',
                       'CO_num','Sex','cos.width','cos.mid')
hist(cos.widths,xlim=c(0,30000),breaks=1000)
abline(v=median(cos.widths),col='red',lty=2,lwd=0.5)

colcen.meth <- read.table(file='DNAmeth_Col_0_allcontexts_10kb_window_smoothN5_genome_10kb.tsv',header=T)
ler.meth <- read.table(file='DNAmeth_Ler_0_allcontexts_10kb_window_smoothN5_genome_10kb.tsv',header=T)

colcen.cenh3 <- read.table(file='Col_0_CenH3_ChIP_Col_0_CenH3_Input_MappedOn_t2t-col.20210610_lowXM_both_sort_norm_binSize10kblog2met1_col_unsmoothed.tsv',header=T)
ler.cenh3 <- read.table(file='Ler_0_CenH3_ChIP_Ler_0_CenH3_Input_MappedOn_Ler-0_110x.ragtag_scaffolds_lowXM_both_sort_norm_binSize10kb_unsmoothed.tsv',header=T)

all.hors <- read.csv(file='66_genomes_HORcounts.csv')
dim(all.hors)
lerhifi.cen178 <- all.hors[which(all.hors[,12]=="Ler-0_110x.ragtag_scaffolds.fa"),]
dim(lerhifi.cen178)
# 67954 13

colcen.cen178 <- read.table(file='ColCEN_CEN180.gff3',header=T)
colcen.cen178[which(colcen.cen178[,1]=='Col-CEN_chr1'),1] <- 'Chr1'
colcen.cen178[which(colcen.cen178[,1]=='Col-CEN_chr2'),1] <- 'Chr2'
colcen.cen178[which(colcen.cen178[,1]=='Col-CEN_chr3'),1] <- 'Chr3'
colcen.cen178[which(colcen.cen178[,1]=='Col-CEN_chr4'),1] <- 'Chr4'
colcen.cen178[which(colcen.cen178[,1]=='Col-CEN_chr5'),1] <- 'Chr5'
# 66131 18

colcen.athila <- read.table(file='ColCEN_ATHILA.gff3',header=F)
# 158 9

colcen.genes <- read.table('ColCEN_GENES_TAIR10.gff3',header=F)
colcen.genes <- colcen.genes[which(colcen.genes[,3]=='gene'),]
# 28782 9

colcen.edta <- read.table(file='t2t-col.20201227.fasta.mod.EDTA.TEanno.gff3')
# 42927 9
ler.edta <- read.table(file='Ler-0_110x.ragtag_scaffolds.fa.mod.EDTA.TEanno.gff3')
# 77549

syri.out <- read.table(file='Col-0_Ler-0.syri.out')
# 555865     12
syri.snps <- syri.out[which(syri.out[,11]=='SNP'),]
# 453714 
syri.syn <- syri.out[which(syri.out[,11]=='SYN'),]
# 214 12
syri.inv <- syri.out[which(syri.out[,11]=='INV'),]
# 214 12
syri.dup <- syri.out[which(syri.out[,11]=='DUP'),]
# 275 12
syri.hdr <- syri.out[which(syri.out[,11]=='HDR'),]
# 4654   12

qichao.snps <- read.table(file='Table.SNP_list.final.txt')
# 334680 4
length(qichao.snps[,1])/(colcen.totlen/1000)
# 2.546 SNPs/kb

#########################################
# define zones of crossover suppression #
#########################################

dim(cos.all)
# 14397    11
tot.cos <- length(cos.all[,1])
tot.chromatids <- 3613
one.cM <- tot.chromatids*0.01
genome.cMMb <- 3.027767

# express crossover data as cM/Mb

par(mfcol=c(5,1))
par(mar=c(2,2,2,2))
for(i in 1:5){
  cosall.chr <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=100000)
  wincolcen.cosall <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    wincolcen.cosall <- c(wincolcen.cosall,length(which(cosall.chr[,11]>colcen.wins[j]&cosall.chr[,11]<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.cosall/3613)
  win.cMMb <- win.cM/0.1
  plot(colcen.wins,c(0,win.cMMb),type='s',xlim=c(0,33000000),ylim=c(0,16))
  abline(h=3.027767,col=2,lty=2,lwd=0.5)
}

# define NRZ and middle of NRZ

chr.zero <- c(15000000,4000000,14000000,5000000,13000000)
chr.mid <- c(16000000,5000000,15000000,5000000,13000000)
all.left.nrz <- NULL
all.right.nrz <- NULL
for(i in 1:5){
  cosall.chr <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=100000)
  wincolcen.cosall <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    wincolcen.cosall <- c(wincolcen.cosall,length(which(cosall.chr[,11]>colcen.wins[j]&cosall.chr[,11]<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.cosall/3613)
  win.cMMb <- win.cM/0.1
  left.nrz <- colcen.wins[max(which(win.cMMb[which(colcen.wins<chr.zero[i])]>0))]
  right.wins <- colcen.wins[-which(colcen.wins<chr.zero[i])]
  right.nrz <- right.wins[min(which(win.cMMb[which(colcen.wins>chr.zero[i])]>0))]
  all.left.nrz <- c(all.left.nrz,left.nrz)
  all.right.nrz <- c(all.right.nrz,right.nrz)
}
all.mid.nrz <- all.left.nrz+((all.right.nrz-all.left.nrz)/2)

# define 1 cM boundaries for LRZ

par(mfcol=c(5,1))
par(mar=c(1.8,1.8,1.8,1.8))
all.left.lrz <- NULL
all.right.lrz <- NULL
for(i in 1:5){
  cosall.chr <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=100000)
  chr.cen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  chr.plus <- chr.cen178[which(chr.cen178[,7]=='+'),]
  chr.minus <- chr.cen178[which(chr.cen178[,7]=='-'),]
  wincolcen.cosall <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    wincolcen.cosall <- c(wincolcen.cosall,length(which(cosall.chr[,11]>colcen.wins[j]&cosall.chr[,11]<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.cosall/3613)
  win.cMMb <- win.cM/0.1
  dat <- cbind(colcen.wins,c(wincolcen.cosall,0),c(win.cMMb,0))
  dat.left <- dat[which(dat[,1]<all.mid.nrz[i]),]
  dat.right <- dat[which(dat[,1]>all.mid.nrz[i]),]
  revdat.left <- cbind(rev(dat.left[,1]),rev(dat.left[,2]),rev(dat.left[,3]))
  orient.dat <- cbind(revdat.left[1:48,],dat.right[1:48,])
  orient.dat <- cbind(orient.dat,cbind(orient.dat[,2]+orient.dat[,5]))
  orient.dat <- cbind(orient.dat,cumsum(orient.dat[,2]),cumsum(orient.dat[,5]),cumsum(orient.dat[,7]))
  colnames(orient.dat) <- c('left.coord','left.cos','left.cmmb','right.coord','right.cos','right.cmmb','lr.cos','left.cum','right.cum','lr.cum')
  onecM.left <- orient.dat[min(which(orient.dat[,8]>one.cM)),1]
  onecM.right <- orient.dat[min(which(orient.dat[,9]>one.cM)),4]
  print(i)
  plot(colcen.wins,c(win.cMMb,0),main=paste('cM/Mb Chr',i,sep=""),type='h',col='dark gray',ylim=c(-1,17))
  abline(h=genome.cMMb,col=1,lty=2,lwd=0.5)
  rug(all.mid.nrz[i],col=1,lty=1,lwd=1,side=3,ticksize=0.15)
  rug(c(all.left.nrz[i],all.right.nrz[i]),col='blue',lwd=1,side=3,ticksize=0.15)
  rug(c(onecM.left,onecM.right),col='purple',lwd=1,side=3,ticksize=0.15)
  rug(chr.plus[,4],col='red',ticksize=0.06)
  rug(chr.minus[,4],col='blue',ticksize=0.06)
  abline(v=cen.rough.start[i],col=1,lwd=1)
  abline(v=cen.rough.end[i],col=1,lwd=1)
  all.left.lrz <- c(all.left.lrz,onecM.left)
  all.right.lrz <- c(all.right.lrz,onecM.right)
}

# summary table of coordinates

all.dat <- NULL
for(i in 1:5){
  dat <- c(chrs[i],colcen.chrlen[i],all.left.lrz[i],all.left.nrz[i],all.mid.nrz[i],all.right.nrz[i],all.right.lrz[i])
  all.dat <- rbind(all.dat,dat)
}
all.dat <- as.data.frame(all.dat)
colnames(all.dat) <- c('Chr','chr.len','LRZ.left','NRZ.left','NRZ.mid','NRZ.right','LRZ.right')
cen.rough.start <- c(14841110,3823792,13597188,4203902,11784131)
cen.rough.end <- c(17559778,6045243,15733925,6977949,14551809)
all.dat <- cbind(all.dat,cen.rough.start,cen.rough.end)
write.csv(all.dat,file='Table_S1_LRZ_NRZ.csv')

all.left.lrz
# 13200001 2400001 11800001 1500001 10600001
all.right.lrz
# 18200001 8000001 17600001 8100001 16200001

all.left.nrz
# 14100001 3100001 13500001 3000001 11000001
all.right.nrz
# 17700001 6400001 16500001 7100001 14800001

lrz <- all.right.lrz-all.left.lrz
# 5000000 5600000 5800000 6600000 5600000
lrz.tot <- sum(lrz[1:5])
lrz.tot <- lrz.tot/1000000
# 28.6

nrz <- all.right.nrz-all.left.nrz
# 3600000 3300000 3000000 4100000 3800000
nrz.tot <- sum(nrz)
nrz.tot <- nrz.tot/1000000
# 17.8

lrz-nrz
# 1400000 2300000 2800000 2500000 1800000 
sum(lrz-nrz)/1000000
# 10.8

cen.rough.start <- c(14841110,3823792,13597188,4203902,11784131)
cen.rough.end <- c(17559778,6045243,15733925,6977949,14551809)

colcen.cen178.tally <- cen.rough.end-cen.rough.start
# 2718668 2221451 2136737 2774047 2767678
colcen.cen178.tot <- sum(colcen.cen178.tally)
colcen.cen178.tot <- colcen.cen178.tot/1000000
# 12.61858

lerhifi.cen178left <- c(14538003,3889715,13992335,5394546,12296069)
lerhifi.cen178right <- c(16448388,5401326,17050927,7765831,15662532)	
lerhifi.cen178.tally <- lerhifi.cen178right-lerhifi.cen178left
# 1910385 1511611 3058592 2371285 3366463
lerhifi.cen178.tot <- sum(lerhifi.cen178.tally)
lerhifi.cen178.tot <- lerhifi.cen178.tot/1000000
# 12.21834

colcen.chrlen <- colcen.chrlen[1:5]
colcen.arms <- colcen.chrlen-lrz 
# 27537951 16517479 19941668 14976983 23877880
colcen.arms.tot <- sum(colcen.arms)
colcen.arms.tot <- colcen.arms.tot/1000000
# 102.852

cenh3chip.col.mb <- c(22705001,11925001,21605001,12525001,20075001)
cenh3chip.ler.mb <- c(22025001,10845001,21815001,12545001,20345001)
mean(c(cenh3chip.col.mb,cenh3chip.ler.mb))
# 17641001
mean(c(colcen.cen178.tally,lerhifi.cen178.tally))
# 2483692

########################
# how many cos in LRZ? #
########################

lrz.tally <- NULL
for(i in 1:5){
  print(i)
  cos.chr <- cos.all[which(cos.all[,3]==chrs[i]),]
  sel.cos <- cos.chr[which(cos.chr[,11]>all.left.lrz[i]&cos.chr[,11]<all.right.lrz[i]),]
  print(dim(cos.chr))
  print(dim(sel.cos))
  lrz.tally <- c(lrz.tally,length(sel.cos[,1]))
}
# 67 118  71  78  66
sum(lrz.tally)
# 400

####################################
# calculate SNPs in LRZ, NRZ, arms #
####################################

chr.snp.dat <- NULL
for(i in 1:5){
  print(i)
  chr.snps <- qichao.snps[which(qichao.snps[,1]==chrs[i]),]
  chr.arm.left <- length(chr.snps[which(chr.snps[,2]>1&chr.snps[,2]<all.left.lrz[i]),1])/(all.left.lrz[i]/1000)
  chr.lrz.left <- length(chr.snps[which(chr.snps[,2]>all.left.lrz[i]&chr.snps[,2]<all.left.nrz[i]),1])/((all.left.nrz[i]-all.left.lrz[i])/1000)
  chr.nrz <- length(chr.snps[which(chr.snps[,2]>all.left.nrz[i]&chr.snps[,2]<all.right.nrz[i]),1])/((all.right.nrz[i]-all.left.nrz[i])/1000)
  chr.lrz.right <- length(chr.snps[which(chr.snps[,2]>all.right.nrz[i]&chr.snps[,2]<all.right.lrz[i]),1])/((all.right.lrz[i]-all.right.nrz[i])/1000)
  chr.arm.right <- length(chr.snps[which(chr.snps[,2]>all.right.lrz[i]&chr.snps[,2]<colcen.chrlen[i]),1])/((colcen.chrlen[i]-all.right.lrz[i])/1000)
  dat <- c(i,chr.arm.left,chr.lrz.left,chr.nrz,chr.lrz.right,chr.arm.right)
  chr.snp.dat <- cbind(chr.snp.dat,dat)
}
rownames(chr.snp.dat) <- c('Chr','arm.left','lrz.left','nrz','lrz.right','arm.right')     
write.csv(chr.snp.dat,file='SNP.density.csv')

#################################################################
# Figure S1 compare crossovers and SNPs, with NRZ LRZ indicated #
#################################################################

par(mfcol=c(5,3))
par(mar=c(1.8,1.8,1.8,1.8))
for(i in 1:5){
  cosall.chr <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=100000)
  chr.cen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  chr.plus <- chr.cen178[which(chr.cen178[,7]=='+'),]
  chr.minus <- chr.cen178[which(chr.cen178[,7]=='-'),]
  wincolcen.cosall <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    wincolcen.cosall <- c(wincolcen.cosall,length(which(cosall.chr[,11]>colcen.wins[j]&cosall.chr[,11]<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.cosall/3613)
  win.cMMb <- win.cM/0.1
  plot(colcen.wins,c(win.cMMb,0),main=paste('cM/Mb Chr',i,sep=""),type='h',col='dark gray',ylim=c(-1,17))
  abline(h=genome.cMMb,col=1,lty=2,lwd=0.5)
  rug(all.mid.nrz[i],col=1,lty=1,lwd=1,side=3,ticksize=0.15)
  rug(c(all.left.nrz[i],all.right.nrz[i]),col='blue',lwd=1,side=3,ticksize=0.15)
  rug(c(all.left.lrz[i],all.right.lrz[i]),col='purple',lwd=1,side=3,ticksize=0.15)
  abline(v=all.mid.nrz[i],col=1,lty=1,lwd=1)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col='blue')
  abline(v=c(onecM.left,onecM.right),col='purple')
  rug(chr.plus[,4],col='red',ticksize=0.06)
  rug(chr.minus[,4],col='blue',ticksize=0.06)
}
for(i in 1:5){
  print(i)
  chr.cen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  chr.plus <- chr.cen178[which(chr.cen178[,7]=='+'),]
  chr.minus <- chr.cen178[which(chr.cen178[,7]=='-'),]
  print(i)
  chr.snps <- qichao.snps[which(qichao.snps[,1]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=100000)
  wincolcen.snps <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    wincolcen.snps <- c(wincolcen.snps,length(which(chr.snps[,2]>colcen.wins[j]&chr.snps[,2]<colcen.wins[j+1])))
  }
  mean.snps <- mean(wincolcen.snps)
  plot(colcen.wins,c(wincolcen.snps,0),main=paste('SNPs Chr',i,sep=""),type='h',col='gray',ylim=c(-80,1000))
  abline(h=mean.snps,col=1,lty=2,lwd=0.5)
  rug(all.mid.nrz[i],col=1,lty=1,lwd=1,side=3,ticksize=0.15)
  rug(c(all.left.nrz[i],all.right.nrz[i]),col='blue',lwd=1,side=3,ticksize=0.15)
  rug(c(all.left.lrz[i],all.right.lrz[i]),col='purple',lwd=1,side=3,ticksize=0.15)
  abline(v=all.mid.nrz[i],col=1,lty=1,lwd=1)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col='blue')
  abline(v=c(onecM.left,onecM.right),col='purple')
  rug(chr.plus[,4],col='red',ticksize=0.06)
  rug(chr.minus[,4],col='blue',ticksize=0.06)
}
for(i in 1:5){
  cosall.chr <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=100000)
  chr.cen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  chr.plus <- chr.cen178[which(chr.cen178[,7]=='+'),]
  chr.minus <- chr.cen178[which(chr.cen178[,7]=='-'),]
  chr.snps <- qichao.snps[which(qichao.snps[,1]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=100000)
  wincolcen.cosall <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    wincolcen.cosall <- c(wincolcen.cosall,length(which(cosall.chr[,11]>colcen.wins[j]&cosall.chr[,11]<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.cosall/3613)
  win.cMMb <- win.cM/0.1
  wincolcen.snps <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    wincolcen.snps <- c(wincolcen.snps,length(which(chr.snps[,2]>colcen.wins[j]&chr.snps[,2]<colcen.wins[j+1])))
  }
  mean.snps <- mean(wincolcen.snps)
  print(i)
  plot(colcen.wins,c(wincolcen.snps,0),type='h',col='blue',ylim=c(0,1000))
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),main=paste('cM/Mb Chr',i,sep=""),type='h',col='red',ylim=c(0,17),yaxt='n')
  axis(side=4)
  rug(all.mid.nrz[i],col=1,lty=1,lwd=1,side=3,ticksize=0.15)
  rug(c(all.left.nrz[i],all.right.nrz[i]),col='blue',lwd=1,side=3,ticksize=0.15)
  rug(c(all.left.lrz[i],all.right.lrz[i]),col='purple',lwd=1,side=3,ticksize=0.15)
  rug(chr.plus[,4],col='dark gray',ticksize=0.01)
  rug(chr.minus[,4],col='black',ticksize=0.01)
}

################################################################
# compare crossovers and SNPs, and TEs, with NRZ LRZ indicated #
################################################################

par(mfcol=c(1,5))
par(mar=c(1.8,1.8,1.8,1.8))
for(i in 1:5){
  cosall.chr <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=100000)
  chr.cen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  chr.plus <- chr.cen178[which(chr.cen178[,7]=='+'),]
  chr.minus <- chr.cen178[which(chr.cen178[,7]=='-'),]
  chr.edta <- colcen.edta[which(colcen.edta[,1]==chrs[i]),]
  copia <- chr.edta[which(chr.edta[,3]=='Copia_LTR_retrotransposon'),]
  copia.width <- copia[,5]-copia[,4]
  gypsy <- chr.edta[which(chr.edta[,3]=='Gypsy_LTR_retrotransposon'),]
  gypsy.width <- gypsy[,5]-gypsy[,4]
  linel1 <- chr.edta[which(chr.edta[,3]=='L1_LINE_retrotransposon'),]
  linel1.width <- linel1[,5]-linel1[,4]
  count.copia <- NULL
  count.gypsy <- NULL
  count.linel1 <- NULL
  bp.copia <- NULL
  bp.gypsy <- NULL
  bp.linel1 <- NULL
  wincolcen.cosall <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    wincolcen.cosall <- c(wincolcen.cosall,length(which(cosall.chr[,11]>colcen.wins[j]&cosall.chr[,11]<colcen.wins[j+1])))
    bp.copia <- c(bp.copia,sum(copia.width[which(copia[,4]>colcen.wins[j]&copia[,4]<colcen.wins[j+1])]))
    bp.gypsy <- c(bp.gypsy,sum(gypsy.width[which(gypsy[,4]>colcen.wins[j]&gypsy[,4]<colcen.wins[j+1])]))
    bp.linel1 <- c(bp.linel1,sum(linel1.width[which(linel1[,4]>colcen.wins[j]&linel1[,4]<colcen.wins[j+1])]))
    count.copia <- c(count.copia,length(which(copia[,4]>colcen.wins[j]&copia[,4]<colcen.wins[j+1])))
    count.gypsy <- c(count.gypsy,length(which(gypsy[,4]>colcen.wins[j]&gypsy[,4]<colcen.wins[j+1])))
    count.linel1 <- c(count.linel1,length(which(linel1[,4]>colcen.wins[j]&linel1[,4]<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.cosall/3613)
  win.cMMb <- win.cM/0.1
  plot(colcen.wins,c(bp.gypsy,0),type='h',col='black',ylim=c(-5000,100000),main='RNA transposons')
  par(new=T)
  plot(colcen.wins,c(bp.copia,0),type='h',col='darkorchid1',ylim=c(-5000,100000))
  par(new=T)
  plot(colcen.wins,c(bp.linel1,0),type='h',col='green',ylim=c(-5000,100000))
  rug(all.mid.nrz[i],col=1,lty=1,lwd=1,side=3,ticksize=0.05)
  rug(c(all.left.nrz[i],all.right.nrz[i]),col='blue',lwd=1,side=3,ticksize=0.05)
  rug(c(all.left.lrz[i],all.right.lrz[i]),col='purple',lwd=1,side=3,ticksize=0.05)
  rug(chr.plus[,4],col='red',ticksize=0.06)
  rug(chr.minus[,4],col='blue',ticksize=0.06)
}

par(mfcol=c(5,2))
par(mar=c(1.8,1.8,1.8,1.8))
for(i in 1:5){
  cosall.chr <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=100000)
  chr.cen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  chr.plus <- chr.cen178[which(chr.cen178[,7]=='+'),]
  chr.minus <- chr.cen178[which(chr.cen178[,7]=='-'),]
  chr.edta <- colcen.edta[which(colcen.edta[,1]==chrs[i]),]
  copia <- chr.edta[which(chr.edta[,3]=='Copia_LTR_retrotransposon'),]
  copia.width <- copia[,5]-copia[,4]
  gypsy <- chr.edta[which(chr.edta[,3]=='Gypsy_LTR_retrotransposon'),]
  gypsy.width <- gypsy[,5]-gypsy[,4]
  linel1 <- chr.edta[which(chr.edta[,3]=='L1_LINE_retrotransposon'),]
  linel1.width <- linel1[,5]-linel1[,4]
  count.copia <- NULL
  count.gypsy <- NULL
  count.linel1 <- NULL
  bp.copia <- NULL
  bp.gypsy <- NULL
  bp.linel1 <- NULL
  wincolcen.cosall <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    wincolcen.cosall <- c(wincolcen.cosall,length(which(cosall.chr[,11]>colcen.wins[j]&cosall.chr[,11]<colcen.wins[j+1])))
    bp.copia <- c(bp.copia,sum(copia.width[which(copia[,4]>colcen.wins[j]&copia[,4]<colcen.wins[j+1])]))
    bp.gypsy <- c(bp.gypsy,sum(gypsy.width[which(gypsy[,4]>colcen.wins[j]&gypsy[,4]<colcen.wins[j+1])]))
    bp.linel1 <- c(bp.linel1,sum(linel1.width[which(linel1[,4]>colcen.wins[j]&linel1[,4]<colcen.wins[j+1])]))
    count.copia <- c(count.copia,length(which(copia[,4]>colcen.wins[j]&copia[,4]<colcen.wins[j+1])))
    count.gypsy <- c(count.gypsy,length(which(gypsy[,4]>colcen.wins[j]&gypsy[,4]<colcen.wins[j+1])))
    count.linel1 <- c(count.linel1,length(which(linel1[,4]>colcen.wins[j]&linel1[,4]<colcen.wins[j+1])))
  }
  mean.cos <- mean(wincolcen.cosall)
  plot(colcen.wins,c(bp.gypsy,0),type='h',col='black',ylim=c(-5000,100000),main='RNA transposons')
  par(new=T)
  plot(colcen.wins,c(bp.copia,0),type='h',col='darkorchid1',ylim=c(-5000,100000))
  par(new=T)
  plot(colcen.wins,c(bp.linel1,0),type='h',col='green',ylim=c(-5000,100000))
  rug(all.mid.nrz[i],col=1,lty=1,lwd=1,side=3,ticksize=0.05)
  rug(c(all.left.nrz[i],all.right.nrz[i]),col='blue',lwd=1,side=3,ticksize=0.05)
  rug(c(all.left.lrz[i],all.right.lrz[i]),col='purple',lwd=1,side=3,ticksize=0.05)
  rug(chr.plus[,4],col='red',ticksize=0.06)
  rug(chr.minus[,4],col='blue',ticksize=0.06)
}
for(i in 1:5){
  cosall.chr <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=100000)
  chr.cen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  chr.plus <- chr.cen178[which(chr.cen178[,7]=='+'),]
  chr.minus <- chr.cen178[which(chr.cen178[,7]=='-'),]
  chr.edta <- colcen.edta[which(colcen.edta[,1]==chrs[i]),]
  copia <- chr.edta[which(chr.edta[,3]=='Copia_LTR_retrotransposon'),]
  copia.width <- copia[,5]-copia[,4]
  gypsy <- chr.edta[which(chr.edta[,3]=='Gypsy_LTR_retrotransposon'),]
  gypsy.width <- gypsy[,5]-gypsy[,4]
  linel1 <- chr.edta[which(chr.edta[,3]=='L1_LINE_retrotransposon'),]
  linel1.width <- linel1[,5]-linel1[,4]
  heli <- chr.edta[which(chr.edta[,3]=='helitron'),]
  heli.width <- heli[,5]-heli[,4]
  mu <- chr.edta[which(chr.edta[,3]=='Mutator_TIR_transposon'),]
  mu.width <- mu[,5]-mu[,4]
  cacta <- chr.edta[which(chr.edta[,3]=='CACTA_TIR_transposon'),]
  cacta.width <- cacta[,5]-cacta[,4]
  pif <- chr.edta[which(chr.edta[,3]=='PIF_Harbinger_TIR_transposon'),]
  pif.width <- pif[,5]-pif[,4]
  hat <- chr.edta[which(chr.edta[,3]=='hAT_TIR_transposon'),]
  hat.width <- hat[,5]-hat[,4]
  tc1 <- chr.edta[which(chr.edta[,3]=='Tc1_Mariner_TIR_transposon'),]
  tc1.width <- tc1[,5]-tc1[,4]
  count.copia <- NULL
  count.gypsy <- NULL
  count.linel1 <- NULL
  bp.copia <- NULL
  bp.gypsy <- NULL
  bp.linel1 <- NULL
  bp.heli <- NULL
  bp.mu <- NULL
  bp.cacta <- NULL
  bp.pif <- NULL
  bp.hat <- NULL
  bp.tc1 <- NULL
  wincolcen.cosall <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    wincolcen.cosall <- c(wincolcen.cosall,length(which(cosall.chr[,11]>colcen.wins[j]&cosall.chr[,11]<colcen.wins[j+1])))
    bp.copia <- c(bp.copia,sum(copia.width[which(copia[,4]>colcen.wins[j]&copia[,4]<colcen.wins[j+1])]))
    bp.gypsy <- c(bp.gypsy,sum(gypsy.width[which(gypsy[,4]>colcen.wins[j]&gypsy[,4]<colcen.wins[j+1])]))
    bp.linel1 <- c(bp.linel1,sum(linel1.width[which(linel1[,4]>colcen.wins[j]&linel1[,4]<colcen.wins[j+1])]))
    bp.heli <- c(bp.heli,sum(heli.width[which(heli[,4]>colcen.wins[j]&heli[,4]<colcen.wins[j+1])]))
    bp.mu <- c(bp.mu,sum(mu.width[which(mu[,4]>colcen.wins[j]&mu[,4]<colcen.wins[j+1])]))
    bp.cacta <- c(bp.cacta,sum(cacta.width[which(cacta[,4]>colcen.wins[j]&cacta[,4]<colcen.wins[j+1])]))
    bp.pif <- c(bp.pif,sum(pif.width[which(pif[,4]>colcen.wins[j]&pif[,4]<colcen.wins[j+1])]))
    bp.hat <- c(bp.hat,sum(hat.width[which(hat[,4]>colcen.wins[j]&hat[,4]<colcen.wins[j+1])]))
    bp.tc1 <- c(bp.tc1,sum(tc1.width[which(tc1[,4]>colcen.wins[j]&tc1[,4]<colcen.wins[j+1])]))
    count.copia <- c(count.copia,length(which(copia[,4]>colcen.wins[j]&copia[,4]<colcen.wins[j+1])))
    count.gypsy <- c(count.gypsy,length(which(gypsy[,4]>colcen.wins[j]&gypsy[,4]<colcen.wins[j+1])))
    count.linel1 <- c(count.linel1,length(which(linel1[,4]>colcen.wins[j]&linel1[,4]<colcen.wins[j+1])))
  }
  mean.cos <- mean(wincolcen.cosall)
  plot(colcen.wins,c(bp.heli,0),type='h',col='black',ylim=c(-2500,50000),main='RNA transposons')
  par(new=T)
  plot(colcen.wins,c(bp.mu,0),type='h',col='darkorchid1',ylim=c(-2500,50000))
  par(new=T)
  plot(colcen.wins,c(bp.cacta,0),type='h',col='green',ylim=c(-2500,50000))
  par(new=T)
  plot(colcen.wins,c(bp.pif,0),type='h',col='red',ylim=c(-2500,50000))
  par(new=T)
  plot(colcen.wins,c(bp.hat,0),type='h',col='blue',ylim=c(-2500,50000))
  par(new=T)
  plot(colcen.wins,c(bp.tc1,0),type='h',col='orange',ylim=c(-2500,50000))
  rug(all.mid.nrz[i],col=1,lty=1,lwd=1,side=3,ticksize=0.05)
  rug(c(all.left.nrz[i],all.right.nrz[i]),col='blue',lwd=1,side=3,ticksize=0.05)
  rug(c(all.left.lrz[i],all.right.lrz[i]),col='purple',lwd=1,side=3,ticksize=0.05)
  rug(chr.plus[,4],col='red',ticksize=0.06)
  rug(chr.minus[,4],col='blue',ticksize=0.06)
}

##################################
# %GC content in NRZs vs outside #
##################################

chrs <- c('Chr1','Chr2','Chr3','Chr4','Chr5')
chrss <- c('chr1','chr2','chr3','chr4','chr5')
lrz.gc <- NULL
nrz.gc <- NULL
arm.gc <- NULL
for(i in 1:5){
  print(i)
  colcen.chrseq <- colcen.seq[[i]]
  lrz.seq <- colcen.chrseq[all.left.lrz[i]:all.right.lrz[i]]
  nrz.seq <- colcen.chrseq[all.left.nrz[i]:all.right.nrz[i]]
  arm.seq <- colcen.chrseq[-all.left.lrz[i]:-all.right.lrz[i]]    
  lrz.gc <- c(lrz.gc,(((length(which(lrz.seq=="g"))+length(which(lrz.seq=="c"))))/length(lrz.seq))*100)
  nrz.gc <- c(nrz.gc,(((length(which(nrz.seq=="g"))+length(which(nrz.seq=="c"))))/length(nrz.seq))*100)
  arm.gc <- c(arm.gc,(((length(which(arm.seq=="g"))+length(which(arm.seq=="c"))))/length(arm.seq))*100)
}
mean(nrz.gc)
# 38.23985
mean(lrz.gc)
# 37.63009
mean(arm.gc)
# 35.83336

wilcox.test(nrz.gc,arm.gc)
# p-value = 0.007937
wilcox.test(lrz.gc,arm.gc)
# p-value = 0.007937
wilcox.test(nrz.gc,lrz.gc)
# p-value = 0.1508

##################################################################
# plot crossovers over Col-CEN CEN178, CENH3 and DNA methylation #
##################################################################

# whole chromosome

par(mfcol=c(3,5))
par(mar=c(2,2,2,2))
chrs <- c('Chr1','Chr2','Chr3','Chr4','Chr5')
chrss <- c('chr1','chr2','chr3','chr4','chr5')
for(i in 1:5){
  print(i)
  chr.snps <- syri.snps[which(syri.snps[,1]==chrs[i]),]
  chr.inv <- syri.inv[which(syri.inv[,1]==chrs[i]),]
  chr.hdr <- syri.hdr[which(syri.hdr[,1]==chrs[i]),]
  chr.syn <- syri.syn[which(syri.syn[,1]==chrs[i]),]
  print(i)
  colcen.chrallcos <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrmeth <- colcen.meth[which(colcen.meth[,1]==chrs[i]),]
  ler.chrmeth <- ler.meth[which(ler.meth[,1]==chrs[i]),]
  colcen.chrcenh3 <- colcen.cenh3[which(colcen.cenh3[,1]==chrs[i]),]
  ler.chrcenh3 <- ler.cenh3[which(ler.cenh3[,1]==chrs[i]),]
  colcen.chrcen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  colcen.chrplus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="+"),]
  colcen.chrminus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="-"),]
  ler.chrcen178 <- lerhifi.cen178[which(lerhifi.cen178[,11]==chrs[i]),]
  ler.chrplus <- ler.chrcen178[which(ler.chrcen178[,7]=="+"),]
  ler.chrminus <- ler.chrcen178[which(ler.chrcen178[,7]=="-"),]
  colcen.chredta <- colcen.edta[which(colcen.edta[,1]==chrs[i]),]
  ler.chredta <- ler.edta[which(ler.edta[,1]==chrs[i]),]
  filt.colcen.edta <- colcen.chredta[-which(colcen.chredta[,3]=='repeat_region'),]
  colcen.edta.widths <- filt.colcen.edta[,5]-filt.colcen.edta[,4]
  filt.colcen.edta <- filt.colcen.edta[-which(colcen.edta.widths<400),]
  filt.ler.edta <- ler.chredta[-which(ler.chredta[,3]=='repeat_region'),]
  ler.edta.widths <- filt.ler.edta[,5]-filt.ler.edta[,4]
  filt.ler.edta <- filt.ler.edta[-which(ler.edta.widths<400),]
  colcen.chrathila <- colcen.athila[which(colcen.athila[,1]==chrs[i]),]
  colcen.chrgenes <- colcen.genes[which(colcen.genes[,1]==chrs[i]),]
  print(i)
  colcen.chrseq <- colcen.seq[[i]]
  lerhifi.chrseq <- lerhifi.seq[[i]]
  print(i)
  colcen.wins <- seq(1,length(colcen.chrseq),by=10000)
  wincolcen.cen178plus <- NULL
  wincolcen.cen178minus <- NULL
  wincolcen.at <- NULL
  wincolcen.gc <- NULL
  wincolcen.allcos <- NULL
  wincolcen.edta <- NULL
  wincolcen.genes <- NULL
  wincolcen.syrisnps <- NULL
  wincolcen.syrihdr <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    print(j)
    wincolcen.allcos <- c(wincolcen.allcos,length(which(colcen.chrallcos[,4]>colcen.wins[j]&colcen.chrallcos[,4]<colcen.wins[j+1])))
    wincolcen.seq <- colcen.chrseq[colcen.wins[j]:colcen.wins[j+1]]
    wincolcen.at <- c(wincolcen.at,(((length(which(wincolcen.seq=="a"))+length(which(wincolcen.seq=="t"))))/length(wincolcen.seq))*100)
    wincolcen.gc <- c(wincolcen.gc,(((length(which(wincolcen.seq=="g"))+length(which(wincolcen.seq=="c"))))/length(wincolcen.seq))*100)
    wincolcen.cen178plus <- c(wincolcen.cen178plus,length(which(colcen.chrplus[,4]>colcen.wins[j]&colcen.chrplus[,4]<colcen.wins[j+1])))
    wincolcen.cen178minus <- c(wincolcen.cen178minus,length(which(colcen.chrminus[,4]>colcen.wins[j]&colcen.chrminus[,4]<colcen.wins[j+1])))
    wincolcen.genes <- c(wincolcen.genes,length(which(colcen.chrgenes[,4]>colcen.wins[j]&colcen.chrgenes[,4]<colcen.wins[j+1])))
    wincolcen.syrisnps <- c(wincolcen.syrisnps,length(which(as.numeric(chr.snps[,2])>colcen.wins[j]&as.numeric(chr.snps[,2])<colcen.wins[j+1])))
    wincolcen.syrihdr <- c(wincolcen.syrihdr,length(which(as.numeric(chr.hdr[,2])>colcen.wins[j]&as.numeric(chr.hdr[,2])<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.allcos/3613)
  win.cMMb <- win.cM/0.01
  colcen.xlim <- c(0,length(colcen.chrseq))
  print(i)
  plot(colcen.wins,c(wincolcen.cen178plus,0),type='h',col='red',yaxt='n',xlim=colcen.xlim,main='cos + cen178 + cenh3',ylim=c(0,70))
  par(new=T)
  plot(colcen.wins,c(wincolcen.cen178minus,0),type='h',col='blue',yaxt='n',xlim=colcen.xlim,ylim=c(0,70))
  rect(chr.syn[,2],65,chr.syn[,3],70,col=4)
  rect(chr.inv[,2],60,chr.inv[,3],65,col=2)
  axis(side=4)
  par(new=T)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=colcen.xlim,yaxt='n',ylim=c(0,6))
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=colcen.xlim,ylim=c(0,44))
  rug(colcen.chrathila[,4],col='red',lwd=2)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)
  print(i)
  plot(colcen.chrmeth[,2],colcen.chrmeth[,4],type='h',col='red',ylim=c(0,1.2),xlim=colcen.xlim,main='cos + methylation')
  par(new=T)
  plot(colcen.chrmeth[,2],colcen.chrmeth[,5],type='h',col='blue',ylim=c(0,1.2),xlim=colcen.xlim)
  par(new=T)
  plot(colcen.chrmeth[,2],colcen.chrmeth[,6],type='h',col='green',ylim=c(0,1.2),xlim=colcen.xlim)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)
  rect(chr.syn[,2],1.1,chr.syn[,3],1.2,col=4)
  rect(chr.inv[,2],1,chr.inv[,3],1.1,col=2)
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=colcen.xlim,yaxt='n',ylim=c(0,44))
  axis(side=4)
  rug(colcen.chrathila[,4],col='red',lwd=2)
  abline(v=all.left.lrz[i],col=6,lwd=2)
  abline(v=all.right.lrz[i],col=6,lwd=2)
  print(i)
  plot(colcen.wins,c(wincolcen.genes,0),type='h',xlim=colcen.xlim,col='dark green',ylim=c(0,10))
  par(new=T)
  plot(colcen.wins,c(wincolcen.gc,0),type='l',xlim=colcen.xlim,col='blue',ylim=c(20,50),yaxt='n')
  rect(chr.syn[,2],48,chr.syn[,3],50,col=4)
  rect(chr.inv[,2],46,chr.inv[,3],48,col=2)
  rug(colcen.chrathila[,4],col='red',lwd=2)
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)
}

# LRZ zoom

par(mfcol=c(3,5))
par(mar=c(2,2,2,2))
chrs <- c('Chr1','Chr2','Chr3','Chr4','Chr5')
chrss <- c('chr1','chr2','chr3','chr4','chr5')
for(i in 1:5){
  print(i)
  chr.snps <- syri.snps[which(syri.snps[,1]==chrs[i]),]
  chr.inv <- syri.inv[which(syri.inv[,1]==chrs[i]),]
  chr.hdr <- syri.hdr[which(syri.hdr[,1]==chrs[i]),]
  chr.syn <- syri.syn[which(syri.syn[,1]==chrs[i]),]
  print(i)
  colcen.chrallcos <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrmeth <- colcen.meth[which(colcen.meth[,1]==chrs[i]),]
  ler.chrmeth <- ler.meth[which(ler.meth[,1]==chrs[i]),]
  colcen.chrcenh3 <- colcen.cenh3[which(colcen.cenh3[,1]==chrs[i]),]
  ler.chrcenh3 <- ler.cenh3[which(ler.cenh3[,1]==chrs[i]),]
  colcen.chrcen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  colcen.chrplus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="+"),]
  colcen.chrminus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="-"),]
  ler.chrcen180 <- lerhifi.cen178[which(lerhifi.cen178[,11]==chrs[i]),]
  ler.chrplus <- ler.chrcen178[which(ler.chrcen178[,7]=="+"),]
  ler.chrminus <- ler.chrcen178[which(ler.chrcen178[,7]=="-"),]
  colcen.chredta <- colcen.edta[which(colcen.edta[,1]==chrs[i]),]
  ler.chredta <- ler.edta[which(ler.edta[,1]==chrs[i]),]
  filt.colcen.edta <- colcen.chredta[-which(colcen.chredta[,3]=='repeat_region'),]
  colcen.edta.widths <- filt.colcen.edta[,5]-filt.colcen.edta[,4]
  filt.colcen.edta <- filt.colcen.edta[-which(colcen.edta.widths<400),]
  filt.ler.edta <- ler.chredta[-which(ler.chredta[,3]=='repeat_region'),]
  ler.edta.widths <- filt.ler.edta[,5]-filt.ler.edta[,4]
  filt.ler.edta <- filt.ler.edta[-which(ler.edta.widths<400),]
  colcen.chrathila <- colcen.athila[which(colcen.athila[,1]==chrs[i]),]
  colcen.chrgenes <- colcen.genes[which(colcen.genes[,1]==chrs[i]),]
  print(i)
  colcen.chrseq <- colcen.seq[[i]]
  lerhifi.chrseq <- lerhifi.seq[[i]]
  print(i)
  colcen.wins <- seq(1,length(colcen.chrseq),by=10000)
  wincolcen.cen178plus <- NULL
  wincolcen.cen178minus <- NULL
  wincolcen.at <- NULL
  wincolcen.gc <- NULL
  wincolcen.allcos <- NULL
  wincolcen.edta <- NULL
  wincolcen.genes <- NULL
  wincolcen.syrisnps <- NULL
  wincolcen.syrihdr <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    print(j)
    wincolcen.allcos <- c(wincolcen.allcos,length(which(colcen.chrallcos[,4]>colcen.wins[j]&colcen.chrallcos[,4]<colcen.wins[j+1])))
    wincolcen.seq <- colcen.chrseq[colcen.wins[j]:colcen.wins[j+1]]
    wincolcen.at <- c(wincolcen.at,(((length(which(wincolcen.seq=="a"))+length(which(wincolcen.seq=="t"))))/length(wincolcen.seq))*100)
    wincolcen.gc <- c(wincolcen.gc,(((length(which(wincolcen.seq=="g"))+length(which(wincolcen.seq=="c"))))/length(wincolcen.seq))*100)
    wincolcen.cen178plus <- c(wincolcen.cen178plus,length(which(colcen.chrplus[,4]>colcen.wins[j]&colcen.chrplus[,4]<colcen.wins[j+1])))
    wincolcen.cen178minus <- c(wincolcen.cen178minus,length(which(colcen.chrminus[,4]>colcen.wins[j]&colcen.chrminus[,4]<colcen.wins[j+1])))
    wincolcen.genes <- c(wincolcen.genes,length(which(colcen.chrgenes[,4]>colcen.wins[j]&colcen.chrgenes[,4]<colcen.wins[j+1])))
    wincolcen.syrisnps <- c(wincolcen.syrisnps,length(which(as.numeric(chr.snps[,2])>colcen.wins[j]&as.numeric(chr.snps[,2])<colcen.wins[j+1])))
    wincolcen.syrihdr <- c(wincolcen.syrihdr,length(which(as.numeric(chr.hdr[,2])>colcen.wins[j]&as.numeric(chr.hdr[,2])<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.allcos/3613)
  win.cMMb <- win.cM/0.01
  cenlim.start <- all.left.lrz[i]-10000
  cenlim.end <- all.right.lrz[i]+10000
  cenlim <- c(cenlim.start,cenlim.end)
  colcen.xlim <- c(0,length(colcen.chrseq))
  print(i)
  plot(colcen.wins,c(wincolcen.cen178plus,0),type='h',col='red',yaxt='n',xlim=cenlim,main='cos + cen178 + cenh3',ylim=c(0,70))
  par(new=T)
  plot(colcen.wins,c(wincolcen.cen178minus,0),type='h',col='blue',yaxt='n',xlim=cenlim,ylim=c(0,70))
  rect(chr.syn[,2],65,chr.syn[,3],70,col=4)
  rect(chr.inv[,2],60,chr.inv[,3],65,col=2)
  axis(side=4)
  par(new=T)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=cenlim,yaxt='n',ylim=c(0,6))
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=cenlim,ylim=c(0,44))
  rug(colcen.chrathila[,4],col='red',lwd=2)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)
  print(i)
  plot(colcen.chrmeth[,2],colcen.chrmeth[,4],type='h',col='red',ylim=c(0,1.2),xlim=cenlim,main='cos + methylation')
  par(new=T)
  plot(colcen.chrmeth[,2],colcen.chrmeth[,5],type='h',col='blue',ylim=c(0,1.2),xlim=cenlim)
  par(new=T)
  plot(colcen.chrmeth[,2],colcen.chrmeth[,6],type='h',col='green',ylim=c(0,1.2),xlim=cenlim)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)
  rect(chr.syn[,2],1.1,chr.syn[,3],1.2,col=4)
  rect(chr.inv[,2],1,chr.inv[,3],1.1,col=2)
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=cenlim,yaxt='n',ylim=c(0,44))
  axis(side=4)
  rug(colcen.chrathila[,4],col='red',lwd=2)
  print(i)
  plot(colcen.wins,c(wincolcen.genes,0),type='h',xlim=cenlim,col='dark green',ylim=c(0,10))
  par(new=T)
  plot(colcen.wins,c(wincolcen.gc,0),type='l',xlim=cenlim,col='blue',ylim=c(20,50),yaxt='n')
  rect(chr.syn[,2],48,chr.syn[,3],50,col=4)
  rect(chr.inv[,2],46,chr.inv[,3],48,col=2)
  rug(colcen.chrathila[,4],col='red',lwd=2)
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)
}

# plotting to obtain axes for CENH3
par(mfcol=c(3,1))
par(mar=c(2,2,2,2))
print(i)
colcen.chrcenh3 <- colcen.cenh3[which(colcen.cenh3[,1]==chrs[i]),]
colcen.xlim <- c(0,length(colcen.chrseq))
plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=colcen.xlim,ylim=c(0,6))

##################################################################################
# plot H3K9me2, H3K27me1, H2A.W6 and H2A.W7 across the centromeres vs crossovers #
##################################################################################

chrom.dat <- read.table(file='ChIP_input_log2ChIPinput_RaGOO_v2.0_genome_norm_coverage_matrix_10kb_unsmoothed.tsv',header=T)

# whole chromosome crossovers vs heterochromatin

par(mfcol=c(4,5))
par(mar=c(2,2,2,2))
chrs <- c('Chr1','Chr2','Chr3','Chr4','Chr5')
chrss <- c('chr1','chr2','chr3','chr4','chr5')
for(i in 1:5){
  print(i)
  chrchrom.dat <- chrom.dat[which(chrom.dat[,1]==chrs[i]),]
  chr.snps <- syri.snps[which(syri.snps[,1]==chrs[i]),]
  chr.inv <- syri.inv[which(syri.inv[,1]==chrs[i]),]
  chr.hdr <- syri.hdr[which(syri.hdr[,1]==chrs[i]),]
  chr.syn <- syri.syn[which(syri.syn[,1]==chrs[i]),]
  print(i)
  colcen.chrallcos <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrcenh3 <- colcen.cenh3[which(colcen.cenh3[,1]==chrs[i]),]
  ler.chrcenh3 <- ler.cenh3[which(ler.cenh3[,1]==chrs[i]),]
  colcen.chrcen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  colcen.chrplus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="+"),]
  colcen.chrminus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="-"),]
  ler.chrcen178 <- lerhifi.cen178[which(lerhifi.cen178[,11]==chrs[i]),]
  ler.chrplus <- ler.chrcen178[which(ler.chrcen178[,7]=="+"),]
  ler.chrminus <- ler.chrcen178[which(ler.chrcen178[,7]=="-"),]
  colcen.chredta <- colcen.edta[which(colcen.edta[,1]==chrs[i]),]
  ler.chredta <- ler.edta[which(ler.edta[,1]==chrs[i]),]
  filt.colcen.edta <- colcen.chredta[-which(colcen.chredta[,3]=='repeat_region'),]
  colcen.edta.widths <- filt.colcen.edta[,5]-filt.colcen.edta[,4]
  filt.colcen.edta <- filt.colcen.edta[-which(colcen.edta.widths<400),]
  filt.ler.edta <- ler.chredta[-which(ler.chredta[,3]=='repeat_region'),]
  ler.edta.widths <- filt.ler.edta[,5]-filt.ler.edta[,4]
  filt.ler.edta <- filt.ler.edta[-which(ler.edta.widths<400),]
  colcen.chrathila <- colcen.athila[which(colcen.athila[,1]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  lerhifi.chrseq <- lerhifi.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=10000)
  wincolcen.allcos <- NULL
  wincolcen.syrisnps <- NULL
  wincolcen.syrihdr <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    print(j)
    wincolcen.allcos <- c(wincolcen.allcos,length(which(colcen.chrallcos[,4]>colcen.wins[j]&colcen.chrallcos[,4]<colcen.wins[j+1])))
    wincolcen.syrisnps <- c(wincolcen.syrisnps,length(which(as.numeric(chr.snps[,2])>colcen.wins[j]&as.numeric(chr.snps[,2])<colcen.wins[j+1])))
    wincolcen.syrihdr <- c(wincolcen.syrihdr,length(which(as.numeric(chr.hdr[,2])>colcen.wins[j]&as.numeric(chr.hdr[,2])<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.allcos/3613)
  win.cMMb <- win.cM/0.01
  colcen.xlim <- c(0,length(colcen.chrseq))
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=colcen.xlim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_WT_H3K9me2_Rep1_ChIP_WT_REC8_Myc_Rep1_input,type='h',col='red',xlim=colcen.xlim,ylim=c(0,3.4),main='H3K9me2')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=colcen.xlim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],44,chr.syn[,3],45,col=4)
  #rect(chr.inv[,2],44,chr.inv[,3],45,col=2)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=colcen.xlim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_H3K27me1_ChIP_set5_input_set5,type='h',col='red',xlim=colcen.xlim,ylim=c(0,1.2),main='H3K27me1')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=colcen.xlim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=colcen.xlim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_H2AW6_ChIP_set1_input_set1,type='h',col='red',xlim=colcen.xlim,ylim=c(0,2.4),main='H2AW6')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=colcen.xlim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=colcen.xlim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_H2AW7_ChIP_set2_input_set2,type='h',col='red',xlim=colcen.xlim,ylim=c(0,2.2),main='H2AW7')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=colcen.xlim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
}

# LRZ zoom crossovers vs heterochromatin

par(mfcol=c(4,5))
par(mar=c(2,2,2,2))
chrs <- c('Chr1','Chr2','Chr3','Chr4','Chr5')
chrss <- c('chr1','chr2','chr3','chr4','chr5')
for(i in 1:5){
  print(i)
  chrchrom.dat <- chrom.dat[which(chrom.dat[,1]==chrs[i]),]
  chr.snps <- syri.snps[which(syri.snps[,1]==chrs[i]),]
  chr.inv <- syri.inv[which(syri.inv[,1]==chrs[i]),]
  chr.hdr <- syri.hdr[which(syri.hdr[,1]==chrs[i]),]
  chr.syn <- syri.syn[which(syri.syn[,1]==chrs[i]),]
  print(i)
  colcen.chrallcos <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrcenh3 <- colcen.cenh3[which(colcen.cenh3[,1]==chrs[i]),]
  ler.chrcenh3 <- ler.cenh3[which(ler.cenh3[,1]==chrs[i]),]
  colcen.chrcen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  colcen.chrplus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="+"),]
  colcen.chrminus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="-"),]
  ler.chrcen178 <- lerhifi.cen178[which(lerhifi.cen178[,11]==chrs[i]),]
  ler.chrplus <- ler.chrcen178[which(ler.chrcen178[,7]=="+"),]
  ler.chrminus <- ler.chrcen178[which(ler.chrcen178[,7]=="-"),]
  colcen.chredta <- colcen.edta[which(colcen.edta[,1]==chrs[i]),]
  ler.chredta <- ler.edta[which(ler.edta[,1]==chrs[i]),]
  filt.colcen.edta <- colcen.chredta[-which(colcen.chredta[,3]=='repeat_region'),]
  colcen.edta.widths <- filt.colcen.edta[,5]-filt.colcen.edta[,4]
  filt.colcen.edta <- filt.colcen.edta[-which(colcen.edta.widths<400),]
  filt.ler.edta <- ler.chredta[-which(ler.chredta[,3]=='repeat_region'),]
  ler.edta.widths <- filt.ler.edta[,5]-filt.ler.edta[,4]
  filt.ler.edta <- filt.ler.edta[-which(ler.edta.widths<400),]
  colcen.chrathila <- colcen.athila[which(colcen.athila[,1]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  lerhifi.chrseq <- lerhifi.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=10000)
  wincolcen.allcos <- NULL
  wincolcen.syrisnps <- NULL
  wincolcen.syrihdr <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    print(j)
    wincolcen.allcos <- c(wincolcen.allcos,length(which(colcen.chrallcos[,4]>colcen.wins[j]&colcen.chrallcos[,4]<colcen.wins[j+1])))
    wincolcen.syrisnps <- c(wincolcen.syrisnps,length(which(as.numeric(chr.snps[,2])>colcen.wins[j]&as.numeric(chr.snps[,2])<colcen.wins[j+1])))
    wincolcen.syrihdr <- c(wincolcen.syrihdr,length(which(as.numeric(chr.hdr[,2])>colcen.wins[j]&as.numeric(chr.hdr[,2])<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.allcos/3613)
  win.cMMb <- win.cM/0.01
  cenlim.start <- all.left.lrz[i]-10000
  cenlim.end <- all.right.lrz[i]+10000
  cen.lim <- c(cenlim.start,cenlim.end)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=cen.lim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_WT_H3K9me2_Rep1_ChIP_WT_REC8_Myc_Rep1_input,type='h',col='red',xlim=cen.lim,ylim=c(0,3.4),main='H3K9me2')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=cen.lim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=cen.lim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_H3K27me1_ChIP_set5_input_set5,type='h',col='red',xlim=cen.lim,ylim=c(0,1.2),main='H3K27me1')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=cen.lim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=cen.lim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_H2AW6_ChIP_set1_input_set1,type='h',col='red',xlim=cen.lim,ylim=c(0,2.4),main='H2AW6')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=cen.lim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=cen.lim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_H2AW7_ChIP_set2_input_set2,type='h',col='red',xlim=cen.lim,ylim=c(0,2.2),main='H2AW7')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=cen.lim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
}

###########################################################################
# plot REC8, ASY1 and SPO11-1-oligos across the centromeres vs crossovers #
###########################################################################

# whole chromosome crossovers vs REC8 ASY1 and SPO11-1

par(mfcol=c(3,5))
par(mar=c(2,2,2,2))
chrs <- c('Chr1','Chr2','Chr3','Chr4','Chr5')
chrss <- c('chr1','chr2','chr3','chr4','chr5')
for(i in 1:5){
  print(i)
  chrchrom.dat <- chrom.dat[which(chrom.dat[,1]==chrs[i]),]
  chr.snps <- syri.snps[which(syri.snps[,1]==chrs[i]),]
  chr.inv <- syri.inv[which(syri.inv[,1]==chrs[i]),]        
  chr.hdr <- syri.hdr[which(syri.hdr[,1]==chrs[i]),]
  chr.syn <- syri.syn[which(syri.syn[,1]==chrs[i]),]
  print(i)
  colcen.chrallcos <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrcenh3 <- colcen.cenh3[which(colcen.cenh3[,1]==chrs[i]),]
  ler.chrcenh3 <- ler.cenh3[which(ler.cenh3[,1]==chrs[i]),]
  colcen.chrcen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  colcen.chrplus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="+"),]
  colcen.chrminus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="-"),]
  ler.chrcen178 <- lerhifi.cen178[which(lerhifi.cen178[,11]==chrs[i]),]
  ler.chrplus <- ler.chrcen178[which(ler.chrcen178[,7]=="+"),]
  ler.chrminus <- ler.chrcen178[which(ler.chrcen178[,7]=="-"),]
  colcen.chredta <- colcen.edta[which(colcen.edta[,1]==chrs[i]),]
  ler.chredta <- ler.edta[which(ler.edta[,1]==chrs[i]),]
  filt.colcen.edta <- colcen.chredta[-which(colcen.chredta[,3]=='repeat_region'),]
  colcen.edta.widths <- filt.colcen.edta[,5]-filt.colcen.edta[,4]
  filt.colcen.edta <- filt.colcen.edta[-which(colcen.edta.widths<400),]
  filt.ler.edta <- ler.chredta[-which(ler.chredta[,3]=='repeat_region'),]
  ler.edta.widths <- filt.ler.edta[,5]-filt.ler.edta[,4]
  filt.ler.edta <- filt.ler.edta[-which(ler.edta.widths<400),]
  colcen.chrathila <- colcen.athila[which(colcen.athila[,1]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  lerhifi.chrseq <- lerhifi.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=10000)
  wincolcen.allcos <- NULL
  wincolcen.syrisnps <- NULL
  wincolcen.syrihdr <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    print(j)
    wincolcen.allcos <- c(wincolcen.allcos,length(which(colcen.chrallcos[,4]>colcen.wins[j]&colcen.chrallcos[,4]<colcen.wins[j+1])))
    wincolcen.syrisnps <- c(wincolcen.syrisnps,length(which(as.numeric(chr.snps[,2])>colcen.wins[j]&as.numeric(chr.snps[,2])<colcen.wins[j+1])))
    wincolcen.syrihdr <- c(wincolcen.syrihdr,length(which(as.numeric(chr.hdr[,2])>colcen.wins[j]&as.numeric(chr.hdr[,2])<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.allcos/3613)
  win.cMMb <- win.cM/0.01
  colcen.xlim <- c(0,length(colcen.chrseq))
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=colcen.xlim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_WT_REC8_HA_Rep2_ChIP_WT_REC8_Myc_Rep1_input,type='h',col='red',xlim=colcen.xlim,main='REC8',ylim=c(0,1))
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=colcen.xlim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=colcen.xlim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_WT_ASY1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,type='h',ylim=c(0,0.95),col='red',xlim=colcen.xlim,main='ASY1')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=colcen.xlim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=colcen.xlim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_WT_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,type='h',col='red',ylim=c(0,0.85),xlim=colcen.xlim,main='SPO11-1')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=colcen.xlim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
  print(i)
}

# LRZ zoom crossovers vs REC8 ASY1 and SPO11-1

par(mfcol=c(3,5))
par(mar=c(2,2,2,2))
chrs <- c('Chr1','Chr2','Chr3','Chr4','Chr5')
chrss <- c('chr1','chr2','chr3','chr4','chr5')
for(i in 1:5){
  print(i)
  chrchrom.dat <- chrom.dat[which(chrom.dat[,1]==chrs[i]),]
  chr.snps <- syri.snps[which(syri.snps[,1]==chrs[i]),]
  chr.inv <- syri.inv[which(syri.inv[,1]==chrs[i]),]        
  chr.hdr <- syri.hdr[which(syri.hdr[,1]==chrs[i]),]
  chr.syn <- syri.syn[which(syri.syn[,1]==chrs[i]),]
  print(i)
  colcen.chrallcos <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrcenh3 <- colcen.cenh3[which(colcen.cenh3[,1]==chrs[i]),]
  ler.chrcenh3 <- ler.cenh3[which(ler.cenh3[,1]==chrs[i]),]
  colcen.chrcen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  colcen.chrplus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="+"),]
  colcen.chrminus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="-"),]
  ler.chrcen178 <- lerhifi.cen178[which(lerhifi.cen178[,11]==chrs[i]),]
  ler.chrplus <- ler.chrcen178[which(ler.chrcen178[,7]=="+"),]
  ler.chrminus <- ler.chrcen178[which(ler.chrcen178[,7]=="-"),]
  colcen.chredta <- colcen.edta[which(colcen.edta[,1]==chrs[i]),]
  ler.chredta <- ler.edta[which(ler.edta[,1]==chrs[i]),]
  filt.colcen.edta <- colcen.chredta[-which(colcen.chredta[,3]=='repeat_region'),]
  colcen.edta.widths <- filt.colcen.edta[,5]-filt.colcen.edta[,4]
  filt.colcen.edta <- filt.colcen.edta[-which(colcen.edta.widths<400),]
  filt.ler.edta <- ler.chredta[-which(ler.chredta[,3]=='repeat_region'),]
  ler.edta.widths <- filt.ler.edta[,5]-filt.ler.edta[,4]
  filt.ler.edta <- filt.ler.edta[-which(ler.edta.widths<400),]
  colcen.chrathila <- colcen.athila[which(colcen.athila[,1]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  lerhifi.chrseq <- lerhifi.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=10000)
  wincolcen.allcos <- NULL
  wincolcen.syrisnps <- NULL
  wincolcen.syrihdr <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    print(j)
    wincolcen.allcos <- c(wincolcen.allcos,length(which(colcen.chrallcos[,4]>colcen.wins[j]&colcen.chrallcos[,4]<colcen.wins[j+1])))
    wincolcen.syrisnps <- c(wincolcen.syrisnps,length(which(as.numeric(chr.snps[,2])>colcen.wins[j]&as.numeric(chr.snps[,2])<colcen.wins[j+1])))
    wincolcen.syrihdr <- c(wincolcen.syrihdr,length(which(as.numeric(chr.hdr[,2])>colcen.wins[j]&as.numeric(chr.hdr[,2])<colcen.wins[j+1])))
  }
  win.cM <- 100*(wincolcen.allcos/3613)
  win.cMMb <- win.cM/0.01
  cenlim.start <- all.left.lrz[i]-10000
  cenlim.end <- all.right.lrz[i]+10000
  cen.lim <- c(cenlim.start,cenlim.end)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=cen.lim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_WT_REC8_HA_Rep2_ChIP_WT_REC8_Myc_Rep1_input,type='h',col='red',xlim=cen.lim,main='REC8',ylim=c(0,1))
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=cen.lim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=cen.lim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_WT_ASY1_Rep1_ChIP_WT_REC8_Myc_Rep1_input,type='h',ylim=c(0,0.95),col='red',xlim=cen.lim,main='ASY1')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=cen.lim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
  print(i)
  plot(colcen.chrcenh3[,2],colcen.chrcenh3[,6],type='h',col='green',xlim=cen.lim,yaxt='n',ylim=c(0,5),yaxt='n')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_WT_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,type='h',col='red',ylim=c(0,0.85),xlim=cen.lim,main='SPO11-1')
  par(new=T)
  plot(colcen.wins,c(win.cMMb,0),type='h',xlim=cen.lim,ylim=c(0,44),yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)    
  #rect(chr.syn[,2],19,chr.syn[,3],20,col=4)
  #rect(chr.inv[,2],19,chr.inv[,3],20,col=2)
  print(i)
}

###################################################
# quantify various factors in arms, LRZs and NRZs #
###################################################

chrom.dat <- read.table(file='ChIP_input_log2ChIPinput_RaGOO_v2.0_genome_norm_coverage_matrix_10kb_unsmoothed.tsv',header=T)
all.arm.spo11 <- NULL
all.lrz.spo11 <- NULL
all.nrz.spo11 <- NULL
all.arm.asy1 <- NULL
all.lrz.asy1 <- NULL
all.nrz.asy1 <- NULL
all.arm.rec8 <- NULL
all.lrz.rec8 <- NULL
all.nrz.rec8 <- NULL
all.arm.h3k9 <- NULL
all.lrz.h3k9 <- NULL
all.nrz.h3k9 <- NULL
all.arm.h3k27 <- NULL
all.lrz.h3k27 <- NULL
all.nrz.h3k27 <- NULL
all.arm.h2aw6 <- NULL
all.lrz.h2aw6 <- NULL
all.nrz.h2aw6 <- NULL
all.arm.h2aw7 <- NULL
all.lrz.h2aw7 <- NULL
all.nrz.h2aw7 <- NULL
all.arm.cg <- NULL
all.lrz.cg <- NULL
all.nrz.cg <- NULL
all.arm.chg <- NULL
all.lrz.chg <- NULL
all.nrz.chg <- NULL
all.arm.chh <- NULL
all.lrz.chh <- NULL
all.nrz.chh <- NULL
all.arm.cenh3 <- NULL
all.lrz.cenh3 <- NULL
all.nrz.cenh3 <- NULL
for(i in 1:5){
  print(i)
  colcen.chrcenh3 <- colcen.cenh3[which(colcen.cenh3[,1]==chrs[i]),]
  colcen.chrmeth <- colcen.meth[which(colcen.meth[,1]==chrs[i]),]
  colcen.chrmeth <- colcen.chrmeth[-length(colcen.chrmeth[,1]),]
  chrchrom.dat <- chrom.dat[which(chrom.dat[,1]==chrs[i]),]
  chr.coords <- chrchrom.dat[,2]
  chr.spo11 <- chrchrom.dat$log2_WT_SPO11oligos_Rep1_WT_gDNA_Rep1_R1
  chr.rec8 <- chrchrom.dat$log2_WT_REC8_HA_Rep2_ChIP_WT_REC8_Myc_Rep1_input
  chr.asy1 <- chrchrom.dat$log2_WT_ASY1_Rep1_ChIP_WT_REC8_Myc_Rep1_input
  chr.h3k9 <- chrchrom.dat$log2_WT_H3K9me2_Rep1_ChIP_WT_REC8_Myc_Rep1_input
  chr.h2aw6 <- chrchrom.dat$log2_H2AW6_ChIP_set1_input_set1
  chr.h2aw7 <- chrchrom.dat$log2_H2AW7_ChIP_set2_input_set2
  chr.h3k27 <- chrchrom.dat$log2_H3K27me1_ChIP_set5_input_set5
  all.arm.cenh3 <- c(all.arm.cenh3,mean(colcen.chrcenh3[c(which(colcen.chrcenh3[,2]<all.left.lrz[i]),which(colcen.chrcenh3[,2]>all.right.lrz[i])),6]))
  all.lrz.cenh3 <- c(all.lrz.cenh3,mean(colcen.chrcenh3[c(which(colcen.chrcenh3[,2]>all.left.lrz[i]&colcen.chrcenh3[,2]<all.left.nrz[i]),which(colcen.chrcenh3[,2]<all.right.lrz[i]&colcen.chrcenh3[,2]>all.right.nrz[i])),6]))
  all.nrz.cenh3 <- c(all.nrz.cenh3,mean(colcen.chrcenh3[which(colcen.chrcenh3[,2]>all.left.nrz[i]&colcen.chrcenh3[,2]<all.right.nrz[i]),6]))
  all.arm.h2aw7 <- c(all.arm.h2aw7,mean(chr.h2aw7[c(which(chr.coords<all.left.lrz[i]),which(chr.coords>all.right.lrz[i]))]))
  all.lrz.h2aw7 <- c(all.lrz.h2aw7,mean(chr.h2aw7[c(which(chr.coords>all.left.lrz[i]&chr.coords<all.left.nrz[i]),which(chr.coords<all.right.lrz[i]&chr.coords>all.right.nrz[i]))]))
  all.nrz.h2aw7 <- c(all.nrz.h2aw7,mean(chr.h2aw7[which(chr.coords>all.left.nrz[i]&chr.coords<all.right.nrz[i])]))
  all.arm.h2aw6 <- c(all.arm.h2aw6,mean(chr.h2aw6[c(which(chr.coords<all.left.lrz[i]),which(chr.coords>all.right.lrz[i]))]))
  all.lrz.h2aw6 <- c(all.lrz.h2aw6,mean(chr.h2aw6[c(which(chr.coords>all.left.lrz[i]&chr.coords<all.left.nrz[i]),which(chr.coords<all.right.lrz[i]&chr.coords>all.right.nrz[i]))]))
  all.nrz.h2aw6 <- c(all.nrz.h2aw6,mean(chr.h2aw6[which(chr.coords>all.left.nrz[i]&chr.coords<all.right.nrz[i])]))
  all.arm.h3k9 <- c(all.arm.h3k9,mean(chr.h3k9[c(which(chr.coords<all.left.lrz[i]),which(chr.coords>all.right.lrz[i]))]))
  all.lrz.h3k9 <- c(all.lrz.h3k9,mean(chr.h3k9[c(which(chr.coords>all.left.lrz[i]&chr.coords<all.left.nrz[i]),which(chr.coords<all.right.lrz[i]&chr.coords>all.right.nrz[i]))]))
  all.nrz.h3k9 <- c(all.nrz.h3k9,mean(chr.h3k9[which(chr.coords>all.left.nrz[i]&chr.coords<all.right.nrz[i])]))
  all.arm.h3k27 <- c(all.arm.h3k27,mean(chr.h3k27[c(which(chr.coords<all.left.lrz[i]),which(chr.coords>all.right.lrz[i]))]))
  all.lrz.h3k27 <- c(all.lrz.h3k27,mean(chr.h3k27[c(which(chr.coords>all.left.lrz[i]&chr.coords<all.left.nrz[i]),which(chr.coords<all.right.lrz[i]&chr.coords>all.right.nrz[i]))]))
  all.nrz.h3k27 <- c(all.nrz.h3k27,mean(chr.h3k27[which(chr.coords>all.left.nrz[i]&chr.coords<all.right.nrz[i])]))
  all.arm.asy1 <- c(all.arm.asy1,mean(chr.asy1[c(which(chr.coords<all.left.lrz[i]),which(chr.coords>all.right.lrz[i]))]))
  all.lrz.asy1 <- c(all.lrz.asy1,mean(chr.asy1[c(which(chr.coords>all.left.lrz[i]&chr.coords<all.left.nrz[i]),which(chr.coords<all.right.lrz[i]&chr.coords>all.right.nrz[i]))]))
  all.nrz.asy1 <- c(all.nrz.asy1,mean(chr.asy1[which(chr.coords>all.left.nrz[i]&chr.coords<all.right.nrz[i])]))
  all.arm.rec8 <- c(all.arm.rec8,mean(chr.rec8[c(which(chr.coords<all.left.lrz[i]),which(chr.coords>all.right.lrz[i]))]))
  all.lrz.rec8 <- c(all.lrz.rec8,mean(chr.rec8[c(which(chr.coords>all.left.lrz[i]&chr.coords<all.left.nrz[i]),which(chr.coords<all.right.lrz[i]&chr.coords>all.right.nrz[i]))]))
  all.nrz.rec8 <- c(all.nrz.rec8,mean(chr.rec8[which(chr.coords>all.left.nrz[i]&chr.coords<all.right.nrz[i])]))
  all.arm.spo11 <- c(all.arm.spo11,mean(chr.spo11[c(which(chr.coords<all.left.lrz[i]),which(chr.coords>all.right.lrz[i]))]))
  all.lrz.spo11 <- c(all.lrz.spo11,mean(chr.spo11[c(which(chr.coords>all.left.lrz[i]&chr.coords<all.left.nrz[i]),which(chr.coords<all.right.lrz[i]&chr.coords>all.right.nrz[i]))]))
  all.nrz.spo11 <- c(all.nrz.spo11,mean(chr.spo11[which(chr.coords>all.left.nrz[i]&chr.coords<all.right.nrz[i])]))
  all.arm.cg <- c(all.arm.cg,mean(colcen.chrmeth[c(which(colcen.chrmeth[,2]<all.left.lrz[i]),which(colcen.chrmeth[,2]>all.right.lrz[i])),4]))
  all.lrz.cg <- c(all.lrz.cg,mean(colcen.chrmeth[c(which(colcen.chrmeth[,2]>all.left.lrz[i]&colcen.chrmeth[,2]<all.left.nrz[i]),which(colcen.chrmeth[,2]<all.right.lrz[i]&colcen.chrmeth[,2]>all.right.nrz[i])),4]))
  all.nrz.cg <- c(all.nrz.cg,mean(colcen.chrmeth[which(colcen.chrmeth[,2]>all.left.nrz[i]&colcen.chrmeth[,2]<all.right.nrz[i]),4]))
  all.arm.chg <- c(all.arm.chg,mean(colcen.chrmeth[c(which(colcen.chrmeth[,2]<all.left.lrz[i]),which(colcen.chrmeth[,2]>all.right.lrz[i])),5]))
  all.lrz.chg <- c(all.lrz.chg,mean(colcen.chrmeth[c(which(colcen.chrmeth[,2]>all.left.lrz[i]&colcen.chrmeth[,2]<all.left.nrz[i]),which(colcen.chrmeth[,2]<all.right.lrz[i]&colcen.chrmeth[,2]>all.right.nrz[i])),5]))
  all.nrz.chg <- c(all.nrz.chg,mean(colcen.chrmeth[which(colcen.chrmeth[,2]>all.left.nrz[i]&colcen.chrmeth[,2]<all.right.nrz[i]),5]))
  all.arm.chh <- c(all.arm.chh,mean(colcen.chrmeth[c(which(colcen.chrmeth[,2]<all.left.lrz[i]),which(colcen.chrmeth[,2]>all.right.lrz[i])),6]))
  all.lrz.chh <- c(all.lrz.chh,mean(colcen.chrmeth[c(which(colcen.chrmeth[,2]>all.left.lrz[i]&colcen.chrmeth[,2]<all.left.nrz[i]),which(colcen.chrmeth[,2]<all.right.lrz[i]&colcen.chrmeth[,2]>all.right.nrz[i])),6]))
  all.nrz.chh <- c(all.nrz.chh,mean(colcen.chrmeth[which(colcen.chrmeth[,2]>all.left.nrz[i]&colcen.chrmeth[,2]<all.right.nrz[i]),6]))
}
mean(all.arm.cenh3) -0.6585342
mean(all.lrz.cenh3) -0.1370754
mean(all.nrz.cenh3) 1.846096
mean(all.arm.h3k9) -0.6213336
mean(all.lrz.h3k9) 1.263242
mean(all.nrz.h3k9) 1.562817
mean(all.arm.h3k27) -0.05817143
mean(all.lrz.h3k27) 0.5513931
mean(all.nrz.h3k27) 0.3462527
mean(all.arm.h2aw6) -0.4869081
mean(all.lrz.h2aw6) 1.311294
mean(all.nrz.h2aw6) 1.037431
mean(all.arm.h2aw7) -0.2450053
mean(all.lrz.h2aw7) 0.9868705
mean(all.nrz.h2aw7) 0.283572
mean(all.arm.rec8) -0.06766098
mean(all.lrz.rec8) 0.3676711
mean(all.nrz.rec8) 0.4193986
mean(all.arm.asy1) -0.07592318
mean(all.lrz.asy1) 0.3165892
mean(all.nrz.asy1) 0.371229
mean(all.arm.spo11) 0.01023579
mean(all.lrz.spo11) -0.226827
mean(all.nrz.spo11) -0.3105742
mean(all.arm.cg) 0.2176132
mean(all.lrz.cg) 0.757213
mean(all.nrz.cg) 0.7950914
mean(all.arm.chg) 0.06233326
mean(all.lrz.chg) 0.3819987
mean(all.nrz.chg) 0.3052455
mean(all.arm.chh) 0.02067925
mean(all.lrz.chh) 0.07577885
mean(all.nrz.chh) 0.0691473

##########################################
# plot genes, SPO11-1-oligos, crossovers #
##########################################

chrom.dat <- read.table(file='ChIP_input_log2ChIPinput_RaGOO_v2.0_genome_norm_coverage_matrix_10kb_unsmoothed.tsv',header=T)

par(mfcol=c(5,1))
par(mar=c(2,2,2,2))
chrs <- c('Chr1','Chr2','Chr3','Chr4','Chr5')
chrss <- c('chr1','chr2','chr3','chr4','chr5')
dat.lrz.genes <- NULL
dat.lrz.spo11 <- NULL
dat.lrz.cos <- NULL
for(i in 1:5){
  print(i)
  chrchrom.dat <- chrom.dat[which(chrom.dat[,1]==chrs[i]),]
  chr.snps <- syri.snps[which(syri.snps[,1]==chrs[i]),]
  chr.inv <- syri.inv[which(syri.inv[,1]==chrs[i]),]        
  chr.hdr <- syri.hdr[which(syri.hdr[,1]==chrs[i]),]
  chr.syn <- syri.syn[which(syri.syn[,1]==chrs[i]),]
  print(i)
  colcen.chrallcos <- cos.all[which(cos.all[,3]==chrs[i]),]
  colcen.chrcenh3 <- colcen.cenh3[which(colcen.cenh3[,1]==chrs[i]),]
  ler.chrcenh3 <- ler.cenh3[which(ler.cenh3[,1]==chrs[i]),]
  colcen.chrcen178 <- colcen.cen178[which(colcen.cen178[,1]==chrs[i]),]
  colcen.chrplus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="+"),]
  colcen.chrminus <- colcen.chrcen178[which(colcen.chrcen178[,7]=="-"),]
  ler.chrcen178 <- lerhifi.cen178[which(lerhifi.cen178[,11]==chrs[i]),]
  ler.chrplus <- ler.chrcen178[which(ler.chrcen178[,7]=="+"),]
  ler.chrminus <- ler.chrcen178[which(ler.chrcen178[,7]=="-"),]
  colcen.chredta <- colcen.edta[which(colcen.edta[,1]==chrs[i]),]
  ler.chredta <- ler.edta[which(ler.edta[,1]==chrs[i]),]
  filt.colcen.edta <- colcen.chredta[-which(colcen.chredta[,3]=='repeat_region'),]
  colcen.edta.widths <- filt.colcen.edta[,5]-filt.colcen.edta[,4]
  filt.colcen.edta <- filt.colcen.edta[-which(colcen.edta.widths<400),]
  filt.ler.edta <- ler.chredta[-which(ler.chredta[,3]=='repeat_region'),]
  ler.edta.widths <- filt.ler.edta[,5]-filt.ler.edta[,4]
  filt.ler.edta <- filt.ler.edta[-which(ler.edta.widths<400),]
  colcen.chrgenes <- colcen.genes[which(colcen.genes[,1]==chrs[i]),]
  colcen.chrathila <- colcen.athila[which(colcen.athila[,1]==chrs[i]),]
  colcen.chrseq <- colcen.seq[[i]]
  lerhifi.chrseq <- lerhifi.seq[[i]]
  colcen.wins <- seq(1,length(colcen.chrseq),by=10000)
  wincolcen.allcos <- NULL
  wincolcen.syrisnps <- NULL
  wincolcen.syrihdr <- NULL
  wincolcen.genes <- NULL
  wincolcen.cen178plus <- NULL
  wincolcen.cen178minus <- NULL
  for(j in 1:(length(colcen.wins)-1)){
    print(j)
    wincolcen.allcos <- c(wincolcen.allcos,length(which(colcen.chrallcos[,4]>colcen.wins[j]&colcen.chrallcos[,4]<colcen.wins[j+1])))
    wincolcen.syrisnps <- c(wincolcen.syrisnps,length(which(as.numeric(chr.snps[,2])>colcen.wins[j]&as.numeric(chr.snps[,2])<colcen.wins[j+1])))
    wincolcen.syrihdr <- c(wincolcen.syrihdr,length(which(as.numeric(chr.hdr[,2])>colcen.wins[j]&as.numeric(chr.hdr[,2])<colcen.wins[j+1])))
    wincolcen.genes <- c(wincolcen.genes,length(which(colcen.chrgenes[,4]>colcen.wins[j]&colcen.chrgenes[,4]<colcen.wins[j+1])))
    wincolcen.cen178plus <- c(wincolcen.cen178plus,length(which(colcen.chrplus[,4]>colcen.wins[j]&colcen.chrplus[,4]<colcen.wins[j+1])))
    wincolcen.cen178minus <- c(wincolcen.cen178minus,length(which(colcen.chrminus[,4]>colcen.wins[j]&colcen.chrminus[,4]<colcen.wins[j+1])))
  }
  cenlim.start <- all.left.lrz[i]-10000
  cenlim.end <- all.right.lrz[i]+10000
  cen.lim <- c(cenlim.start,cenlim.end)
  print(i)
  plot(colcen.wins,c(wincolcen.genes,0),type='h',xlim=cen.lim,col='dark green',ylim=c(0,10))
  #par(new=T)
  #plot(colcen.wins,c(wincolcen.cen178plus,0),type='h',col='dark gray',yaxt='n',xlim=cen.lim,main='cos + SPO11-1 + genes',ylim=c(0,70))
  #par(new=T)
  #plot(colcen.wins,c(wincolcen.cen178minus,0),type='h',col='black',yaxt='n',xlim=cen.lim,ylim=c(0,70))
  par(new=T)
  plot(colcen.wins,c(wincolcen.allcos,0),type='h',xlim=cen.lim,ylim=c(0,20),yaxt='n',col='purple')
  par(new=T)
  plot(chrchrom.dat[,2],chrchrom.dat$log2_WT_SPO11oligos_Rep1_WT_gDNA_Rep1_R1,type='h',col='red',ylim=c(0,0.85),xlim=cen.lim,yaxt='n')
  axis(side=4)
  abline(v=c(all.left.lrz[i],all.right.lrz[i]),col=6,lwd=2)
  abline(v=c(all.left.nrz[i],all.right.nrz[i]),col=1,lwd=2)
  dat.lrz.genes <- c(dat.lrz.genes,wincolcen.genes[which(colcen.wins>all.left.lrz[i]&colcen.wins<all.right.lrz[i])])
  spo11 <- chrchrom.dat$log2_WT_SPO11oligos_Rep1_WT_gDNA_Rep1_R1
  dat.lrz.spo11 <- c(dat.lrz.spo11,spo11[which(colcen.wins>all.left.lrz[i]&colcen.wins<all.right.lrz[i])])
  dat.lrz.cos <- c(dat.lrz.cos,wincolcen.allcos[which(colcen.wins>all.left.lrz[i]&colcen.wins<all.right.lrz[i])])
}

###############################################
# correlation of values in CTL3.9 map windows #
###############################################

data <- read.csv(file='map_data_ColCEN.csv')
ctl.markers <- data[,2]
ctl.cMMb <- data[,5]

chrom.dat <- read.table(file='ChIP_input_log2ChIPinput_RaGOO_v2.0_genome_norm_coverage_matrix_10kb_unsmoothed.tsv',header=T)
chrom.chr <- chrom.dat[which(chrom.dat[,1]=='Chr3'),]

colcen.genes <- read.table('ColCEN_GENES_TAIR10.gff3',header=F)
colcen.genes <- colcen.genes[which(colcen.genes[,3]=='gene'),]
colcen.genes <- colcen.genes[which(colcen.genes[,1]=='Chr3'),]

colcen.edta <- read.table(file='t2t-col.20201227.fasta.mod.EDTA.TEanno.gff3')
colcen.edta <- colcen.edta[which(colcen.edta[,1]=='Chr3'),]
colcen.gypsy <- colcen.edta[which(colcen.edta[,3]=='Gypsy_LTR_retrotransposon'),]

deep.cg <- read.delim(file='Deepsignal_30kb90%_CG_call_mods_frequency_T2T_100621.tsv',header=F)
deep.cg <- deep.cg[which(deep.cg[,1]=='Chr3'),]
deep.chg <- read.delim(file='Deepsignal_30kb90%_CHG_call_mods_frequency_T2T_100621.tsv',header=F)
deep.chg <- deep.chg[which(deep.chg[,1]=='Chr3'),]
deep.chh <- read.delim(file='Deepsignal_30kb90%_CHH_call_mods_frequency_T2T_100621.tsv',header=F)
deep.chh <- deep.chh[which(deep.chh[,1]=='Chr3'),]

ctl.cg <- NULL
ctl.chg <- NULL
ctl.chh <- NULL
ctl.genes <- NULL
ctl.gypsy <- NULL
ctl.h2az <- NULL
ctl.h3k4 <- NULL
ctl.h3k9 <- NULL
ctl.rec8 <- NULL
ctl.asy1 <- NULL
ctl.cenh3 <- NULL
ctl.spo11wt <- NULL
for(k in 1:length(ctl.markers)){
  print(k)
  ctl.cg <- c(ctl.cg,mean(deep.cg[which(deep.cg[,2]>ctl.markers[k]&deep.cg[,2]<=ctl.markers[k+1]),10]))  
  ctl.chg <- c(ctl.chg,mean(deep.chg[which(deep.chg[,2]>ctl.markers[k]&deep.chg[,2]<=ctl.markers[k+1]),10]))  
  ctl.chh <- c(ctl.chh,mean(deep.chh[which(deep.chh[,2]>ctl.markers[k]&deep.chh[,2]<=ctl.markers[k+1]),10]))
  ctl.genes <- c(ctl.genes,length(which(colcen.genes[,4]>ctl.markers[k]&colcen.genes[,4]<=ctl.markers[k+1])))
  ctl.gypsy <- c(ctl.gypsy,length(which(colcen.gypsy[,4]>ctl.markers[k]&colcen.gypsy[,4]<=ctl.markers[k+1])))
  ctl.h2az <- c(ctl.h2az,mean(chrom.chr[which(chrom.chr[,2]>ctl.markers[k]&chrom.chr[,2]<=ctl.markers[k+1]),18])) 
  ctl.h3k4 <- c(ctl.h3k4,mean(chrom.chr[which(chrom.chr[,2]>ctl.markers[k]&chrom.chr[,2]<=ctl.markers[k+1]),36])) 
  ctl.h3k9 <- c(ctl.h3k9,mean(chrom.chr[which(chrom.chr[,2]>ctl.markers[k]&chrom.chr[,2]<=ctl.markers[k+1]),60])) 
  ctl.rec8 <- c(ctl.rec8,mean(chrom.chr[which(chrom.chr[,2]>ctl.markers[k]&chrom.chr[,2]<=ctl.markers[k+1]),81])) 
  ctl.asy1 <- c(ctl.asy1,mean(chrom.chr[which(chrom.chr[,2]>ctl.markers[k]&chrom.chr[,2]<=ctl.markers[k+1]),84])) 
  ctl.cenh3 <- c(ctl.cenh3,mean(chrom.chr[which(chrom.chr[,2]>ctl.markers[k]&chrom.chr[,2]<=ctl.markers[k+1]),6])) 
  ctl.spo11wt <- c(ctl.spo11wt,mean(chrom.chr[which(chrom.chr[,2]>ctl.markers[k]&chrom.chr[,2]<=ctl.markers[k+1]),99]))
}



bp.copia <- c(bp.copia,sum(copia.width[which(copia[,4]>colcen.wins[j]&copia[,4]<colcen.wins[j+1])]))
bp.gypsy <- c(bp.gypsy,sum(gypsy.width[which(gypsy[,4]>colcen.wins[j]&gypsy[,4]<colcen.wins[j+1])]))

par(mfcol=c(5,2))
par(mar=c(2,2,2,2))
plot(ctl.cg,ctl.cMMb,main='CG CHG CHH',col='red',xlim=c(0,1))
par(new=T)
plot(ctl.chg,ctl.cMMb,col='blue',xlim=c(0,1))
par(new=T)
plot(ctl.chh,ctl.cMMb,col='green',xlim=c(0,1))
plot(ctl.gypsy,ctl.cMMb,main='gypsy',col='blue')
plot(ctl.h3k9,ctl.cMMb,main='h3k9',col='blue')
plot(ctl.cenh3,ctl.cMMb,main='cenh3',col='blue')
plot(ctl.rec8,ctl.cMMb,main='rec8',col='blue')
plot(ctl.asy1,ctl.cMMb,main='asy1',col='blue')
plot(ctl.h2az,ctl.cMMb,main='h2az',col='red')
plot(ctl.genes,ctl.cMMb,main='genes',col='red')
plot(ctl.h3k4,ctl.cMMb,main='h3k4',col='red')
plot(ctl.spo11wt,ctl.cMMb,main='spo11',col='red')

cor.test(ctl.cg,ctl.cMMb)
-0.5529176 p-value = 6.244e-09
cor.test(ctl.chg,ctl.cMMb)
-0.5517173 6.837e-09
cor.test(ctl.chh,ctl.cMMb)
-0.4894811  4.808e-07

cor.test(ctl.h3k9,ctl.cMMb)
-0.537668 1.927e-08

cor.test(ctl.rec8,ctl.cMMb)
-0.559852  3.67e-09
cor.test(ctl.asy1,ctl.cMMb)
-0.5316429 2.962e-08

cor.test(ctl.gypsy,ctl.cMMb)
-0.3338206 0.000888
cor.test(ctl.cenh3,ctl.cMMb)
-0.1332343 0.198

cor.test(ctl.genes,ctl.cMMb)
0.222659 0.02922
cor.test(ctl.spo11wt,ctl.cMMb)
0.4938681 3.66e-07
cor.test(ctl.h3k4,ctl.cMMb)
0.50958 1.335e-07

cor.test(ctl.h2az,ctl.cMMb)
0.5795472 7.576e-10

hs.code <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,2,0,0,0,2,2,0,2,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0)
data.add <- cbind(hs.code,data,ctl.cg,ctl.chg,ctl.chh,ctl.h3k9,ctl.h3k4,ctl.h2az,ctl.asy1,ctl.rec8,ctl.cenh3,ctl.spo11wt)
data.add <- data.add[-length(data.add[,1]),]
data.hot <- data.add[which(data.add[,1]==1),]
data.cold <- data.add[which(data.add[,1]==2),]
data.rest <- data.add[which(data.add[,1]==0),]

par(mfcol=c(2,5))
par(mar=c(2,2,2,2))

plot(rep(1,length(data.hot[,1])),data.hot[,13],ylim=c(0,1),xlim=c(0,4),col='red',main='CG')
par(new=T)
plot(rep(2,length(data.rest[,1])),data.rest[,13],ylim=c(0,1),xlim=c(0,4))
par(new=T)
plot(rep(3,length(data.cold[,1])),data.cold[,13],ylim=c(0,1),xlim=c(0,4),col='blue')

wilcox.test(data.hot[,13],data.rest[,13])
0.07256
wilcox.test(data.hot[,13],data.cold[,13])
0.001998

plot(rep(1,length(data.hot[,1])),data.hot[,17],ylim=c(-1.4,0.5),xlim=c(0,4),col='red',main='H3K4me3')
par(new=T)
plot(rep(2,length(data.rest[,1])),data.rest[,17],ylim=c(-1.4,0.5),xlim=c(0,4))
par(new=T)
plot(rep(3,length(data.cold[,1])),data.cold[,17],ylim=c(-1.4,0.5),xlim=c(0,4),col='blue')

wilcox.test(data.hot[,17],data.rest[,17])
0.0286
wilcox.test(data.hot[,17],data.cold[,17])
0.001998

plot(rep(1,length(data.hot[,1])),data.hot[,14],ylim=c(0,0.6),xlim=c(0,4),col='red',main='CHG')
par(new=T)
plot(rep(2,length(data.rest[,1])),data.rest[,14],ylim=c(0,0.6),xlim=c(0,4))
par(new=T)
plot(rep(3,length(data.cold[,1])),data.cold[,14],ylim=c(0,0.6),xlim=c(0,4),col='blue')

wilcox.test(data.hot[,14],data.rest[,14])
0.0667
wilcox.test(data.hot[,14],data.cold[,14])
0.001998

plot(rep(1,length(data.hot[,1])),data.hot[,18],ylim=c(-1.7,0.6),xlim=c(0,4),col='red',main='H2AZ')
par(new=T)
plot(rep(2,length(data.rest[,1])),data.rest[,18],ylim=c(-1.7,0.6),xlim=c(0,4))
par(new=T)
plot(rep(3,length(data.cold[,1])),data.cold[,18],ylim=c(-1.7,0.6),xlim=c(0,4),col='blue')

wilcox.test(data.hot[,18],data.rest[,18])
0.03801
wilcox.test(data.hot[,18],data.cold[,18])
0.001998

plot(rep(1,length(data.hot[,1])),data.hot[,15],ylim=c(0,0.14),xlim=c(0,4),col='red',main='CHH')
par(new=T)
plot(rep(2,length(data.rest[,1])),data.rest[,15],ylim=c(0,0.14),xlim=c(0,4))
par(new=T)
plot(rep(3,length(data.cold[,1])),data.cold[,15],ylim=c(0,0.14),xlim=c(0,4),col='blue')

wilcox.test(data.hot[,15],data.rest[,15])
0.06302
wilcox.test(data.hot[,15],data.cold[,15])
0.007992

plot(rep(1,length(data.hot[,1])),data.hot[,16],ylim=c(-0.8,2.2),xlim=c(0,4),col='red',main='H3K9me2')
par(new=T)
plot(rep(2,length(data.rest[,1])),data.rest[,16],ylim=c(-0.8,2.2),xlim=c(0,4))
par(new=T)
plot(rep(3,length(data.cold[,1])),data.cold[,16],ylim=c(-0.8,2.2),xlim=c(0,4),col='blue')

wilcox.test(data.hot[,16],data.rest[,16])
0.05616
wilcox.test(data.hot[,16],data.cold[,16])
0.003996

plot(rep(1,length(data.hot[,1])),data.hot[,19],ylim=c(-0.5,0.8),xlim=c(0,4),col='red',main='ASY1')
par(new=T)
plot(rep(2,length(data.rest[,1])),data.rest[,19],ylim=c(-0.5,0.8),xlim=c(0,4))
par(new=T)
plot(rep(3,length(data.cold[,1])),data.cold[,19],ylim=c(-0.5,0.8),xlim=c(0,4),col='blue')

wilcox.test(data.hot[,19],data.rest[,19])
0.1084
wilcox.test(data.hot[,19],data.cold[,19])
0.001998

plot(rep(1,length(data.hot[,1])),data.hot[,20],ylim=c(-0.6,1.0),xlim=c(0,4),col='red',main='REC8')
par(new=T)
plot(rep(2,length(data.rest[,1])),data.rest[,20],ylim=c(-0.6,1.0),xlim=c(0,4))
par(new=T)
plot(rep(3,length(data.cold[,1])),data.cold[,20],ylim=c(-0.6,1.0),xlim=c(0,4),col='blue')

wilcox.test(data.hot[,20],data.rest[,20])
0.04847
wilcox.test(data.hot[,20],data.cold[,20])
0.001998

plot(rep(1,length(data.hot[,1])),data.hot[,21],ylim=c(-1.1,0.4),xlim=c(0,4),col='red',main='CENH3')
par(new=T)
plot(rep(2,length(data.rest[,1])),data.rest[,21],ylim=c(-1.1,0.4),xlim=c(0,4))
par(new=T)
plot(rep(3,length(data.cold[,1])),data.cold[,21],ylim=c(-1.1,0.4),xlim=c(0,4),col='blue')

wilcox.test(data.hot[,21],data.rest[,21])
0.4734
wilcox.test(data.hot[,21],data.cold[,21])
0.007992

plot(rep(1,length(data.hot[,1])),data.hot[,22],ylim=c(-1.0,0.6),xlim=c(0,4),col='red',main='SPO11-1-oligos')
par(new=T)
plot(rep(2,length(data.rest[,1])),data.rest[,22],ylim=c(-1.0,0.6),xlim=c(0,4))
par(new=T)
plot(rep(3,length(data.cold[,1])),data.cold[,22],ylim=c(-1.0,0.6),xlim=c(0,4),col='blue')

wilcox.test(data.hot[,22],data.rest[,22])
0.06124
wilcox.test(data.hot[,22],data.cold[,22])
0.001998


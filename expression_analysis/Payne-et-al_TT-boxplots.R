## Payne-et-al_TT-boxplots.R
##
## Script used to plot expression boxplots in Payne et al 2021
##
## cyp I-2022

## xbirch reference
brain_data <- read.csv("input_files/TT-brain-22c-xbirch-gtf_DGE_lfc-shr_all.csv",header=T)
liver_data <- read.csv("input_files/TT-liver-22c-xbirch-gtf_DGE_lfc-shr_all.csv",header=T)

## visual settings
malcol_22=rgb(0/255,0/255,175/255)
hetcol_22=rgb(100/255,0/255,175/255)
bircol_22=rgb(150/255,0/255,0/255)
malcol_33=rgb(50/255,150/255,255/255)
hetcol_33=rgb(200/255,0/255,175/255)
bircol_33=rgb(255/255,0/255,0/255)

## error bar function
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

## plotting function
TT_boxplot <- function(gene, tissue, data2plot, ymin, ymax, legend_pos_x, legend_pos_y, save_plot) {
  # pull expression data
  if(tissue=="brain") {
    n<-subset(brain_data,Gene==gene)
    mal_22c <-c(n$S1,n$S25,n$S48)
    bir_22c <-c(n$S5,n$S29,n$S52)
    f1_22c  <-c(n$S9,n$S33)
    mal_33c <-c(n$S13,n$S36,n$S60)
    bir_33c <-c(n$S17,n$S40,n$S64)
    f1_33c  <-c(n$S21,n$S44,n$S68)
  } else if(tissue=="liver") {
    n <- subset(liver_data,Gene==gene)
    mal_22c <- c(n$S4,n$S28,n$S51)
    bir_22c  <- c(n$S8,n$S32)
    f1_22c   <- c(n$S12,n$S59)
    mal_33c <- c(n$S16,n$S63)
    bir_33c  <- c(n$S20,n$S43,n$S67)
    f1_33c   <- c(n$S24,n$S47,n$S71)
  }
  
  if(save_plot) { 
    outfile = paste0(data2plot,"_TT-exp-plot_",tissue,"_",n$Gene,".pdf")
    pdf(outfile,width=4,height=6) 
  }
  if(data2plot == "both") {
    plot(0.4:2.4,c(mean(mal_22c),mean(f1_22c),mean(bir_22c)),col=c(malcol_22,hetcol_22,bircol_22),ylim=c(ymin,ymax),xlab="Genotype",ylab="Normalized counts",xaxt="n",pch=20,cex=2,xlim=c(0.1,3))
    error.bar(0.4:2.4,c(mean(mal_22c),mean(f1_22c),mean(bir_22c)),c(sd(mal_22c),sd(f1_22c),sd(bir_22c)),col=c(malcol_22,hetcol_22,bircol_22),lwd=2)
    #BB_22.5C<-bir_22c; MB_22.5C<-f1_22c; MM_22.5C<-mal_22c;
    points(0.7:2.7,c(mean(mal_33c),mean(f1_33c),mean(bir_33c)),col=c(malcol_33,hetcol_33,bircol_33),pch=18,cex=1.5)
    error.bar(0.7:2.7,c(mean(mal_33c),mean(f1_33c),mean(bir_33c)),c(sd(mal_33c),sd(f1_33c),sd(bir_33c)),col=c(malcol_33,hetcol_33,bircol_33),lwd=2)
    mtext(c("MM","MB","BB"),at=0.6:2.6,side=1,line=1)
    # legend
    lcol=c("black","dark gray")
    legend(legend_pos_x, legend_pos_y, legend=c(expression(paste("22.5",degree,"C")),expression(paste("33.5",degree,"C"))), col=lcol, pch = c(20, 18), lty = c(1, 1), pt.cex = c(1.8, 1.4), cex=1)
  } else { 
    if(data2plot == "22") {
      BB<-bir_22c; MB<-f1_22c; MM<-mal_22c;
      malcol=rgb(0/255,0/255,175/255,alpha=0.3)
      hetcol=rgb(100/255,0/255,175/255,alpha=0.6)
      bircol=rgb(150/255,0/255,0/255,alpha=0.6)
      plot(1:3,c(mean(MM),mean(MB),mean(BB)),col=c(malcol_22,hetcol_22,bircol_22),ylim=c(ymin,ymax),xlab="Genotype",ylab="Normalized counts",xaxt="n",pch=20,cex=2,xlim=c(0.5,3.5))
      error.bar(1:3,c(mean(MM),mean(MB),mean(BB)),c(sd(MM),sd(MB),sd(BB)),col=c(malcol_22,hetcol_22,bircol_22),lwd=2)
      
      # add transparent points
      noise<-runif(length(BB),0.2,0.3)
      points(rep(3,length(BB))+noise,BB,pch=20,cex=1.8,col=bircol)
      noise<-runif(length(MB),0.2,0.3)
      points(rep(2,length(MB))+noise,MB,pch=20,cex=1.8,col=hetcol)
      noise<-runif(length(MM),0.2,0.3)
      points(rep(1,length(MM))+noise,MM,pch=20,cex=1.8,col=malcol)
    } else { 
      BB<-bir_33c; MB<-f1_33c; MM<-mal_33c;
      malcol=rgb(50/255,150/255,255/255,alpha=0.3)
      hetcol=rgb(200/255,0/255,175/255,alpha=0.6)
      bircol=rgb(255/255,0/255,0/255,alpha=0.6)
      plot(1:3,c(mean(MM),mean(MB),mean(BB)),col=c(malcol_33,hetcol_33,bircol_33),ylim=c(ymin,ymax),xlab="Genotype",ylab="Normalized counts",xaxt="n",pch=18,cex=2,xlim=c(0.5,3.5))
      error.bar(1:3,c(mean(MM),mean(MB),mean(BB)),c(sd(MM),sd(MB),sd(BB)),col=c(malcol_33,hetcol_33,bircol_33),lwd=2)
    
      # add transparent points
      noise<-runif(length(BB),0.2,0.3)
      points(rep(3,length(BB))+noise,BB,pch=18,cex=1.8,col=bircol)
      noise<-runif(length(MB),0.2,0.3)
      points(rep(2,length(MB))+noise,MB,pch=18,cex=1.8,col=hetcol)
      noise<-runif(length(MM),0.2,0.3)
      points(rep(1,length(MM))+noise,MM,pch=18,cex=1.8,col=malcol)
    }
    
    mtext(c("MM","MB","BB"),at=1:3,side=1)
  }
  
  title(n$Gene)

  if(save_plot) { dev.off() }
}

## paper plots
TT_boxplot("g4364", "liver", "33", 0, 525, 0.1, 515, TRUE) # tnfaip3
TT_boxplot("g4429", "brain", "22", 390, 740, 0, 0, TRUE) # p4ha1
TT_boxplot("g15226", "liver", "both", 180, 1360, 0.1, 1370, FALSE) # nfkbia
TT_boxplot("g15263", "brain", "22", 715, 975, 0, 0, FALSE) # snw1
TT_boxplot("g4403", "liver", "both", 25, 160, 1.8, 160, FALSE) # sf3b5
#TT_boxplot("g4403", "liver", "both", 0, 150, 1.8, 160, FALSE) # sf3b5
TT_boxplot("g4403", "liver", "both", 25, 160, 1.8, 160, FALSE) # sf3b5


TT_boxplot("g252", "brain", "both", 1000, 2000, 0, 0, FALSE) # p4ha1
#TT_boxplot("g8451", "liver", "both", 0, 100, 0, 0, FALSE) # p4ha1
TT_boxplot("g19980", "brain", "both", 1500, 3000, 0, 0, TRUE) # xbp1
TT_boxplot("g19980", "liver", "both", 6000, 31500, 0.1, 31500, FALSE) # xbp1
TT_boxplot("g28231", "brain", "both", 0, 500, 0.1, 31500, FALSE) # nfkb1

TT_boxplot("g15226", "liver", "both", 180, 1360, 0.1, 1380, TRUE) # nfkbia
TT_boxplot("g4364", "liver", "both", 0, 510, 0.1, 520, TRUE) # tnfaip3
#TT_boxplot("g15263", "brain", "both", 715, 975, 1.8, 980, FALSE) # snw1
#TT_boxplot("g19980", "brain", "both", 1500, 3000, 1.8, 3030, FALSE) # xbp1
TT_boxplot("g15263", "brain", "both", 715, 975, 0,0, TRUE) # snw1
TT_boxplot("g19980", "brain", "both", 1500, 3000, 0,0, TRUE) # xbp1


TT_boxplot("g4383", "brain", "both", 8, 100, 0.1, 1380, FALSE) # zfp62
TT_boxplot("g4403", "liver", "both", 28, 158, 0.1, 1380, TRUE) # sf3b5
TT_boxplot("g4424", "brain", "both", 395, 780, 0.1, 1380, TRUE) # ndnf
TT_boxplot("g4395", "liver", "both", 75, 420, 0.1, 1380, TRUE) # ifngr1l
TT_boxplot("g30037", "brain", "both", 490, 2310, 0, 0, TRUE) # nr1d2a
TT_boxplot("g15351", "brain", "both", 27, 170, 0.1, 1380, TRUE) # cipcb

TT_boxplot("g8984", "brain", "both", 1200, 2930, 0, 0, TRUE) # nr1d2b



## all genes of interest
gene_list <- c("g4383") #zfp62 in liver
gene_list <- c("g4403") #sf3b5 in both
gene_list <- c("g4403","g15196") # sf3b5, sf3b6
gene_list <- c("g4424") # ndnf in brain
gene_list <- c("g4429") # p4ha1 in brain * = y: 335,715 - legend(0.1, 715 - pdf: 4.3x6.6
gene_list <- c("g4395") # ifngr1l in liver
gene_list <- c("g4364") # tnfaip3 in liver * = y: 0,510 - legend(0.1, 515
gene_list <- c("g4364", "g15226") # tnfaip3, nfkbia in liver
gene_list <- c("g11253","g29988","g28231", "g5161", "g15733", "g4364", "g15226") # p65, tnfb, nfkb1, nfkb2, nfkbie, tnfaip3, nfkbia in liver
gene_list <- c("g15248") # nrxn3b in brain
gene_list <- c("g16226") # nr1d1 in brain
gene_list <- c("g521","g13770") # roraa, rorab in brain
gene_list <- c("g18195") # arntl2 in brain
TT_boxplot("g4371", "brain", "both", 0, 1000, 0.1, 31500, FALSE) # akt3
TT_boxplot("g8984", "brain", "both", 1000, 3100, 0,0, FALSE) # nr1d2b
TT_boxplot("g18195", "brain", "both", 0, 300, 0,0, FALSE) # arntl2
TT_boxplot("g900", "brain", "both", 0, 1000, 0,0, FALSE) # arntl1a



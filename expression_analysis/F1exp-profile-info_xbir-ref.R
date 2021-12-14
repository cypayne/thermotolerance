## Outputs DGE csv file with information about F1 expression profile
## used to find genes with significantly transgressive response in F1s
## input file must have lfc and padj values for both parents and F1s
## cyp V-2020

## XXX Configuration:
infile  <- "input_files/TT-brain-22c-xbirch-gtf_DGE_lfc-shr_all.csv"
#infile  <- "input_files/TT-brain-33c-xbirch-gtf_DGE_lfc-shr_all.csv"
#infile  <- "input_files/TT-liver-22c-xbirch-gtf_DGE_lfc-shr_all.csv"
#infile  <- "input_files/TT-liver-33c-xbirch-gtf_DGE_lfc-shr_all.csv"

outfile <- paste(infile,"_with-F1info.csv",sep="")
thresh  <- 0.05 # set padj threshold
cond    <- "33" # set the condition: 22c="22", 33c="33", 33cV22c="delta"

data<-read.csv(infile,sep=",",head=TRUE)
#nuc <- data[1:19106,1:length(colnames(data))]
#mito <- data[19106:19143,1:length(colnames(data))]
#data <- mito

## Function: test_transgression(x) asigns an F1 expression profile: 
## high=          sig higher than both parents
## low=           sig lower than both parents
## intermediate=  intermediate, sig diff from one or both parents
test_transgression <- function(x,...) {
  f1trans.high = 0
  f1trans.low  = 0
  f1exp        = "NA"
  # pull lfc and padj info
  if(cond == "22") {
    bir.lfc   <- x$LFC_res.22CbirVF1  # f1 higher if -
    bir.padj  <- x$padj_res.22CbirVF1 
    mal.lfc   <- x$LFC_res.22CF1Vmal  # f1 higher if +
    mal.padj  <- x$padj_res.22CF1Vmal 
  } else if(cond == "33") {
    bir.lfc   <- x$LFC_res.33CbirVF1  # f1 higher if -
    bir.padj  <- x$padj_res.33CbirVF1 
    mal.lfc   <- x$LFC_res.33CF1Vmal  # f1 higher if +
    mal.padj  <- x$padj_res.33CF1Vmal 
  } else if(cond == "delta") {
    bir.lfc   <- x$LFC_res.birVF1.33cV22c  # f1 higher if -
    bir.padj  <- x$padj_res.birVF1.33cV22c 
    mal.lfc   <- x$LFC_res.F1Vmal.33cV22c  # f1 higher if +
    mal.padj  <- x$padj_res.F1Vmal.33cV22c 
  }
  
  # debug: sprintf("bir-f1: lfc %f, padj %f; f1-mal: lfc %f, padj %f", bir.lfc, bir.padj, mal.lfc, mal.padj)
  
  # if padj values are NA, skip and set output vals to NA
  if(!is.na(bir.padj) && !is.na(mal.padj)) { 
    # only interested in bir-F1 and F1-mal comparisons that are significant
    # i.e. padj less than padj threshold
    if(bir.padj<thresh && mal.padj<thresh) {
      if(bir.lfc<0 && mal.lfc>0) {
        # when f1 upreg: bir-F1 = -, f1-mal = +
        f1trans.high <- 1 
        f1exp <- "high"
      } else if(bir.lfc>0 && mal.lfc<0) {
        # when f1 downreg: bir-F1 = +, f1-mal = -
        f1trans.low <- 1
        f1exp <- "low"
      } else { f1exp <- "intermediate" }
    }
  } else { f1trans.high = "NA"; f1trans.low = "NA" }
  
  # return output row to be added to final csv
  new.row <- c(f1exp,f1trans.high,f1trans.low)
  return(new.row)
}

int <- lapply(1:nrow(data), function(x) test_transgression(data[x,]) )
transtest.out <- do.call(rbind,int)
# run test_transgression function on all rows of data
#debug: transtest.out <- do.call(rbind, lapply(1:nrow(data), function(x) test_transgression(data[x,]) ))
tt.df <- as.data.frame(transtest.out)
colnames(tt.df) <- c("sig-F1_exp_profile","F1_transgress.high","F1_transgress.low")
# create combined csv with new information
new.data <- data.frame(data,tt.df)
write.csv(new.data,outfile)

high.22c <- new.data[new.data$F1_transgress.high == 1,]
nrow(high.22c)
low.22c <- new.data[new.data$F1_transgress.low == 1,]
nrow(low.22c)

high.33c <- new.data[new.data$F1_transgress.high == 1,]
nrow(high.33c)
low.33c <- new.data[new.data$F1_transgress.low == 1,]
nrow(low.33c)

### extra exploration
## brain 22C: 
##    thresh=0.05: 132 high, 159 low
##    thresh=0.10: 246 high, 275 low
## brain 33C: 
##    thresh=0.05: 66 high, 24 low
##    thresh=0.10: 92 high, 52 low
## liver 22C: 
##    thresh=0.05: 2 high, 5 low
##    thresh=0.10: 3 high, 7 low
## liver 33C: 
##    thresh=0.05: 38 high, 16 low
##    thresh=0.10: 59 high, 37 low

## overlapping genes between sets (padj0.1)
# brain high 33c and 22c: "g14782" "g23748" "g23749" "g27396" "g28186" "g29149" "g29240" "g359"   "g451" 
# liver high 33c and 22c: none
overlap <- intersect(high.22c$Gene,high.33c$Gene)
overlap
# brain low 33c and 22c: "g14795" "g15411" "g19249" "g20044"
# liver low 33c and 22c: "g23367"
overlap <- intersect(low.22c$Gene,low.33c$Gene)
overlap
# brain low 22c, high 33c: "g1871"  "g21588" "g30342" "g728"
# brain high 22c, low 33c: none
# liver low 22c, high 33c: none
# liver high 22c, low 33c: none
overlap <- intersect(low.22c$Gene,high.33c$Gene)
overlap
overlap <- intersect(high.22c$Gene,low.33c$Gene)
overlap


## overlapping genes between sets (padj0.05)
# brain high 33c and 22c: "g14782" "g23748" "g23749" "g29240" "g359"   "g451"
# liver high 33c and 22c: none
overlap <- intersect(high.22c$Gene,high.33c$Gene)
overlap
# brain low 33c and 22c: "g14795" "g20044"
# liver low 33c and 22c: "g23367"
overlap <- intersect(low.22c$Gene,low.33c$Gene)
overlap
# brain low 22c, high 33c: "g21588" "g728" 
# brain high 22c, low 33c: none
# liver low 22c, high 33c: none
# liver high 22c, low 33c: none
overlap <- intersect(low.22c$Gene,high.33c$Gene)
overlap
overlap <- intersect(high.22c$Gene,low.33c$Gene)
overlap


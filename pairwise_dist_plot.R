library("spider")
library("ape")
library("seqinr")
library("gridExtra")
library(hexbin)
library(grid)
library(ggplot2)
library(sjstats)
library(gplots)

# Function plots pairwise nucleotide distance comparison plot.
# Each dot corresponds to a pair of nucleotide distances between 
# the same pair of genomes in two genomic regions - odd sites and even sites of genome (see axis).
# Returns list with ggplot, matrices of distances between pairs of seqences calculated for 
# odd and even sites

plot_control = function(dna_object){
  
  # subalignments for odd and even sites
  al_odd=dna_object[1:length(dna_object[,1]), seq(from = 1, to = length(dna_object[1,]), by=2)]
  al_even=dna_object[1:length(dna_object[,1]), seq(from = 2, to = length(dna_object[1,]), by=2)]
  
  # distance matrices for each region
  dna_sl_dist1 <-dist.gene(al_odd, method = "percentage",  pairwise.deletion = TRUE)
  dna_sl_dist2 <-dist.gene(al_even, method = "percentage",  pairwise.deletion = TRUE)
  
  # adding random noise to distances matrices' values
  dist1= as.vector(dna_sl_dist1) + rnorm(length(dna_sl_dist1),mean = 0,sd= 0.001)
  dist2= as.vector(dna_sl_dist2) + rnorm(length(dna_sl_dist2),mean = 0,sd= 0.0001)
  
  # pairwise nucleotide distance comparison plot
  dist_plot=ggplot(data.frame(dist1,dist2),aes(dist1,dist2))+stat_bin2d(binwidth = 0.001)+ 
    scale_fill_gradientn(colours=c("blue","red"))+ theme(legend.justification=c(1,0), legend.position=c(1,0))+
    xlab("odd sites")+ylab("even sites")
  #+  geom_smooth(method='lm',formula=y~x)
  
  return(list(dist_plot, dna_sl_dist1, dna_sl_dist2))
  
}

# Function plots pairwise nucleotide distance comparison plot.
# Each dot corresponds to a pair of nucleotide distances between 
# the same pair of genomes in two genomic regions - st1-e1 and st2-e2 (see axis).
# Returns list with ggplot, matrices of distances between pairs of seqences calculated for 
# st1-e1 and st2-e2 regions

plot_dist_test = function(dna_object, st1,e1,st2,e2){
  
  # subalignments for st1-e1 and st2-e2 regions
  dna_sl1=dna_object[1:length(dna_object[,1]), seq(from = st1, to = e1, by=1)]
  dna_sl2=dna_object[1:length(dna_object[,1]), seq(from = st2, to = e2, by=1)]
  
  # distance matrices for each region
  dna_sl_dist1 <-dist.gene(dna_sl1, method = "percentage",  pairwise.deletion = TRUE)
  dna_sl_dist2 <-dist.gene(dna_sl2, method = "percentage",  pairwise.deletion = TRUE)
  #HepadnaDist1 <-dist.dna(Hepadna1, as.matrix = TRUE, model = "JC69")
  #HepadnaDist4 <-dist.dna(Hepadna4, as.matrix = TRUE, model = "JC69")
  
  # adding random noise to distances matrices' values
  dist1= as.vector(dna_sl_dist1) + rnorm(length(dna_sl_dist1),mean = 0,sd= 0.001)
  dist2= as.vector(dna_sl_dist2) + rnorm(length(dna_sl_dist2),mean = 0,sd= 0.001)
  
  #Hepadna1=abs(rnorm(HepadnaDist1[lower.tri(HepadnaDist1,diag = FALSE)], mean=HepadnaDist1[lower.tri(HepadnaDist1,diag = FALSE)], sd = 0))
  #Hepadna4=abs(rnorm(HepadnaDist4[lower.tri(HepadnaDist4,diag = FALSE)], mean=HepadnaDist4[lower.tri(HepadnaDist4,diag = FALSE)], sd = 0))
  
  
  # pairwise nucleotide distance comparison plot
  dist_plot=ggplot(data.frame(dist1,dist2),aes(dist1,dist2))+stat_bin2d(binwidth = 0.003)+ 
    scale_fill_gradientn(colours=c("blue","red"))+ theme(legend.justification=c(1,0), legend.position=c(1,0))+
    xlab(paste(toString(st1),toString(e1),sep=":"))+ylab(paste(toString(st2),toString(e2),sep=":"))
  #+  geom_smooth(method='lm',formula=y~x)
  
  return(list(dist_plot, dna_sl_dist1, dna_sl_dist2))
  
}

# Prints names of sequences pairs which pairwise nucleotide distances lay between 
# val11-val12 in genomic region 1 and val21-val22 in genomic region 2/
# distM1 - distance matrix built for region 1
# distM2 - distance matrix built for region 2
# Prints strings "name1\tname2\tdistance_in_region1\tdistance_in_region2"
find_recomb_names <- function(distM1, val11, val12, distM2, val21, val22){
  
  #positions of matrix for region 1 with values between val11 and val12
  b1 = find_dist_slice(distM1,val11, val12)
  #positions of matrix for region 2 with values between val21 and val22
  b2 = find_dist_slice(distM2,val21, val22)
  #intersection of b1 and b2
  b= intersect(b1,b2)
  len_b = length(distM1[1,])
  
  #returns positions of rows in matrix
  r = sapply(b,function(x){if(x%%len_b==0){return(len_b)} else{return(x%%len_b)}})
  #returns positions of columns in matrix
  c = sapply(b,function(x){if(x%%len_b==0){return(x%/%len_b)} else{return(x%/%len_b+1)}})
  
  #values of rows and columns
  c_names <- sapply(c,function(x) {colnames(distM1)[x]})
  r_names <- sapply(r,function(x) {rownames(distM1)[x]})
  
  #print(c_names)
  #print(r_names)
  df = cbind(c_names,r_names)
  
  #sorts names in alphabetical order
  sort_str = function(x){
    if (x["c_names"]>x["r_names"]) return(as.vector(x))
    else {return(c(x["r_names"],x["c_names"]))}
  }
  
  sorted = data.frame(unique(t(apply(df,1,sort_str))),stringsAsFactors = FALSE)
  colnames(sorted) = c(1,2)
  #  values in matrices
  x = apply(sorted, 1, function(x){distM1[which(rownames(distM1) == x[1]),which(colnames(distM1) == x[2])]})
  y = apply(sorted, 1, function(x){distM2[which(rownames(distM1) == x[1]),which(colnames(distM1) == x[2])]})
  sorted = cbind(sorted,x,y)

  
  return(sorted)
  
}

#returns positions in matrix which values are between val1 and val2
find_dist_slice<-function(distM, val1, val2){
  b = which((distM>=val1)&(distM<=val2))
  #print(which((distM>val1)&(distM<val2),arr.ind=T))
  return(b)
}

#plots pairwise distance comparison plots for different pairs of regions 
#dna_object - list of DNA sequences (class DNAbin)
#step - step for dna regions/ window for regions is equal to step
#method - method of calculation distances ("pdist", "JC", "Kimura", "TN")
#fig_dir - path to directory for saving figures
#name_fig - prefix for figure name

#created directory fig_dir/step and saves pairwise distance comparison plot there

plot_dist = function(dna_object, step, method, fig_dir, name_fig){
  
  length_aln = length(dna_object[1,]) #length of alignment
  num_seq = length(dna_object[,1]) # number od sequences in alignment
  starts = seq(from=1, to = length_aln, by = step) #start positions of regions 
  ends = seq(from=step, to = length_aln, by = step) #end positions of regions
  
  
  if (length_aln%%step>0){ends=c(ends,length_aln)}
  
  #positions of genomes regions
  df_intervals = cbind(starts,ends)
  
  #list with distance matrices for each region
  dist_matrices = list()
  for (i in 1:nrow(df_intervals)){
    #slice of alignment
    slice = dna_object[1:num_seq, seq(from = df_intervals[i,"starts"], to = df_intervals[i,"ends"], by=1)]
    #dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "JC69")
    
    #calculation of distance matrix
    if (method == "pdist"){dist_matrices[[i]] = dist.gene(slice, method = "percentage",  pairwise.deletion = TRUE)}
    else {
      if (method == "JC"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "JC69")}
      if (method == "Kimura"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "K80")}
      if (method == "TN"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "K80")}
      #else{print("Unknown method")}
    }
    
  }
  
  #list with ggplots of distance comparison plots for different pairs of regions
  ggplots = list()
  k=1
  for(i in 1:(nrow(df_intervals)-1)){
    for(j in (i+1):nrow(df_intervals)){
      #xlab name
      rname1 = paste(toString(df_intervals[i,"starts"]),"-",toString(df_intervals[i,"ends"]),sep="")
      #ylab name
      rname2 = paste(toString(df_intervals[j,"starts"]),"-",toString(df_intervals[j,"ends"]),sep="")
      #dist1 = dist_matrices[[i]][lower.tri(dist_matrices[[i]],diag = FALSE)]+ rnorm(length(dist_matrices[[i]][lower.tri(dist_matrices[[i]],diag = FALSE)]), mean=0, sd = 0.0001)
      #dist2 = dist_matrices[[j]][lower.tri(dist_matrices[[j]],diag = FALSE)]+ rnorm(length(dist_matrices[[j]][lower.tri(dist_matrices[[j]],diag = FALSE)]), mean=0, sd = 0.0001)
      
      #adding random noise to values 
      if (method == "pdist"){
        dist1 = as.vector(dist_matrices[[i]]+ rnorm(length(dist_matrices[[i]]), mean=0, sd = 0.0001))
        dist2 = as.vector(dist_matrices[[j]] + rnorm(length(dist_matrices[[j]]), mean=0, sd = 0.0001))
      }
      else{
        dist1 = dist_matrices[[i]][lower.tri(dist_matrices[[i]],diag = FALSE)]+ rnorm(length(dist_matrices[[i]][lower.tri(dist_matrices[[i]],diag = FALSE)]), mean=0, sd = 0.0001)
        dist2 = dist_matrices[[j]][lower.tri(dist_matrices[[j]],diag = FALSE)]+ rnorm(length(dist_matrices[[j]][lower.tri(dist_matrices[[j]],diag = FALSE)]), mean=0, sd = 0.0001)
        
      }
      
      ggplots[[k]] = ggplot(data.frame(dist1,dist2),aes(dist1,dist2))+
        stat_bin2d(binwidth = 0.005)+scale_fill_gradientn(colours=c("blue","red"))+
        theme(legend.justification=c(1,0), legend.position=c(1,0))+
        xlab(rname1)+ylab(rname2)
      ggsave(ggplots[[k]],file= paste(fig_dir,name_fig,rname1,"_",rname2,".png",sep="") )
      xirk=k+1
      
    }
  }
}

# plots heatmap with RMSE in pairwise distance comparison plot for each pair of genomic regions
# dna_object -  list of DNA sequences (class DNAbin)
# step 
# window - length of genomic regions to compare
# method - method of calculation distances ("pdist", "JC", "Kimura", "TN")
# modification - pairwise deletion of positions with gaps or not

#returns matrix with rmse values for each pair f=of genomic regions

plot_rmse = function(dna_object, step,window, method, modification=NA){

  length_aln = length(dna_object[1,]) #length of alignment
  num_seq = length(dna_object[,1]) # number of sequences in alignment
  
  starts = seq(from=0, to =length_aln-window, by = step) # start positions of genomic regions
  starts[1]=1
  ends = seq(from=window, to = length_aln, by = step) # end positions of genomic regions
  if (length_aln%%step>step){ends=c(ends,length_aln)}
  
  df_intervals = cbind(starts,ends) #intervals
  
  #names = apply(df_intervals, 1, function(x){paste(toString(x[1]),toString(x[2]),sep="_")})
  
  #dataframe to store RMSE values of each comparison
  rmse_df = data.frame(matrix(ncol=length(starts), nrow = length(starts))) 
  colnames(rmse_df)=starts
  rownames(rmse_df)=starts
  
  #list of distance matrices for each pair of genomic regions
  dist_matrices = list()
  for (i in 1:nrow(df_intervals)){
    slice = dna_object[1:num_seq, seq(from = df_intervals[i,"starts"], to = df_intervals[i,"ends"], by=1)]
    
    #dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "JC69")
    if (method == "pdist"){
      if (modification=="pairwise"){
        dist_matrices[[i]] = dist.gene(slice, method = "percentage",  pairwise.deletion = TRUE)}
      else {      
        dist_matrices[[i]] = dist.gene(slice, method = "percentage",  pairwise.deletion = FALSE)}
    }
    
    else {
      if (method == "JC"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "JC69")}
      if (method == "Kimura"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "K80")}
      if (method == "TN"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "K80")}
      #else{print("Unknown method")}
    }
    
  }
  print(length(dist_matrices))
  for (i in 1:nrow(df_intervals)){
    for (j in 1:nrow(df_intervals)){
      #fits pairwise distance comparison plots linear model, calculates rmse
      rmse_df[i,j] =rmse(lm(dist_matrices[[j]]~dist_matrices[[i]]))
    }
  }
  colnames(rmse_df)
  return(rmse_df)
  
}


homedir = "D:\\kursach2\\8"
figdir = "D:\\kursach2\\8\\2500-3000"
setwd(homedir)

#genotype / prefix of input filename
gt="A1"

#alignment for genotype
alignment = read.dna(paste(gt,"fasta",sep="."), format="fasta")

#start and end positions of fragments of interest
st1 = 1
e1 = 260
st2 = 261
e2 =513

#list with pairwise distance comparison plot, distance matrices for regions st1-e1 and st2-e2
l =  plot_dist_test(alignment, st1, e1, st2, e2)
al_plot = l[[1]]+ggtitle(gt)
al_plot# + geom_smooth(method = 'lm', formula=y~x)
ggsave(filename= paste(figdir,gt,"-",toString(st1),"-",toString(e1),"_",toString(st2),"-",toString(e2),".png",sep=""), al_plot)

# pairwise distance comparison plot of odd and eve sites
contr = plot_control(alignment)
c_plot = contr[[1]] + ggtitle(paste(gt,"control",sep=" "))
c_plot 
ggsave(filename= paste(figdir,gt,"-control-",toString(st1),"-",toString(e1),"_",toString(st2),"-",toString(e2),".png",sep=""), c_plot)

#distance matrices for two regions
alDistM1 = as.matrix(l[[2]])
alDistM2 = as.matrix(l[[3]])

rec_names = find_recomb_names(alDistM1,0.05, 0.1, alDistM2 ,0.15, 0.25)
rec_names
write.table(rec_names, paste(figdir, gt, "-",toString(st1),"-",toString(e1),"_",toString(st2),"-",toString(e2),".txt",sep="") ,sep = "\t", quote=FALSE,row.names = FALSE,col.names = FALSE)

step = 100
name_fig = gt
#only "pdist","JC","Kimura", "TN"
method = "pdist"
#fig_dir = paste(figdir,gt,"\\step",toString(step),"\\",method,"\\",sep="")

if (dir.exists(paste(figdir,"\\", gt, sep="")) == FALSE){
  dir.create(paste(figdir,"\\", gt, sep=""))
} 

fig_dir = paste(figdir,"\\", gt, "\\step", toString(step), "\\", sep="")
if (dir.exists(fig_dir) == FALSE){
  dir.create(fig_dir)
}

print(step)
plot_dist(alignment, step, method, fig_dir, name_fig)


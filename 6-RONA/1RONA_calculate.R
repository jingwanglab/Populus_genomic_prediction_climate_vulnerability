#! /usr_storage/syp/.conda/envs/R/bin/Rscript --no-save --no-restore
.libPaths("/usr_storage/syp/.conda/envs/R/lib/R/library")
library(data.table)
library(plotrix)
REF_freq=read.table("/usr_storage/syp/pk_paper/new_RONA/1799_REF_freq.txt",head=T,row.names=1)
input_unLD_dir="/usr_storage/syp/pk_paper/new_RONA/step1_adaption-loci/unLD/"
now_env=read.csv("/usr_storage/syp/NC/RONA/NOW24.csv",head=T,row.names=1)
future_env_dir="/usr_storage/syp/NC/future_data/future_point_24"
outdir_SE="/usr_storage/syp/NC/RONA/SE"
outdir_weighted="/usr_storage/syp/NC/RONA/weighted"
outdir_weighted_SE="/usr_storage/syp/NC/RONA/weighted_SE"
now_env=now_env[,-c(1:2)]
setwd("/usr_storage/syp/NC/RONA")

#### RONA R function ####
rona_pred <- function(gen, present, future, significance, weighted){
  #present=as.data.frame(lapply(present,as.numeric))
  #future=as.data.frame(lapply(future,as.numeric))
  rona_output <- matrix(NA, nrow=nrow(gen), ncol=ncol(gen)) # empty matrix to save risk values
  allfreq_output <- matrix(NA, nrow=nrow(gen), ncol=ncol(gen)) # empty matrix to save future allele frequencies predicted
  reg_output <- matrix(NA, nrow=ncol(gen), ncol=5) # empty matrix to save regression parameters
  colnames(rona_output) <- colnames(gen) # track the locus id in the output
  for (j in 1:ncol(gen)) { # loop for each locus
    for (i in 1:nrow(gen)) { # loop for each population
      reg <- lm(gen[,j] ~  unlist(present)) # compute the linear regression
      rona_output[i,j] <- abs((coef(reg)[2] * future[i,] + coef(reg)[1]) - gen[i,j]) # calculate the risk of non-adaptedness
      allfreq_output[i,j] <- coef(reg)[2] * future[i,] + coef(reg)[1]
      reg_output[j,1] <- summary(reg)$fstatistic[1] # save parameters of the linear regression for each locus
      reg_output[j,2] <- summary(reg)$coefficients[2,4]
      reg_output[j,3] <- summary(reg)$r.squared
      reg_output[j,4] <- summary(reg)$adj.r.squared
      reg_output[j,5] <- summary(reg)$sigma
    }
  }
  SEM=matrix(NA, nrow=nrow(gen), ncol=1)
  for (i in 1:nrow(rona_output)) {
		SEM[i,]=std.error(rona_output[i,],na.rm)
	}
  rona_output <- as.data.frame(rona_output, row.names=rownames(gen)) # convert matrix to dataframe for downstream data maniputation and calculation
  colnames(allfreq_output) <- colnames(gen)
  rownames(allfreq_output) <- rownames(gen)
  colnames(reg_output) <- c("statistic", "pvalue", "Rsquared", "adjRsquared", "sigma") # rename output columns
  row.names(reg_output) <- colnames(gen) # rename output rows
  reg_output <- as.data.frame(reg_output) # convert matrix to dataframe for downstream data maniputation and calculation
  
  if (significance == "TRUE") {
    reg_output <- subset(reg_output, subset=reg_output$pvalue < 0.05) # significance set at 0.05
    rona_output <- rona_output[, row.names(reg_output)]
    allfreq_output <- allfreq_output[, row.names(reg_output)]
  } 
  
  avg_rona_output <- as.data.frame(rowMeans(rona_output, na.rm=TRUE)) # calculate the averaged unweighted RONA across the loci
  colnames(avg_rona_output) <- c("Avg_unweighted_RONA")
  
  if (weighted == "TRUE") {
    tmp <- matrix(NA, nrow=nrow(gen), ncol=1) # create a temporary matrix
    for (i in 1:nrow(gen)) {
      tmp[i] <- weighted.mean(rona_output[i,], reg_output[,"Rsquared"], na.rm=TRUE) # weighting the loci per R squared
    }
    avg_rona_output[,1] <- tmp
    colnames(avg_rona_output) <- c("Avg_weighted_RONA")
  }
  
  colnames(rona_output) <- paste("rona", colnames(rona_output), sep="_") # rename loci to avoid any confusion with input data
  colnames(allfreq_output) <- paste("pred", colnames(allfreq_output), sep="_")
  
  final_output <- list(reg_output, rona_output, avg_rona_output, allfreq_output,SEM) # create a list object for the final output
  return(final_output)
}

future_model=dir(path="/usr_storage/syp/NC/future_data/future_point_24",full.names = F)

for (model in future_model){
	future_env_file=list.files(path=paste(future_env_dir,model,sep="/"),full.names = F,pattern = "csv")
	for (i in 1:length(future_env_file)){
		a=unlist(strsplit(future_env_file[i],split="-SSP"))
	    outfile_weighted=paste(outdir_weighted,"/",a[1],"/",a[1],"_ssp",gsub("_24.csv","",a[2]),"NC_RONA_weighted.csv",sep="")
		outfile_SE=paste(outdir_SE,"/",a[1],"/",a[1],"_ssp",gsub("_24.csv","",a[2]),"NC_RONA_SE.csv",sep="")
		outfile_weighted_SE=paste(outdir_weighted_SE,"/",a[1],"/",a[1],"_ssp",gsub("_24.csv","",a[2]),"NC_RONA_weightedSE.csv",sep="")
		result_weighted=as.data.frame(matrix(nrow=24))
		result_weighted$POP=rownames(now_env)
		result_SE=as.data.frame(matrix(nrow=24))
		result_SE$POP=rownames(now_env)
		result_weightedSE=as.data.frame(matrix(nrow=24))
		result_weightedSE$POP=rownames(now_env)
		future_env=read.table(paste(future_env_dir,model,future_env_file[i],sep="/"),head=T,row.names=1)
		future_env=future_env[,-c(1:2)]
		for  (j in 1:19){
			BIO=paste("ssp",gsub("_24.csv","",a[2]),"_BIO",j,sep="")
			unLD_ID=read.table(paste(input_unLD_dir,"1779ID_BIO",j,".txt",sep=""),head=F,row.names=1)
			temp_BIO=REF_freq[,colnames(REF_freq)%in%rownames(unLD_ID)]
			temp_unLD=rona_pred(temp_BIO,now_env[j],future_env[j],"FALSE","TRUE")
			weight_rona=data.frame(matrix(unlist(temp_unLD[3]), nrow=24, byrow=T),stringsAsFactors=FALSE)
			se=data.frame(matrix(unlist(temp_unLD[5]), nrow=24, byrow=T),stringsAsFactors=FALSE)
			rona_se=data.frame(a=0)
			for (m in 1:nrow(weight_rona)){
				rona_se[m,]=paste(round(weight_rona[m,],4),"±(",round(se[m,],4),")",sep="")
			}
			result_weighted[,BIO]=weight_rona
			result_SE[,BIO]=se
			result_weightedSE[,paste(BIO,"(SE)",sep="")]=rona_se
		}
		write.csv(result_weighted[,-1],file=outfile_weighted,quote=F,row.names=F)
		write.csv(result_SE[,-1],file=outfile_SE,quote=F,row.names=F)
		write.csv(result_weightedSE[,-1],file=outfile_weighted_SE,quote=F,row.names=F)
	}
}

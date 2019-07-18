#!/bin/Rscript
args<-commandArgs(TRUE)
setwd(args[2]) 

sids.df = read.table(args[1],header=TRUE,sep="\t")
sids.df$Coordinates <- as.character(sids.df$Coordinates)
rownames(sids.df) = sids.df$Coordinates
sids.df = sids.df[,-1]

# store the column and row names
c.names <- names(sids.df)

# create a list of data frames (for storing output of the for loop)
tissue.df.list <- list()

# Iterate through each column in the matrix
for (i in c.names)
{
  # extract a tissue
  mytissue<-sids.df[,i, drop = FALSE]
  
  # Add a bin column (0 to place all non-expressed IDs in bin 0)
  mytissue$bin<-0
  
  # Identify the limits of the deciles
  expressed.ids<-mytissue[mytissue[[i]] != 0,, drop=F]
  decile.limits.df <- data.frame(quantile(expressed.ids[[i]], probs=seq(.1,.9,.1)))
  
  # Place all expressed IDs in bins 1 to 10
  mytissue[mytissue[[i]] > 0, c('bin')] <- 1
  mytissue[mytissue[[i]] > decile.limits.df['10%',], c('bin')] <- 2
  mytissue[mytissue[[i]] > decile.limits.df['20%',], c('bin')] <- 3
  mytissue[mytissue[[i]] > decile.limits.df['30%',], c('bin')] <- 4
  mytissue[mytissue[[i]] > decile.limits.df['40%',], c('bin')] <- 5
  mytissue[mytissue[[i]] > decile.limits.df['50%',], c('bin')] <- 6
  mytissue[mytissue[[i]] > decile.limits.df['60%',], c('bin')] <- 7
  mytissue[mytissue[[i]] > decile.limits.df['70%',], c('bin')] <- 8
  mytissue[mytissue[[i]] > decile.limits.df['80%',], c('bin')] <- 9
  mytissue[mytissue[[i]] > decile.limits.df['90%',], c('bin')] <- 10
  
  
  # add the my tissue data frame to the list of data frams
  tissue.df.list[[i]] <- mytissue
} 

# prepare outfiles
dir.create("Results")
setwd("Results")

for (i in seq(tissue.df.list))
{
  # store the dataframe in an easier variable
  df<-tissue.df.list[[i]] 
  # store the names of the columns
  columns<-names(df) 
  # remove the column with the qn values
  df[1]<-NULL 
  # rename the bin column to the tissue name
  colnames(df)[1]<-columns[1] 
  # name the dataframe
  assign(paste(columns[1]), df) 
  # write the dataframe to an outfile
  write.table(df, file = columns[1], sep = "\t", dec = ".", row.names = TRUE, col.names = NA, quote = FALSE) 
}



multimerge <- function(x)                                                                                    # Merges files in a specified folder
{
  filenames=list.files(path=x, full.names=TRUE)                                                              # lists the file names in a folder
  datalist = lapply(filenames, function(x){read.delim(file = x, header = TRUE, dec = ".", fill = FALSE)})    # reads the file contents into a list of data frames
  Reduce(function(x,y) {merge(x,y)}, datalist)                                                               # merges all of the data frames
}

# merge all dataframes to create tau profile
tau.in<-multimerge(getwd())
tau.in$X <- as.character(tau.in$X)                              # changes the IDs from factors to characters
row.names(tau.in) <- tau.in$X                                   # assigning the row names as the IDs
tau.in <- tau.in[,c(2:ncol(tau.in))]                      # removes 1st column

write.table(tau.in, file = "tau_inputFile.txt", sep = "\t", dec = ".", row.names = TRUE, col.names = NA, quote = FALSE)

cat ("Input file written for tau_score.pl")
#################################################################
# Apply TAU algorithm ###########################################
#################################################################

#  TAU Algorithm
#  T = (SUM (1-(sample intensity/max intensity)))/n-1

# Worked example: 
# For gene X expression profile ‘0 8 0 0 0 2 0 2 0 0 0 0’, Tau = 0.95 because:
#  ((9(1-(0/8)))	+	(1(1-(8/8)))	+	(2(1-(2/8))))/11

# For each row you need:
#  1)  n (the total number of columns)
#  2)  max intensity (the biggest value in the row)

###############################################################
####### This part has been written in perl for speed ##########
####### tau_score.pl ##########################################
###############################################################

# # Create a dataframe to store tau value and tissue of max expression
# tau.file<-data.frame(tau = 0,maxTissue = names(tau.in)[apply(tau.in,1,which.max)])

# # get row names of tau.in and apply them to tau.file
# r.names <- row.names(tau.in)
# row.names(tau.file)<-r.names

# # Calculate a tau value for every row
# n <-length(tau.in)
# for (x in r.names)
# {
#   max.intensity <- max(tau.in[x,])                   # gets the max expression of the gene
#   max.tissue <- names(which.max(tau.in[x,]))         # gets the tissue in which max expression occurs
#   sum <-0                                                
#   for(y in c.names)
#   {
#     sample.intensity <- tau.in[x,y]                  # gets the expression of the gene in the tissue
#     norm <- 1-(sample.intensity/max.intensity)            # adjusts for max expression
#     sum <- sum+norm                                       
#   }
#   tau <- sum/(n-1)                                        # adjusts for number of tisses
#   tau.file[x,'tau'] = tau                               # populates the tau column
# }


# write.table(tau.file, file = "tau_results.txt", sep = "\t", dec = ".", row.names = TRUE, col.names = NA, quote = FALSE)
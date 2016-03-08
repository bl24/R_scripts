###### 2016-01-25
###### Script for combining multiple Excel files. Intended for use on VarSeq output.
###### This script assumes that VarSeq reports were exported as text files.

# Set the working directory (where files are saved)
setwd() #Insert directory here

# Get all the right file names
file_names = list.files(getwd())
file_names = file_names[grepl(".txt",file_names)]

# Print file_names vector
file_names

# get the read.csv function working
files = read.delim("", header=F, stringsAsFactors = T) #Insert name of file in the brackets

# use lapply to apply the read function to all values of file_names
files = lapply(file_names, read.delim, header=F, stringsAsFactors = F)
files = do.call(rbind,files)

# create new file using output
write.csv(files, file = "merged.csv", na = "") #Customize file name when appropriate

###### To remove duplicate rows:

# read in merged file
merged <- read.csv("merged.csv", header = T)

# remove duplicates
dedup <- merged[!duplicated(merged), ]

# create file for deduped data
write.csv(dedup, file = "merged_dedup.csv", na = "")
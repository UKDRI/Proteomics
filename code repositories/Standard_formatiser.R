#
# Libraries
#

library(odbc)
library(DBI)
library(RPostgres)
library(plyr)
# 
# WD check
#

print(getwd())

#
# Determines what type of search was used
#

if(identical(character(0),list.files(pattern = "proteinGroups.txt"))==FALSE){
 # convert to standard format 
  print("MaxQuant")
  combined_protein <- read.delim("proteinGroups.txt")
  temp_format <- combined_protein
  temp_format <- temp_format[,c(2,8,8,10)]
  colnames(temp_format) <- c("ID","Uniprot_Entry","Species","Total_Peptides")
  temp_format$Species <- gsub(".*_","",gsub(" .*","",gsub(".*@","",gsub("|","@",temp_format$Species, fixed =T))))
  temp_format$Uniprot_Entry <- gsub("_.*","",gsub(".*@","",gsub("|","@",temp_format$Uniprot_Entry, fixed =T)))
  temp_format$Type <- "MaxQuantDDA"
  temp_format$SearchDate <- strsplit(as.character(file.info("proteinGroups.txt")[5][1,1]),split = " ")[[1]][1]
  temp_format$Xcode <- strsplit(getwd(), split = "/")[[1]][grep("X",strsplit(getwd(), split = "/")[[1]])]
  temp_format$prikey <- paste(temp_format$ID,temp_format$Xcode,temp_format$SearchDate,sep="")
  temp_format <- temp_format[,c(ncol(temp_format),1:(ncol(temp_format)-1))]
  #align column names with DB
  colnames(temp_format) <- c("prikey","id","uniprot","species","peptides","type","date","xcode")
  standard.format <- temp_format
  
}else if(identical(character(0),list.files(pattern = "combined_protein.tsv"))==FALSE){
  
  print("FragPipe")
  
  # import file
  
  combined_protein <- read.delim("combined_protein.tsv")
  
  # Determine whether DDA or DIA
  
  if((sum(combined_protein[,grep("Intensity",colnames(combined_protein))[1]]) == 0) == TRUE){
    #DDA
    print("DDA")
    
    temp_format <- combined_protein
    temp_format <- temp_format[,c(1,3,6,12)]
    
    colnames(temp_format) <- c("ID","Uniprot_Entry","Species","Total_Peptides")
    
    temp_format$Type <- "FragPipeDDA"
    
    temp_format$SearchDate <- strsplit(as.character(file.info("combined_protein.tsv")[5][1,1]),split = " ")[[1]][1]
    
    temp_format$Xcode <- strsplit(getwd(), split = "/")[[1]][grep("X",strsplit(getwd(), split = "/")[[1]])]
    temp_format$prikey <- paste(temp_format$ID,temp_format$Xcode,temp_format$SearchDate,sep="")
    temp_format <- temp_format[,c(ncol(temp_format),1:(ncol(temp_format)-1))]
    colnames(temp_format) <- c("prikey","id","uniprot","species","peptides","type","date","xcode")
    standard.format <- temp_format
    
  }else{
    #DIA
    print("DIA")
    
    temp_format <- combined_protein
    temp_format <- temp_format[,c(1,3,6,12)]
    
    colnames(temp_format) <- c("ID","Uniprot_Entry","Species","Total_Peptides")
    
    temp_format$Type <- "FragPipeDDALFQ"
    
    temp_format$SearchDate <- strsplit(as.character(file.info("combined_protein.tsv")[5][1,1]),split = " ")[[1]][1]
    
    temp_format$Xcode <- strsplit(getwd(), split = "/")[[1]][grep("X",strsplit(getwd(), split = "/")[[1]])]
    temp_format$prikey <- paste(temp_format$ID,temp_format$Xcode,temp_format$SearchDate,sep="")
    temp_format <- temp_format[,c(ncol(temp_format),1:(ncol(temp_format)-1))]
    colnames(temp_format) <- c("prikey","id","uniprot","species","peptides","type","date","xcode")
    standard.format <- temp_format
  }
  
  
  # convert to standard format 
  
  
  
  
}else if(identical(character(0),list.files(pattern = "report.tsv"))==FALSE){
  # convert to standard format
  print("Dia-NN")
  
  # import file
  
  combined_protein <- read.delim("report.tsv")
  temp_format <- combined_protein
  temp_format <- temp_format[,c(3,5,5)]
  temp_format$Tot <- NA
  colnames(temp_format) <- c("ID","Uniprot_Entry","Species","Total_Peptides")
  temp_format$Species <- gsub(".*_","",gsub(" .*","",gsub(".*@","",gsub("|","@",temp_format$Species, fixed =T))))
  temp_format$Uniprot_Entry <- gsub("_.*","",gsub(".*@","",gsub("|","@",temp_format$Uniprot_Entry, fixed =T)))
  temp_format$Type <- "DiaNNDIA"
  temp_format$SearchDate <- strsplit(as.character(file.info("report.tsv")[5][1,1]),split = " ")[[1]][1]
  
  if(identical(character(0),strsplit(getwd(), split = "/")[[1]][grep("X",strsplit(getwd(), split = "/")[[1]])])==FALSE){
    temp_format$Xcode <- strsplit(getwd(), split = "/")[[1]][grep("X",strsplit(getwd(), split = "/")[[1]])]
    
  }else{
    
    temp_format$Xcode <- paste(strsplit(getwd(), split = "/")[[1]][4],strsplit(getwd(), split = "/")[[1]][5],sep = "")
  }
  
  
  
  temp_format <- temp_format[order(temp_format$Uniprot_Entry),]
  
  temp_format$Total_Peptides <- count(temp_format$Uniprot_Entry)[match(temp_format$Uniprot_Entry,count(temp_format$Uniprot_Entry)[,1]),2]
  
  temp_format <- unique(temp_format)
  
  temp_format$prikey <- paste(temp_format$ID,temp_format$Xcode,temp_format$SearchDate,sep="")
  temp_format <- temp_format[,c(ncol(temp_format),1:(ncol(temp_format)-1))]
  #align column names with DB
  colnames(temp_format) <- c("prikey","id","uniprot","species","peptides","type","date","xcode")
  standard.format <- temp_format
  
  
}else if(identical(character(0),list.files(pattern = "spectro.tsv"))==FALSE){
  # convert to standard format
  print("spectronaut")
}else{
  message("Unsupported input files / input files missing")
}
  
#
# Print out standard format
#

exp_name <- strsplit(dirname(getwd()), split = "/")[[1]][4]
sub_name <- strsplit(getwd(), split = "/")[[1]][(length(strsplit(getwd(), split = "/")[[1]]))]

write.csv(standard.format ,paste(exp_name,"_",sub_name,".csv",sep = ""))

# Connect to the default postgres database

con <- dbConnect(Postgres(), user = "postgres", password = "119583rescue")

dbAppendTable(con,"xresultsdb",standard.format)











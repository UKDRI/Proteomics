library(readxl)
library(readr)
library(reshape2)
library(ggplot2)
library(limma)
library(pheatmap)
library(mice)
library(preprocessCore)
library(corrplot)
set.seed(500)
library(odbc)
library(DBI)
library(RPostgres)
library(dplyr)
library(Rtsne)
library(rrcovNA)
beth.colours.reverse <- colorRampPalette( c("yellow","darkgrey","blue"), space="rgb")(100)
beth.colours.colour <- colorRampPalette( c("#D52D00","#EF7627","#FF9A56","white","#D162A4","#B55690","#A30262"), space="rgb")(100)
source("~/Scripts/BETH_functions.r")

# Upload of DIA-NN output file

Aonly <- read_delim("report.pg_matrix.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

Aonly <- Aonly[grep(";",Aonly$Protein.Group,invert = T),]

# Removing excessive filenaming

colnames(Aonly) <- gsub(".*DRI","DRI",colnames(Aonly))
colnames(Aonly) <- gsub("_S2.*","",colnames(Aonly))
colnames(Aonly) <- gsub(".raw","",colnames(Aonly))
colnames(Aonly) <- gsub(".2024.*","",colnames(Aonly))
colnames(Aonly) <- gsub("_202507.*","",colnames(Aonly))
# Removing excessive columns

# removing plasma and headers

m.Aonly <- data.matrix(Aonly[,-c(1:6,61,62)])

# Fixing column naming

colnames(m.Aonly) <- gsub("DRI_","DRI",colnames(m.Aonly))

# Adding condition info now

#Condition <- mapper$X10[match(colnames(m.Aonly),mapper$X1)]
Condition <- mapper$X10[match(colnames(m.Aonly),mapper$X1)]
sample.types <- mapper$X11[match(colnames(m.Aonly),mapper$X1)]
types <- unique(mapper$X11)

unique(Condition)



#m.Aonly <- m.Aonly[,c(1:6)]
rownames(m.Aonly) <- gsub(";.*","",Aonly$Protein.Group)

saved.con <- Condition
saved.mAonly <- m.Aonly


mainDir <- getwd()




for(k in 1:length(unique(types))){
  
  m.Aonly <- saved.mAonly
  Condition <- saved.con
  subdiv.con <- Condition[grep(unique(types)[k],fixed=T,invert=F,sample.types)]

  
  
for(l in 1:1){

setwd(mainDir)

m.Aonly <- saved.mAonly
Condition <- saved.con
    
    
m.Aonly <- m.Aonly[,grep(paste(unique(sample.types)[k],unique(subdiv.con)[l],sep="|"), invert=F, sample.types)]
Condition <- Condition[grep(paste(unique(sample.types)[k],unique(subdiv.con)[l],sep="|"), invert=F,sample.types)]


SubDir = paste(paste(unique(Condition),collapse  = "_vs_"))
dir.create(file.path(mainDir, SubDir))
setwd(file.path(mainDir, SubDir))    
    
m.Aonly <- m.Aonly[rowSums(is.na(m.Aonly)) != ncol(m.Aonly),]


#m.Aonly <- m.Aonly[,which(Condition == unique(saved.con[i]),saved.con[j])



#m.Aonly <- which(mapper$X5[match(colnames(m.Aonly),mapper$X1)]=="APP23")

#
# QC plots
#


X20220915 <- m.Aonly

# Plots to show metrics

na_count <-sapply(X20220915, function(y) sum(length(which(is.na(y)))))

na.df <- as.data.frame(nrow(X20220915) - na_count)
colnames(na.df) <- "ProteinIDs"
na.df$Sample <- rownames(na.df)

# Melt for ggplot purposes

me.data <- melt(X20220915)
me.data <- me.data[complete.cases(me.data),]
me.data$Condition <- gsub(".","",fixed = T, gsub("[0-9]","",me.data$Var2))



pdf("proteins_per_fraction.pdf")

ggplot(me.data,aes(Var2, fill = Condition))+
  geom_bar()+
  #scale_fill_manual(values = c('#ffe119', '#4363d8', '#f58231', '#dcbeff', '#800000', '#000075', '#a9a9a9',  '#000000'))+
  ylab("Proteins per fraction")+
  ggtitle("Amount of proteins in each sample")+
  xlab("Fraction name")+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

rownas <- as.data.frame(rowSums(is.na(X20220915)))
rownas$Protein.Group <- rownames(X20220915)
colnames(rownas)[1] <- "Na.Count"

rownas$IDin <- (ncol(X20220915) - rownas$Na.Count)

pdf("ID_overlap.pdf")

ggplot(rownas, aes(IDin))+
  geom_bar(fill = "#ffe119")+
  ylab("Number of proteins")+
  xlab("Identified in number of samples")+
  ggtitle("Protein identifications overlap")+
  theme_bw()

dev.off()

me.data <- melt(X20220915)

rownas$Na.Count[match(me.data$Var1,rownas$Protein.Group)]>0

me.data$missing <- rownas$Na.Count[match(me.data$Var1,rownas$Protein.Group)]>1

me.data <- me.data[complete.cases(me.data),]
me.data$value <- log(me.data$value,2)

pdf("Missing_values_density.pdf")

ggplot(me.data,aes(value,col = missing))+
  geom_density()+
  scale_color_manual(values = c('#4363d8', '#f58231', '#dcbeff', '#800000', '#000075', '#a9a9a9',  '#000000'))+
  xlab("Log2 Intensity")+
  ylab("Density")+
  theme_bw()

dev.off()

# Density plot of samples with and without missing values

me.data <- melt(X20220915)

rownas$Na.Count[match(me.data$Var1,rownas$Protein.Group)]>0

me.data$missing <- rownas$Na.Count[match(me.data$Var1,rownas$Protein.Group)]>0

me.data <- me.data[complete.cases(me.data),]
me.data$value <- log(me.data$value,2)

pdf("density_plots.pdf")

ggplot(me.data,aes(value,col = Var2))+
  geom_density(show.legend = F)+
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  xlab("Log2 Intensity")+
  facet_wrap(vars(Var2))+
  ylab("Density")+
  theme_bw()

dev.off()

# Normalisation


library(preprocessCore)


norm.dat <- as.data.frame(normalizeMedianValues(X20220915))

norm.dat$Protein <- rownames(X20220915)
colnames(norm.dat) <- c(colnames(X20220915),"Protein.Group")

me.nor.dat <- melt(norm.dat)
me.data <- melt(X20220915)

me.nor.dat$Norm <- "Normalised"
me.data$Norm <- "Not Normalised"

colnames(me.data) <- colnames(me.nor.dat)

comb.me.data <- rbind(me.nor.dat,me.data)
comb.me.data$value <- log(comb.me.data$value,2)


pdf("Normalisation_plot.pdf")

ggplot(comb.me.data,aes(value,variable, fill = Norm))+
  facet_grid(cols=vars(Norm))+
  scale_fill_manual(values = c('#4363d8', '#f58231', '#dcbeff', '#800000', '#000075', '#a9a9a9',  '#000000'))+
  xlab("")+
  ylab("")+
  geom_boxplot()+
  theme_bw()

dev.off()


#
# Median normalisation
#

name.save <- colnames(m.Aonly)
rowname.save <- rownames(m.Aonly)

m.Aonly <- normalizeMedianValues(m.Aonly)

rownames(m.Aonly) <- rowname.save
colnames(m.Aonly) <- name.save

#
# Imputation
#


m.aonly.save <- m.Aonly


df_data =m.Aonly

library(rrcovNA)


# This will work out proteins with enough coverage that will allow for decent imputation

zeroProp <- 0.25
nSample=ncol(df_data)

vect_zeroProp<-apply(df_data,1,(function(col) length(col[col==0])/nSample))
index = which(as.numeric(vect_zeroProp) <= zeroProp)
df_data2<-m.Aonly[which(as.numeric(vect_zeroProp) <= zeroProp),]
N=nrow(df_data2)
print(paste(N," proteins with more than ",(1 -zeroProp)*100,"% of ind with raw exp > 0"))
print(nrow(df_data2))


c.mmao <- impSeqRob(df_data2)$x


# This returns the proteins that did not have enough coverage to the imputed dataset

seq2 <- 1:nrow(m.Aonly)
comb.AO <- rbind(c.mmao,m.Aonly[seq2[!seq2 %in% match(rownames(c.mmao),rownames(m.Aonly))],])

#comb.AO <- df_data2


#
# The experimental condition !! Change this depending on which experimental groups you wish to showcase !!
#



#Condition <- gsub(".*_","",colnames(comb.AO))
#Condition <- mapper$X5[match(colnames(comb.AO),mapper$X1)]
#Condition <- gsub(" ","_",Condition)

#Condition <- c(rep("Sham",3),rep("AAV",12))
#Condition <- c(rep("APP23_young",3),rep("Neg_young",3))

# Removing any stray NAs as 1 (not the best solution but a necessary evil)


comb.AO[comb.AO<0] <- 1
comb.AO[is.na(comb.AO)] <- 1

new_condition <- Condition

con.combao <- comb.AO
colnames(con.combao) <- Condition


# Correlation plaots
cor_pdf_name <- paste(gsub("_.*","",strsplit(getwd(), split = "/")[[1]][length(strsplit(getwd(), split = "/")[[1]])]),"cor_plots.pdf",sep = "_")

pdf(cor_pdf_name)
corrplot(cor(comb.AO), method = 'color', order = 'AOE')
corrplot(cor(con.combao), method = 'color', order = 'AOE')
dev.off()


df_data2 <- comb.AO

write.csv(df_data2, file="processed_normalised_filtered.csv",row.names=TRUE, quote=FALSE)

m.Aonly <- df_data2

colnames(m.Aonly) <- gsub("DRIX_099_","",gsub("_S2.*","",colnames(m.Aonly)))

colnames(m.Aonly) <- paste(Condition,colnames(m.Aonly),sep=".")

csv_out_name <- paste(gsub("_.*","",strsplit(getwd(), split = "/")[[1]][length(strsplit(getwd(), split = "/")[[1]])]),"_raw.csv",sep = "_")

write.csv(m.Aonly,csv_out_name)

m.Aonly <- log(m.Aonly,2)



#
# Statistical Comparisons, heat maps, PCA, tSNE, output for volcanos / Curtain
#

m.Aonly <- m.Aonly[,order(Condition)]
Condition <- Condition[order(Condition)]

m.save <- m.Aonly
con.save <- Condition



# Looping statistical comparisons

m.Aonly <- m.save
Condition <- con.save

zeroProps <- c(0,0.5,1)
pros <- c("100","50","0")

tryCatch(
  expr = {

for(b in 1:3){
  
  
  for(i in 1:1){
    
    m.Aonly <- m.save
    Condition <- con.save
    sub.con <- Condition[grep(unique(Condition)[i],fixed=T,invert=T,Condition)]
    
    for(j in 1:i){
      
      message(paste("Running cutoff:", pros[b]," for ",unique(Condition)[i]))
      
      m.Aonly <- m.save
      Condition <- con.save
      
      
      #m.Aonly <- m.Aonly[,grep(paste(unique(Condition)[i],unique(sub.con)[j],sep="|"), invert=F, Condition)]
     # Condition <- Condition[grep(paste(unique(Condition)[i],unique(sub.con)[j],sep="|"), invert=F,Condition)]
      
      
      
      m.Aonly[m.Aonly==0] <- NA
      m.Aonly[m.Aonly=="-Inf"] <- NA
      #m.Aonly <- m.Aonly[complete.cases(m.Aonly),]
      
      zeroProp <- zeroProps[b]
      nSample=ncol(m.Aonly)
      
      vect_zeroProp<-apply(m.Aonly,1,(function(col) length(col[col==0])/nSample))
      index = which(as.numeric(vect_zeroProp) <= zeroProp)
      m.Aonly<-m.Aonly[which(as.numeric(vect_zeroProp) <= zeroProp),]
      N=nrow(m.Aonly)
      print(paste(N," proteins with more than ",(1 -zeroProp)*100,"% of ind with raw exp > 0"))
      print(nrow(m.Aonly))
      
      
      new_condition <- Condition
      
      
      df_data2 = m.Aonly
      
      df_data2[is.na(df_data2)] <- 1
      df_data2[which(df_data2 == "-Inf")] <- 1
      df_data2 <- df_data2[!(rowSums(df_data2) == ncol(df_data2)),]
      
      fac.con <- factor(Condition)
      t.df_d2 <- t(df_data2)
      
      # anova getting p.values

      anova_results <- apply(df_data2, 1, function(protein_values) {
        fit <- aov(protein_values ~ fac.con)
        summary(fit)[[1]][["Pr(>F)"]][1]
      })
      
      aov.out <- as.data.frame(anova_results)
      
      # post hoc tukey HSD
     
      
      
      anova_posthoc_results <- lapply(rownames(df_data2), function(prot) {
        values <- as.numeric(df_data2[prot, ])
        
        # Skip if all values are NA or constant
        if (all(is.na(values)) || var(values, na.rm = TRUE) == 0) {
          return(list(Protein = prot, ANOVA_P = NA, PostHoc = NULL))
        }
        
        fit <- try(aov(values ~ fac.con), silent = TRUE)
        if (inherits(fit, "try-error")) {
          return(list(Protein = prot, ANOVA_P = NA, PostHoc = NULL))
        }
        
        smry <- summary(fit)[[1]]
        anova_p <- smry[["Pr(>F)"]][1]
        
        # Run Tukey HSD regardless (you can change to anova_p < 0.05 if you want to limit)
        tukey <- try(TukeyHSD(fit), silent = TRUE)
        if (inherits(tukey, "try-error")) {
          tukey_df <- NULL
        } else {
          tukey_df <- as.data.frame(tukey$fac.con)
          tukey_df$Contrast <- rownames(tukey_df)
          tukey_df$Protein <- prot
        }
        
        return(list(Protein = prot, ANOVA_P = anova_p, PostHoc = tukey_df))
      })
      
      # ANOVA summary
      anova_summary <- data.frame(
        Protein = sapply(anova_posthoc_results, `[[`, "Protein"),
        ANOVA_P = sapply(anova_posthoc_results, `[[`, "ANOVA_P")
      )
      
      # Adjust for multiple testing (FDR)
      anova_summary$Adj_P <- p.adjust(anova_summary$ANOVA_P, method = "fdr")
      
      # Combine Tukey results
      posthoc_df <- do.call(rbind, lapply(anova_posthoc_results, function(x) x$PostHoc))
      
      write.csv(anova_summary, paste(pros[b],i,j,"ANOVA_summary.csv",sep="_"), row.names = FALSE)
      write.csv(posthoc_df, paste(pros[b],i,j,"Tukey_posthoc_results.csv",sep = "_"), row.names = FALSE)
      
      beth.colours.reverse <- colorRampPalette( c("yellow","darkgrey","blue"), space="rgb")(100)
      
      
      
      pdf(paste(pros[b],i,j,"ANOVA_comparison_PCAs_heatmap.pdf",sep="_"))
      
      print(ggbiplot(prcomp(t(df_data2), scale=F), scale = 1, obs.scale = 0, var.scale = 0, var.axes = FALSE,  varname.size = 0.75, ellipse.prob = 0.45, groups = Condition,  ellipse = TRUE, circle = FALSE))
      
      print(ggbiplot(prcomp(t(df_data2), scale=F), scale = 1, obs.scale = 0, var.scale = 0,labels = colnames(m.Aonly), labels.size = 2, var.axes = FALSE,  varname.size = 0.75, ellipse.prob = 0.45, groups = Condition,  ellipse = TRUE, circle = FALSE))
      
      print(ggbirepelplot(prcomp(t(df_data2), scale=F), scale = 1, obs.scale = 0, var.scale = 0, var.axes = FALSE,  varname.size = 0.75, ellipse.prob = 0.45, groups = Condition,  ellipse = TRUE, circle = FALSE))
      
      print(ggbirepelplot(prcomp(t(df_data2), scale=F), scale = 1, obs.scale = 0, var.scale = 0,labels = colnames(m.Aonly), labels.size = 2, var.axes = FALSE,  varname.size = 0.75, ellipse.prob = 0.45, groups = Condition,  ellipse = TRUE, circle = FALSE))
      
      mat_data = df_data2
      index = which(mat_data == "-Inf")
      mat_data[index] = 1
      patient_group = Condition
      Experiment_type= factor(patient_group)
      annotation_col = data.frame(Experiment_type)
      rownames(annotation_col)=colnames(mat_data)
      colnames(annotation_col) <- " "
      mat_data[]
      # Specify colors
      ann_colors = list(" " = c("#999999","#E69F00","#56B4E9","#009E73"))#,"#001100"))
      names(ann_colors[[1]]) <- unique(Condition)
      
      mat_data <- mat_data[which(rowSums(mat_data)!=ncol(mat_data)),]
      
      pheatmap(mat_data, annotation_colors = ann_colors, color = beth.colours.reverse,show_colnames = T, show_rownames = F,clustering_distance_rows ="euclidean", clustering_distance_cols ="euclidean",  annotation_col = annotation_col,main = paste0(" "),scale = "row")
      
      
      dev.off()
      
      
    }
    
    
    
  }
}

},
error = function(e){
  message('Error!')
  print(e)
}
)


}
}

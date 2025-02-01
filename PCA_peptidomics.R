# R code used to generate supplementary figures S1-S4 in Bera et al (PLOS One under revision) 
# and to prepare data input files for plogo visualisations presented in figures 2,3 and 4 
# This is an example of the code used to analyse garbanzo soaked non tryptic output file from MaxQuant excluding peptides below 0.95 PEP (see example input file 'peptides.txt' provided with this code). 
# other input files (see PRIDE:PXD048224) can be analysed by changing the input file here. 

setwd("C:/Users/indrani/Desktop/seqlogo&PCA/soaked")
getwd()

#....
library(readr)
library(dplyr)
library(stringr)
library(fastDummies)
library(plsgenomics)
library(ggplot2)
library(stringr)
library(FactoMineR)
library(factoextra)
library(factoextra)
library(ggfortify)




# Read the text file into a data frame
peptides <- read.table("peptides.txt", header = TRUE, sep = "\t")

# Sort the data frame based on the 'PEP' column
peptides <- peptides[order(peptides$PEP), ]

# Filter the data frame to include only rows where the 'PEP' score is above 0.95
peptides_sorted <- peptides[peptides$PEP > 0.95, ]

# Print the filtered data frame
print(peptides_sorted)

# Read the csv file
peptides_sorted = read.csv("peptides_soaked_sorted_pep.csv", header=TRUE)

# select sequences and nterminal window
peptides_nterminal=peptides_sorted[ , c(1, 2)]
peptides_cterminal=peptides_sorted[ , c(1, 3)]

# write a new csv file
write.csv(peptides_nterminal, file = 'soaked_nterminal.csv')
write.csv(peptides_cterminal, file = 'soaked_cterminal.csv')

# select only n-terminal and c-terminal window
nterm_cleav=data.frame(peptides_nterminal[,c(2)])
write.table(nterm_cleav, file = 'nterm_soaked.csv',append = FALSE,sep = ",", eol = "\n", na = "NA", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))

cterm_cleav=data.frame(peptides_cterminal[,c(2)])
write.table(cterm_cleav, file = 'cterm_soaked.csv',append = FALSE,sep = ",", eol = "\n", na = "NA", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))

# select subset of the string
nd <- read.csv("nterm_soaked.csv", header=TRUE, colClasses="character")
nd1=data.frame(nd)

library(stringr)
nd1$substring = substring(nd1$peptides_nterminal...c.2..,12,19)
nd1
write.table(nd1, file = 'soaked_nterm.csv',append = FALSE,sep = ",", eol = "\n", na = "NA", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))

#---------
cd <- read.csv("cterm_soaked.csv", header=TRUE, colClasses="character")
cd1=data.frame(cd)

cd1$substring = substring(cd1$peptides_cterminal...c.2..,12,19)
cd1
write.table(cd1, file = 'soaked_cterm.csv',append = FALSE,sep = ",", eol = "\n", na = "NA", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))

nterminal = read.csv(file = "soaked_nterm.csv")

cterminal = read.csv(file = "soaked_cterm.csv")

nt <- filter(nterminal, str_detect(nterminal$substring, "\\*", negate=TRUE))
ct <- filter(cterminal, str_detect(cterminal$substring, "\\*", negate=TRUE))

y <- nt %>% group_by(substring)
y1 <- ct %>% group_by(substring)

#first step is expand out the peptide sequence into seven columns
y <- y %>%mutate(P4= str_sub(substring,1,1) )
y <- y %>%mutate(P3= str_sub(substring,2, 2) )
y <- y %>%mutate(P2= str_sub(substring,3, 3) )
y <- y %>%mutate(P1= str_sub(substring,4, 4) )
y <- y %>%mutate(p1= str_sub(substring,5, 5) )
y <- y %>%mutate(p2= str_sub(substring,6, 6) )
y <- y %>%mutate(p3= str_sub(substring,7, 7) )
y <- y %>%mutate(p4= str_sub(substring,8, 8) )

y1 <- y1 %>%mutate(P4= str_sub(substring,1,1) )
y1 <- y1 %>%mutate(P3= str_sub(substring,2, 2) )
y1 <- y1 %>%mutate(P2= str_sub(substring,3, 3) )
y1 <- y1 %>%mutate(P1= str_sub(substring,4, 4) )
y1 <- y1 %>%mutate(p1= str_sub(substring,5, 5) )
y1 <- y1 %>%mutate(p2= str_sub(substring,6, 6) )
y1 <- y1 %>%mutate(p3= str_sub(substring,7, 7) )
y1 <- y1 %>%mutate(p4= str_sub(substring,8, 8) )


# then expand the seven columns into twenty dummies
y<- dummy_cols(y, select_columns = 'P4')
y<- dummy_cols(y, select_columns = 'P3')
y<- dummy_cols(y, select_columns = 'P2')
y<- dummy_cols(y, select_columns = 'P1')
y<- dummy_cols(y, select_columns = 'p1')
y<- dummy_cols(y, select_columns = 'p2')
y<- dummy_cols(y, select_columns = 'p3')
y<- dummy_cols(y, select_columns = 'p4')

y1<- dummy_cols(y1, select_columns = 'P4')
y1<- dummy_cols(y1, select_columns = 'P3')
y1<- dummy_cols(y1, select_columns = 'P2')
y1<- dummy_cols(y1, select_columns = 'P1')
y1<- dummy_cols(y1, select_columns = 'p1')
y1<- dummy_cols(y1, select_columns = 'p2')
y1<- dummy_cols(y1, select_columns = 'p3')
y1<- dummy_cols(y1, select_columns = 'p4')


write.table(y,file='n_soaked.csv',append = FALSE,sep = ",", eol = "\n", na = "NA", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))

# add in the amino acid composition as additional predictors
y <- mutate(y,countA = str_count(y$substring,"A"))
y <- mutate(y,countC = str_count(y$substring,"C"))
y <- mutate(y,countD = str_count(y$substring,"D"))
y <- mutate(y,countE = str_count(y$substring,"E"))
y <- mutate(y,countF = str_count(y$substring,"F"))
y <- mutate(y,countG = str_count(y$substring,"G"))
y <- mutate(y,countH = str_count(y$substring,"H"))
y <- mutate(y,countI = str_count(y$substring,"I"))
y <- mutate(y,countK = str_count(y$substring,"K"))
y <- mutate(y,countL = str_count(y$substring,"L"))
y <- mutate(y,countM = str_count(y$substring,"M"))
y <- mutate(y,countN = str_count(y$substring,"N"))
y <- mutate(y,countP = str_count(y$substring,"P"))
y <- mutate(y,countQ = str_count(y$substring,"Q"))
y <- mutate(y,countR = str_count(y$substring,"R"))
y <- mutate(y,countS = str_count(y$substring,"S"))
y <- mutate(y,countT = str_count(y$substring,"T"))
y <- mutate(y,countV = str_count(y$substring,"V"))
y <- mutate(y,countW = str_count(y$substring,"W"))
y <- mutate(y,countY = str_count(y$substring,"Y"))

y1 <- mutate(y1,countA = str_count(y$substring,"A"))
y1 <- mutate(y1,countC = str_count(y$substring,"C"))
y1 <- mutate(y1,countD = str_count(y$substring,"D"))
y1 <- mutate(y1,countE = str_count(y$substring,"E"))
y1 <- mutate(y1,countF = str_count(y$substring,"F"))
y1 <- mutate(y1,countG = str_count(y$substring,"G"))
y1 <- mutate(y1,countH = str_count(y$substring,"H"))
y1 <- mutate(y1,countI = str_count(y$substring,"I"))
y1 <- mutate(y1,countK = str_count(y$substring,"K"))
y1 <- mutate(y1,countL = str_count(y$substring,"L"))
y1 <- mutate(y1,countM = str_count(y$substring,"M"))
y1 <- mutate(y1,countN = str_count(y$substring,"N"))
y1 <- mutate(y1,countP = str_count(y$substring,"P"))
y1 <- mutate(y1,countQ = str_count(y$substring,"Q"))
y1 <- mutate(y1,countR = str_count(y$substring,"R"))
y1 <- mutate(y1,countS = str_count(y$substring,"S"))
y1 <- mutate(y1,countT = str_count(y$substring,"T"))
y1 <- mutate(y1,countV = str_count(y$substring,"V"))
y1 <- mutate(y1,countW = str_count(y$substring,"W"))
y1 <- mutate(y1,countY = str_count(y$substring,"Y"))

# Creating input files for plogo analysis , in capitals and refer to figures

write.table(y,file='n_soaked_new.csv',append = FALSE,sep = ",", eol = "\n", na = "NA", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))
write.table(y1,file= 'c_soaked_new.csv',append = FALSE,sep = ",", eol = "\n", na = "NA", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))

#### PCA PLOTS######## refer to figures
###########################

library(readr)
library(dplyr)
library(stringr)
library(fastDummies)
library(plsgenomics)
library(ggplot2)
library(stringr)
library(FactoMineR)
library(factoextra)
library(ggfortify)

dat1 <- read.csv("n_soaked_new.csv") %>% data.frame
dat2 <- read.csv("c_soaked_new.csv") %>% data.frame
nterm<-dat1 %>% select(2,12:174)
cterm<-dat2 %>% select(2,12:174)

write.table(nterm,file='nterm.csv',append = FALSE,sep = ",", eol = "\n", na = "NA", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))
write.table(cterm,file="cterm.csv",append = FALSE,sep = ",", eol = "\n", na = "NA", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))

### read data in a data frame ####
ndat <- read.csv("nterm.csv") %>% data.frame
cdat <- read.csv("cterm.csv") %>% data.frame

###remove duplicate rows####
ndat <- ndat %>% distinct(substring, .keep_all = TRUE)
cdat <- cdat %>% distinct(substring, .keep_all = TRUE)

### make substring at 1st column as the rowname ###
row.names(ndat) <- ndat[, 1]
row.names(cdat) <- cdat[, 1]

ndat = subset(ndat, select = -c(P3__,P2__,P1__,p1__,p2__,p3__,p4__) )
cdat = subset(cdat, select = -c(P3__,P2__,P1__,p1__,p2__,p3__,p4__) )

####
write.table(ndat,file="ndat.csv",append = FALSE,sep = ",", eol = "\n", na = "NA", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))
write.table(cdat,file="cdat.csv",append = FALSE,sep = ",", eol = "\n", na = "NA", row.names = FALSE,col.names = TRUE, qmethod = c("escape", "double"))

nterm.pca <- prcomp(ndat[,-1],center = TRUE)
summary(nterm.pca)

cterm.pca <- prcomp(cdat[,-1], center = TRUE)
summary(cterm.pca)


#Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
fviz_eig(nterm.pca)

#Graph of variables. Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.
# options(ggrepel.max.overlaps = Inf)

jpeg(file="pca1&2_nterm_2days.jpeg",width=6, height=4, units="in", res=600)
scree1.plot<-fviz_pca_var(nterm.pca, labelsize = 3, pointsize = 1, geom = c("point", "text"),axes = c(1, 2), col.circle = "grey70",
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
print(scree1.plot)
dev.off()


jpeg(file="pca1&2_cterm_2days.jpeg",width=6, height=4, units="in", res=600)
scree1.plot<-fviz_pca_var(cterm.pca, labelsize = 3, pointsize = 1, geom = c("point", "text"),axes = c(1, 2), col.circle = "grey70",
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
print(scree1.plot)
dev.off()

jpeg(file="indi1&2_nterm_2days.jpeg",width=6, height=4, units="in", res=600)
scree1.plot<-fviz_pca_ind(nterm.pca, col.ind = "contrib",geom = c("point", "text"),axes = c(1, 2),
gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
repel = TRUE # Avoid text overlapping (slow if many points)
)
print(scree1.plot)
dev.off()


jpeg(file="indi1&2_cterm_2days.jpeg",width=6, height=4, units="in", res=600)
scree1.plot<-fviz_pca_ind(cterm.pca, col.ind = "contrib",geom = c("point", "text"),axes = c(1, 2),
gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
repel = TRUE # Avoid text overlapping (slow if many points)
)
print(scree1.plot)
dev.off()


jpeg(file="Dim1_nterm_soaked.jpeg",width=4, height=3, units="in", res=600)
scree1.plot<-fviz_contrib(nterm.pca, choice = "var", axes = 1, top = 10)
print(scree1.plot)
dev.off()

jpeg(file="Dim2_nterm_soaked.jpeg",width=4, height=3, units="in", res=600)
scree1.plot<-fviz_contrib(nterm.pca, choice = "var", axes = 2, top = 10)
print(scree1.plot)
dev.off()


jpeg(file="Dim1_cterm_soaked.jpeg",width=4, height=3, units="in", res=600)
scree1.plot<-fviz_contrib(nterm.pca, choice = "var", axes = 1, top = 10)
print(scree1.plot)
dev.off()

jpeg(file="Dim2_cterm_soaked.jpeg",width=4, height=3, units="in", res=600)
scree1.plot<-fviz_contrib(nterm.pca, choice = "var", axes = 2, top = 10)
print(scree1.plot)
dev.off()


# Eigenvalues
eig.val <- get_eigenvalue(nterm.pca)
eig.val
  


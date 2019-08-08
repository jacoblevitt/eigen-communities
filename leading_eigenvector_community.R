# import libraries
library(igraph)
library(reticulate)
library(readr)
library(ggplot2)
library(plyr)
library(ggpubr)
library(gridExtra)
library(ggraph)
# make edgelists
source_python('06122019_1.py')
# measure eigenvector centralities
files = dir(pattern = '*.csv')
files = files[file.size(files) > 0]
x = c('df_')
files = files[grepl(paste0(x, collapse = '|'), files)]
x = c('-')
files = files[!grepl(paste0(x, collapse = '|'), files)]
lapply(files, function(x) {
    df2 <- read.table(x, header=FALSE, sep=',') # load file
    colnames(df2) = c('v1', 'v2', 'weight')
    g1 = graph_from_data_frame(df2, directed = FALSE, vertices = NULL)
    e = eigen_centrality(g1, directed = FALSE)
    s = paste(toString(x), ',', toString(e), sep = '')
    print(s)
})
# separate catalytic and non-catalytic eigenvector centrality scores for each amino acid
source_python('06122019_2.py')
# make plots and determine statistical significance
df = read.csv('ARG.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pARG = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8)) + scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Arginine,'~italic('p = 0.000873'))



compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('HIS.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pHIS = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8)) + scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Histidine,'~italic('p = 0.0000133'))



compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('LYS.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pLYS = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Lysine,'~italic('p = 0.474'))




compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('ASP.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pASP = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Asparagine,'~italic('p = 0.000251'))



compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('GLU.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pGLU = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Glutamine,'~italic('p = 0.0979'))


compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('SER.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pSER = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Serine,'~italic('p = 0.0000298'))


compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('THR.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pTHR = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Threonine,'~italic('p = 0.0807'))


compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('ASN.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pASN = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Asparagine,'~italic('p = 0.00167'))

compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('GLN.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pGLN = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Glutamine,'~italic('p = 0.0273'))


compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('CYS.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pCYS = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Cysteine,'~italic('p = 0.0000924'))


compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('GLY.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pGLY = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Glycine,'~italic('p = 0.0419'))

          
compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('PRO.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pPRO = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Proline,'~italic('p = 0.529'))

compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('ALA.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pALA = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Alanine,'~italic('p = 0.0195'))


compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('VAL.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pVAL = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Valine,'~italic('p = 0.721'))


compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('ILE.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pILE = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Isoleucine,'~italic('p = 0.828'))


compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('LEU.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pLEU = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Leucine,'~italic('p = 0.928'))

compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('MET.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pMET = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01))+ ggtitle('Methionine,'~italic('p = 0.000169'))


compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('PHE.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pPHE = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01)) + ggtitle('Phenylalanine,'~italic('p = 0.0220'))

compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('TYR.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pTYR = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01)) + ggtitle('Tyrosine,'~italic('p = 0.00000000791'))

compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

df = read.csv('TRP.csv', sep = '=', col.names = c('Amino_Acid', 'Eigen'))
dfm <- ddply(df, "Amino_Acid", summarise, Eigen.mean=mean(Eigen))
pTRP = ggplot(df, aes(x=Eigen, fill=Amino_Acid)) +
    geom_density(alpha=.3) +
    geom_vline(data=dfm, aes(xintercept=Eigen.mean,  colour=Amino_Acid), size=0.5) +
    theme(legend.position="none") + theme(axis.title.x=element_blank()) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))+ scale_y_continuous(
  labels = scales::number_format(accuracy = 0.01)) + scale_x_continuous(
  labels = scales::number_format(accuracy = 0.01)) + ggtitle('Tryptophan,'~italic('p = 0.0000000239'))


compare_means(Eigen ~ Amino_Acid, df, method = 't.test')

grid.arrange(pTRP, pTYR, pCYS, pASP, pHIS, pGLN, pMET, pGLY, pARG, pSER, pPHE, pVAL, pILE, pPRO, pTHR, pALA, pGLU, pASN, pLEU, pLYS, nrow = 5, left = "Density", bottom = "Eigenvector Centrality")
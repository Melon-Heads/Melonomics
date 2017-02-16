######################
#  INSTALL PACKAGES  #
######################

install.packages("plotly", repos="https://cran.ma.imperial.ac.uk/")
library("plotly")

install.packages("DT", repos="https://cran.ma.imperial.ac.uk/")
library("DT")

install.packages("gplots", repos="https://cran.ma.imperial.ac.uk/")
library("gplots")

install.packages("dendextend", repos="https://cran.ma.imperial.ac.uk/")
library("dendextend")

source("http://bioconductor.org/biocLite.R") 
library("limma")

install.packages("RColorBrewer", repos="https://cran.ma.imperial.ac.uk/")
library("RColorBrewer")

#clear R environment
rm(list = ls())

#setting working directory
setwd(~/ ) 

####################
#  START ANALYSIS  #
####################

#import data
sample <- read.table("/Users/modoupehbetts/Documents/Software_development/app/data/Gene_FPKM_Sample.csv", header=TRUE, sep=",")


#select gene names 
geneNames <- as.character(sample$X)


# select required data and convert to matrix
cols <- ncol(sample)
sample <- sample[,2:cols]
sample <- as.matrix(sample)


#edit row and column names to genes names and class vector 
classVector <- as.character(sample[1,])
rownames(sample) <- geneNames
colnames(sample) <- classVector


#remove first row of data (unnecessary)
sample1 <- sample[-1,] #(keeping sample file with line of data too because volcano plot doesn't work without it)
pClass <- as.factor(colnames(sample1))


#select classes 
levels(pClass)[levels(pClass)=="1"] <- "One"
levels(pClass)[levels(pClass)=="2"] <- "Two"
levels(pClass)[levels(pClass)=="3"] <- "Three"


#index classes
index1 <- which((pClass=="One")==TRUE)
index2 <- which((pClass=="Two")==TRUE)
index3 <- which((pClass=="Three")==TRUE)


#add colour to classes
colour <- NULL
colour[index1] <- "red"
colour[index2] <- "navy"
colour[index3] <- "darkgreen"



##################################
#  HIERACHICAL CLUSTER ANALYSIS  #
##################################

hc <- hclust(dist(t(sample1)))
dend <- as.dendrogram(hc, xlab="GROUP")
colorCodes <- c(Zero="red", One="navy", Two="darkgreen")
labels_colors(dend) <- colorCodes[pClass][order.dendrogram(dend)]
pdf('HCA.pdf')
plot(dend, xlab="GROUP", ylab="HEIGHT", col.axis = "blue")
dev.off()



#############
#  HEATMAP  #
#############

pdf('heatmap.pdf')
heatmap.2(sample1, bg=colour, col=colorRampPalette(brewer.pal(9,"Blues"))(100), colsep=1:ncol(sample1), rowsep=1:nrow(sample1), sepcolor="white", sepwidth = c(0.005,0.005), scale="row", key=T, keysize=1.2, density.info="none", trace="none",cexCol=1, cex=0.85, srtCol = 0, colCol=colour, xlab='GROUP', ylab='GENE', margins = c(3.2,8.6))
dev.off()



#########
#  PCA  #
#########

#calculate PCs, 
log.sample <- log(sample1+0.1) # +0.1 needed to log transform samples of equal to 0.
sample.pca <- prcomp(t(log.sample), center = TRUE)


#summarise PCA for further analysis
s<-summary(sample.pca)


#calculate cumulative variance
cumVar <- s$importance[3,] * 100
cumVar <- as.data.frame(cumVar)


#calculate explained variance 
expVar <- s$importance[2,] * 100
expVar <- as.data.frame(expVar)


#count for number of principle components
no.rows <- nrow(expVar)


#create data frame with all variance info, including column for number of components
expVar["cumVar"] <- NA
expVar["Component"] <- c(1:no.rows)
expVar$cumVar <- cumVar$cumVar
var.rows <- nrow(expVar)



#######################
#  PCA VARIANCE PLOT  #
#######################

#Interactive Variance Plot
#set x and y axis parameters
x <- list(
  autotick=FALSE,
  tick0=0,
  dtick=2,
  tickangle=0
)
y <- list(
  title="Variance(%)",
  autotick=FALSE,
  tick0=0,
  dtick=10
)
var.plot <- plot_ly(expVar, x=~Component)%>%
  add_trace(y=~cumVar, type="scatter", name="Cumulative PCA", mode="lines", line=list(color="rgb(205, 12, 24)"))%>%
  add_trace(y=~expVar, type="scatter", name="PCA", mode="lines", line=list(color ="rgb(22, 96, 167)"))%>%
  layout(xaxis=x, yaxis=y,legend=list(x=0.7,y=0.15))
ggplotly(var.plot)


#save interactive variance plot 
htmlwidgets::saveWidget(var.plot, "variance.html", selfcontained=FALSE)



#####################
#  PCA SCORES PLOT  #
#####################

#create data frame with PCA scores
Xscores <- sample.pca$x 
Xscores <- as.data.frame(Xscores)


#graph(s) showing PC1 vs 2,3,4,5
scores12 <- plot_ly(Xscores,x=~PC1, y=~PC2,  type="scatter", mode="markers", text=~paste('GROUP: ', classVector), color=classVector, showlegend=TRUE)
scores13 <- plot_ly(Xscores, x=~PC1, y=~PC3, type="scatter", mode="markers", text=~paste('GROUP: ', classVector), color=classVector, showlegend=FALSE)
scores14 <- plot_ly(Xscores, x=~PC1, y=~PC4, type="scatter", mode="markers", text=~paste('GROUP: ', classVector), color=classVector, showlegend=FALSE)
scores15 <- plot_ly(Xscores, x=~PC1, y=~PC5, type="scatter", mode="markers", text=~paste('GROUP: ', classVector), color=classVector, showlegend=FALSE)
scores1 <- subplot(scores12, scores13, scores14, scores15, nrows=4, shareX =TRUE, shareY = TRUE)
ggplotly(scores1)

htmlwidgets::saveWidget(scores1, "scores1.html", selfcontained=FALSE)


#graph(s) showing PC2 vs 3,4,5
scores23 <- plot_ly(Xscores, x=~PC2, y=~PC3, type="scatter", mode="markers", text=~paste('GROUP: ', classVector), color=classVector, showlegend=TRUE)
scores34 <- plot_ly(Xscores, x=~PC2, y=~PC4, type="scatter", mode="markers", text=~paste('GROUP: ', classVector), color=classVector, showlegend=FALSE)
scores45 <- plot_ly(Xscores, x=~PC2, y=~PC5, type="scatter", mode="markers", text=~paste('GROUP: ', classVector), color=classVector, showlegend=FALSE)
scores2 <- subplot(scores23, scores34, scores45, nrows=3, shareX =TRUE, shareY = TRUE)
ggplotly(scores2)

htmlwidgets::saveWidget(scores2, "scores2.html", selfcontained=FALSE)


#graphs(s) showing PC3 vs 4,5
scores34 <- plot_ly(Xscores, x=~PC3, y=~PC4, type="scatter", mode="markers", text=~paste('GROUP: ', classVector), color=classVector, showlegend=TRUE)
scores45 <- plot_ly(Xscores, x=~PC3, y=~PC5, type="scatter", mode="markers", text=~paste('GROUP: ', classVector), color=classVector, showlegend=FALSE)
scores3 <- subplot(scores34, scores45, nrows=2, shareX =TRUE, shareY = TRUE)
ggplotly(scores3)

htmlwidgets::saveWidget(scores3, "scores3.html", selfcontained=FALSE)


#graphs(s) showing PC4 vs PC5 
scores4 <- plot_ly(Xscores, x=~PC4, y=~PC5, type="scatter", mode="markers", text=~paste('GROUP: ', classVector), color=classVector)
ggplotly(scores4)

htmlwidgets::saveWidget(scores4, "scores4.html", selfcontained=FALSE)



#########################
#  TESTS FOR TOP GENES  #
#########################

#carry out tests to figure out top genes and for volcano plots
design <- model.matrix(~0+pClass)


#test for top genes
fit1 <- lmFit(sample1, design)
contrasts1 <- makeContrasts(pClassOne - pClassTwo - pClassThree, levels = design)
fit1 <- contrasts.fit(fit1, contrasts1)
fit1 <- eBayes(fit1)


#create table of top genes
toptable <- topTable(fit1, sort.by="p", number=20)



########################
#  TABLE OF TOP GENES  #
########################

#Interactive top genes table
datatable(toptable)
topgenes <- DT::datatable(toptable)

#save interactive table
DT::saveWidget(topgenes, "toptable.html", selfcontained = FALSE)



##################
#  VOLCANO PLOT  #
##################

#test for top genes
fit <- lmFit(sample, design)
contrasts <- makeContrasts(pClassOne - pClassTwo - pClassThree, levels = design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)


#create table of top genes
toptable <-topTable(fit, sort.by="p", number=250)


#create log function for using in volcano plot
lg <- -log10(toptable$P.Value)


#Interactive Volcano plot
volcano <- plot_ly(toptable, x=~logFC, y=~lg, type="scatter", mode="markers", text = ~paste('GENE: ', geneNames), marker=list(color="hotpink"))%>%
  layout(xaxis=list(range=c(-2,2)), yaxis=list(title="-log10(P Value)"))
ggplotly(volcano)

#save volcano plot
htmlwidgets::saveWidget(volcano, "volcanoplot.html", selfcontained=FALSE)


#N.B. tests for top genes were done twice because editted data set did not work for volcano plot, while uneditted version produced an innaccurate datatable for top genes. 

#clear R environment
rm(list = ls())

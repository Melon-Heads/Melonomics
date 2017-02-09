# INSTALL PACKAGES
install.packages("pheatmap")
install.packages("plotly")
install.packages("dendextend")
source("http://bioconductor.org/biocLite.R") 
library("limma")

#import data
sample <- read.table("sample_file.csv", header=TRUE, sep=",")

#edit column names  (just using for sample file, will change later)
colnames(sample)[2] <- "Ohour"
colnames(sample)[3] <- "24hour"
colnames(sample)[4] <- "8hour"

geneNames <- as.character(sample$X)

# select required data and convert to matrix
sample <- sample[,2:4]
sample <- as.matrix(sample)

#change row names to genes names
rownames(sample) <- geneNames
View(sample)

##################################
#  HIERACHICAL CLUSTER ANALYSIS  #
##################################

hc <- hclust(dist(sample))
pdf('HCA.pdf')
plot(hc, hang = -1)
dev.off()


#############
#  HEATMAP  #
#############

pdf('heatmap.pdf')
pheatmap(sample, color=colorRampPalette(c("yellow","yellow2","greenyellow","green","aquamarine","turquoise","steelblue","navy"))(100))
dev.off()

#########
#  PCA  #
#########

#calculate PCs, 
log.sample <- log(sample+0.1) # +0.1 needed to log transform samples of equal to 0.
sample.gene <- sample[,1]
sample.pca <- prcomp(log.sample, center = TRUE, scale. = TRUE)

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

#plot combined variance graph using plotly
x <- list(
  autotick=FALSE,
  tick0=0,
  dtick=5,
  tickangle=0
)
y <- list(
  title="Variance(%)",
  autotick=FALSE,
  tick0=0,
  dtick=10
)
var.plot <- plot_ly(expVar, x=~Component)%>%
  add_trace(y=~cumVar, type="scatter", name="Cumulative Variance", mode="lines", line=list(color="rgb(205, 12, 24)"))%>%
  add_trace(y=~expVar, type="scatter", name="Explained Variance", mode="lines", line=list(color ="rgb(22, 96, 167)"))%>%
  layout(xaxis=x, yaxis=y,legend=list(x=0.7,y=0.15))
ggplotly(var.plot)
htmlwidgets::saveWidget(var.plot, "variance.html")

#plot graph comparing PCs
Xscores <- sample.pca$x 
Xscores <- as.data.frame(Xscores)
scores <- plot_ly(Xscores, x=~PC1, y=~PC2, type="scatter", mode="markers")
ggplotly(scores)
htmlwidgets::saveWidget(scores, "PCAscores.html")

#plot graph showing PC loadings
Xloadings <- sample.pca$rotation
Xloadings <- as.data.frame(Xloadings)
loadings <- plot_ly(Xloadings, x=~PC1, y=~PC2, type="scatter")
ggplotly(loadings)
htmlwidgets::saveWidget(loadings, "PCAloadings.html")

########################
#  TABLE OF TOP GENES  #
########################

#currently work in progress, waiting for something from Nadim 


##################
#  VOLCANO PLOT  #
##################

#see above




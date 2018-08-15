##### CODE TO GENERATE NETWORK PLOTS

### TO RUN: CHANGE DIRECTORY NAME TO RELEVENT LOCATION AND BATCH_NAME TO RELEVENT NAME (i.e. NAME OF GROUP OF SAMPLES), THEN RUN CODE

dir = '/USER/DIR/'
sample = 'Sample_1'

###########################################################################################################################################################
concat = function(v) {
	res = ""
	for (i in 1:length(v)){
	res = paste(res,v[i],sep="")
	}
	res
}

igraphgeneration <- function(sample,dir){
	n=concat(c(dir,"Plot_ids_", sample,".txt"))
	e=concat(c(dir,"Edges_", sample,".txt"))
	nodes <- read.csv(n, head=FALSE, sep="\t")
	edge <- read.csv(e, head=FALSE, sep="\t")
	library(igraph)
	g <- graph.empty( n=0, directed=FALSE)
	freq<-nodes[,3]
	max_freq=(sum(freq))
	frequency<-freq*50/max_freq
	colour1 <- heat.colors(max_freq/2+2)
	colour <-rev(colour1)
	g <- igraph::add.vertices(g, length(nodes[,2]), name=as.character(nodes[, 2]),color = colour[freq/2+2])
	names <- V(g)$name
	ids <- 1:length(names)
	names(ids) <- names
	from <- as.character(edge[,1])
	to <- as.character(edge[,2])
	edges <- matrix(c(ids[from], ids[to]), nc=2)
	g <- add.edges(g, t(edges), weight=1)
	V(g)$label <- V(g)$name
	V(g)$size<-freq
	V(g)$label.cex<-0.0001
	del_ids<-intersect(which(degree(g)==0), which(freq==1))
	g1<-delete.vertices(g,ids[del_ids])
	g1
}
fileout = "Fileout.jpeg"
jpeg(file=fileout, height=7016, width=7016, res=600)
par(mfrow= c(1,1), mar = c(5,5,5,5))
title1=sample
g1<-igraphgeneration(sample,dir)
layout = layout_with_graphopt(g1,niter = 800)
### layout=layout.fruchterman.reingold(g1, niter = 400) ## used with the old version of igraph
plot(g1, layout=layout, edge.color="darkblue", main=title1, edge.width=0.2)
dev.off()

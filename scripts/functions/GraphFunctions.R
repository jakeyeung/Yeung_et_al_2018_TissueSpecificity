# 2016-11-04
# Graph functions

showSigNodes2 <-
  function(DAG, sigTerm, sigTerm_Local, sigTerm_Global, dagTermInfo) {
    
    require('Rgraphviz') || stop('package Rgraphviz is required')
    
    graphAttrs <- getDefaultAttrs(layoutType = 'dot')
    graphAttrs$cluster <- NULL
    graphAttrs$node$shape <- 'ellipse'
    graphAttrs$node$fontsize <- '20'
    
    nodeAttrs <- list()
    edgeAttrs <- list()
    
    allTerm <- as.character(dagTermInfo[,1])
    nodeAttrs$label[allTerm] <- allTerm
    
    rmLocalTerm <- setdiff(sigTerm, sigTerm_Local)
    nodeAttrs$color[rmLocalTerm] <- rep('red', length(rmLocalTerm))
    nodeAttrs$shape[rmLocalTerm] <- rep('circle', length(rmLocalTerm))
    
    rmGlobalTerm <- setdiff(sigTerm_Local, sigTerm_Global)
    nodeAttrs$color[rmGlobalTerm] <- rep('red', length(rmGlobalTerm))
    nodeAttrs$shape[rmGlobalTerm] <- rep('box', length(rmGlobalTerm))
    nodeAttrs$height[rmGlobalTerm] <- rep('0.7', length(rmGlobalTerm))
    nodeAttrs$width[rmGlobalTerm] <- rep('0.7', length(rmGlobalTerm))
    
    nodeAttrs$color[sigTerm_Global] <- rep('red', length(sigTerm_Global))
    nodeAttrs$shape[sigTerm_Global] <- rep('rectangle', length(sigTerm_Global))
    nodeAttrs$height[sigTerm_Global] <- rep('0.7', length(sigTerm_Global))
    nodeAttrs$width[sigTerm_Global] <- rep('1.1', length(sigTerm_Global))
    
    dagTermInfo[dagTermInfo[,5]<2.2E-16,5] <- 2.2E-16;
    dagTermInfo[dagTermInfo[,6]<2.2E-16,6] <- 2.2E-16;
    dagTermInfo$colorran <- round(log10(dagTermInfo[,6])-range(log10(dagTermInfo[,6]))[1] + 1)
    mm <- max(dagTermInfo$colorran)
    colorMap <- heat.colors(mm)
    nodeAttrs$fillcolor[allTerm] <- unlist(lapply(dagTermInfo$colorran, function(x) return(colorMap[x])))
    
    weightsList <- edgeWeights(DAG)
    to <- lapply(weightsList, names)
    from <- nodes(DAG)
    edge.names <- paste(rep(from, listLen(to)), unlist(to), sep = "~")
    edge.weights <- unlist(weightsList)
    names(edge.weights) <- edge.names
    ##    0 for a is_a relation,  1 for a part_of relation
    edgeAttrs$color <- ifelse(edge.weights == 0, 'black', 'red')
    plot(DAG, attrs = graphAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs)
  }
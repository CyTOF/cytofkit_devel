## loading pacakge
require(cytofkit)
require(ggplot2)
require(reshape2)

## load library for cytof_expressionTrends
library(VGAM)
library(plyr)


## Main function for scatter plot
scatterPlot <- function(obj, plotMethod, plotFunction, pointSize=1, 
                      addLabel=TRUE, labelSize=1, sampleLabel = TRUE,
                      FlowSOM_k = 40, selectSamples, facetPlot = FALSE, 
                      colorPalette = "bluered", labelRepel = FALSE, removeOutlier = TRUE){
    
    data <- cbind(obj$expressionData, 
                  obj$dimReducedRes[[plotMethod]], 
                  do.call(cbind, obj$clusterRes))
    data <- as.data.frame(data)
    xlab <- colnames(obj$dimReducedRes[[plotMethod]])[1]
    ylab <- colnames(obj$dimReducedRes[[plotMethod]])[2]
    row.names(data) <- row.names(obj$expressionData)
    
    clusterMethods <- names(obj$clusterRes)
    samples <- sub("_[0-9]*$", "", row.names(obj$expressionData))
    data <- data[samples %in% selectSamples, ]
    nsamples <- samples[samples %in% selectSamples]
    data$sample <- nsamples
    sample_num <- length(unique(nsamples))

    if(plotFunction == "DensityPlot"){
        colPalette <- colorRampPalette(c("blue", "turquoise", "green", 
                                         "yellow", "orange", "red"))
        densCol <- densCols(data[, c(xlab, ylab)], colramp = colPalette)
        data$densCol <- densCol
        gp <- ggplot(data, aes_string(x=xlab, y=ylab)) + 
            geom_point(colour=densCol, size = pointSize) + ggtitle("Density Plot") +
            theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
    }else if(plotFunction == "DotPlot"){
        gp <- ggplot(data, aes_string(x=xlab, y=ylab)) + 
            geom_point(size = pointSize) + ggtitle("Dot Plot") +
            xlab(xlab) + ylab(ylab) + theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
    }else if(plotFunction == "ColorBySample"){
        size_legend_row <- ceiling(sample_num/4)
        sample <- "sample"
        gp <- ggplot(data, aes_string(x=xlab, y=ylab, colour = sample)) +
            geom_point(size = pointSize) + ggtitle("Color By Sample") +
            xlab(xlab) + ylab(ylab) + theme_bw() + theme(legend.position = "bottom") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold")) +
            guides(colour = guide_legend(nrow = size_legend_row, override.aes = list(size = 4)))
    }else if(plotFunction == "FacetByMarker"){
        gp <- cytof_wrap_colorPlot(data = data, 
                              xlab = xlab, 
                              ylab = ylab, 
                              markers = colnames(obj$expressionData), 
                              colorPalette = colorPalette,
                              pointSize = pointSize, 
                              removeOutlier = TRUE)
        
    }else if(plotFunction %in% clusterMethods){
        gp <- cytof_clusterPlot(data = data, 
                                xlab = xlab, 
                                ylab = ylab, 
                                cluster = plotFunction, 
                                sample = "sample",
                                title = plotFunction, 
                                type = ifelse(facetPlot, 2, 1),
                                point_size = pointSize, 
                                addLabel = addLabel, 
                                labelSize = labelSize, 
                                sampleLabel = sampleLabel,
                                labelRepel = labelRepel,
                                fixCoord = FALSE)
    }else{
        gp <- cytof_colorPlot(data = data, 
                              xlab = xlab, 
                              ylab = ylab, 
                              zlab = plotFunction, 
                              colorPalette = colorPalette,
                              pointSize = pointSize, 
                              removeOutlier = TRUE)
    }
    
    return(gp)
}

## Facet wrap plot of marker exporession
cytof_wrap_colorPlot <- function(data, xlab, ylab, markers, 
                            colorPalette = c("bluered", "topo", "heat", "terrain", "cm"), 
                            pointSize=1, 
                            removeOutlier = TRUE){
    
    remove_outliers <- function(x, na.rm = TRUE, ...) {
        qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
        H <- 1.5 * IQR(x, na.rm = na.rm)
        y <- x
        y[x < (qnt[1] - H)] <- qnt[1] - H
        y[x > (qnt[2] + H)] <- qnt[2] + H
        y
    }
    
    data <- as.data.frame(data)
    title <- "Marker Expression Level Plot"
    data <- data[,c(xlab, ylab, markers)]
    
    if(removeOutlier){
        for(m in markers){
            data[[m]] <- remove_outliers(data[ ,m])
        }
    }
    
    data <- melt(data, id.vars = c(xlab, ylab), 
                 measure.vars = markers,
                 variable.name = "markers", 
                 value.name = "Expression")

    colorPalette <- match.arg(colorPalette)
    switch(colorPalette,
           bluered = {
               myPalette <- colorRampPalette(c("blue", "white", "red"))
           },
           topo = {
               myPalette <- colorRampPalette(topo.colors(50))
           },
           heat = {
               myPalette <- colorRampPalette(heat.colors(50))
           },
           terrain = {
               myPalette <- colorRampPalette(terrain.colors(50))
           },
           cm = {
               myPalette <- colorRampPalette(cm.colors(50))
           }
    )
    zlength <- nrow(data)
    grid_col_num <- round(sqrt(length(markers)))
    gp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = "Expression")) + 
        facet_wrap(~markers, ncol = grid_col_num, scales = "fixed") +
        scale_colour_gradientn(name = "Expression", colours = myPalette(zlength)) +
        geom_point(size = pointSize) + theme_bw() + coord_fixed() +
        theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + ggtitle(title) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
    
    return(gp)
}

## Heat Map
heatMap <- function(data, clusterMethod = "DensVM", type = "mean", selectSamples,
                    cex_row_label = 1, cex_col_label = 1, scaleMethod = "none") {
    exprs <- data$expressionData
    samples <- sub("_[0-9]*$", "", row.names(exprs))
    mySamples <- samples %in% selectSamples
    exprs <- exprs[mySamples, , drop = FALSE]
    dataj <- data$clusterRes[[clusterMethod]][mySamples]
    exprs_cluster <- data.frame(exprs, cluster = dataj, check.names = FALSE )
    
    cluster_stat <- cytof_clusterStat(data = exprs_cluster,
                             cluster = "cluster", 
                             statMethod = type)
    
    cytof_heatmap(data = as.matrix(cluster_stat), 
                  baseName = paste(clusterMethod, type), 
                  scaleMethod = scaleMethod, 
                  cex_row_label = cex_row_label, 
                  cex_col_label = cex_col_label,
                  margins = c(8, 8), 
                  keysize = 1, 
                  key.par=list(mgp=c(2, 1, 0), mar=c(4, 3, 4, 0))) 
}


## Combined marker expression trend
cytof_expressionTrends <- function(data, markers, clusters, 
                                  orderCol="isomap_1", 
                                  clusterCol = "cluster", 
                                  reverseOrder = FALSE,
                                  addClusterLabel = TRUE,
                                  clusterLabelSize = 5,
                                  segmentSize = 0.5,
                                  min_expr = NULL, 
                                  trend_formula="expression ~ sm.ns(Pseudotime, df=3)"){
    
    if(!is.data.frame(data)) data <- data.frame(data, check.names = FALSE)
    if(!all(markers %in% colnames(data))) stop("Unmatching markers found!")
    if(!(length(orderCol)==1 && orderCol %in% colnames(data)))
        stop("Can not find orderCol in data!")
    if(!(length(clusterCol)==1 && clusterCol %in% colnames(data)))
        stop("Can not find clusterCol in data!")
    if(!missing(clusters)){
        if(!all(clusters %in% data[[clusterCol]]))
            stop("Wrong clusters selected!")
        data <- data[data[[clusterCol]] %in% clusters, , drop=FALSE]
    }
    
    if(reverseOrder){
        newOrderCol <- paste0(orderCol, "(reverse)")
        data[[newOrderCol]] <- -data[[orderCol]]
        orderCol <- newOrderCol
    }
    orderValue <- data[[orderCol]]
    data <- data[order(orderValue), c(markers, clusterCol)]
    data$Pseudotime <- sort(orderValue)
    
    mdata <- melt(data, id.vars = c("Pseudotime", clusterCol), 
                  variable.name = "markers", value.name= "expression")
    colnames(mdata) <- c("Pseudotime", clusterCol, "markers", "expression")
    mdata$markers <- factor(mdata$markers)
    mdata[[clusterCol]] <- factor(mdata[[clusterCol]])
    min_expr <- min(mdata$expression)
    
    ## tobit regression
    vgamPredict <- ddply(mdata, .(markers), function(x) { 
        fit_res <- tryCatch({
            vg <- suppressWarnings(vgam(formula = as.formula(trend_formula), 
                                        family = VGAM::tobit(Lower = min_expr, lmu = "identitylink"), 
                                        data = x, maxit=30, checkwz=FALSE))
            res <- VGAM::predict(vg, type="response")
            res[res < min_expr] <- min_expr
            res
        }
        ,error = function(e) {
            print("Error!")
            print(e)
            res <- rep(NA, nrow(x))
            res
        }
        )
        expectation = fit_res
        data.frame(Pseudotime=x[["Pseudotime"]], expectation=expectation)
    })
    
    color_by <- clusterCol
    plot_cols <- round(sqrt(length(markers)))
    cell_size <- 1
    x_lab <- orderCol
    y_lab <- "Expression"
    legend_title <- "Cluster"
    
    ## copied from monocle package
    monocle_theme_opts <- function(){
        theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
            #theme(panel.border = element_blank(), axis.line = element_line()) +
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
            theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
            theme(panel.background = element_rect(fill='white')) +
            theme(legend.position = "right") +
            theme(axis.title = element_text(size = 15)) +
            theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))}
    
    q <- ggplot(data=vgamPredict, aes_string(x="Pseudotime", y="expectation", col="markers")) + geom_line(size = 1.5)
    q <- q + ylab(y_lab) + xlab(x_lab) + theme_bw()
    q <- q + guides(colour = guide_legend(title = legend_title, override.aes = list(size = cell_size*3)))
    q <- q + monocle_theme_opts() 
    
    # if(addClusterLabel){
    #     # edata <- data[ ,c("Pseudotime", clusterCol)]
    #     # colnames(edata) <- c('x', "z")
    #     # center <- aggregate(x ~ z, data = edata, median)
    #     # center$y <- -0.5 ## add to the botom
    #     # q <- q + geom_text_repel(data=center, aes(x=x, y=y, label=z), parse=TRUE)
    #     mdata$cluster <- mdata[[clusterCol]]
    #     center <- aggregate(cbind(Pseudotime, expression) ~ cluster + markers, data = mdata, median)
    #     q <- q + geom_text_repel(data=center, aes(x=Pseudotime, y=expression, label=cluster),
    #                              size = clusterLabelSize, fontface = 'bold',
    #                              box.padding = unit(0.5, 'lines'),
    #                              point.padding = unit(1.6, 'lines'),
    #                              segment.color = '#555555',
    #                              segment.size = segmentSize,
    #                              arrow = arrow(length = unit(0.02, 'npc')))
    # }
    
    q
}

## function for opening the results directory
opendir <- function(dir = getwd()){
    if (.Platform['OS.type'] == "windows"){
        shell.exec(dir)
    } else {
        system(paste(Sys.getenv("R_BROWSER"), dir))
    }
}



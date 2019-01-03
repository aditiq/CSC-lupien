#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.4.1
# mordor
# Objective : Scripts edited for aesthetics for cisTopic
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

library(cisTopic)

cellTopicHeatmap2 <-  function (object, method = "Zscore", colorBy = NULL, colVars = list(),
  col.low = "dodgerblue", col.mid = "floralwhite", col.high = "brown1",
  ...){

  if (method == "Zscore") {
      topic.mat <- scale(object@selected.model$document_expects,
          center = TRUE, scale = TRUE)
  }
  else if (method == "probability") {
      alpha <- 50/length(cisTopicObject@selected.model$topic_sums)
      topic.mat <- apply(object@selected.model$document_sums,
          2, function(x) {
              (x + alpha)/sum(x + alpha)
          })
  }
  rownames(topic.mat) <- paste("Topic", seq(1, nrow(topic.mat)))
  colnames(topic.mat) <- object@cell.names
  cl.cells <- hclust.vector(t(topic.mat), method = "ward",
      metric = "euclidean")
  dd.cells <- as.dendrogram(cl.cells)
  colorPal <- grDevices::colorRampPalette(c(col.low, col.mid,
      col.high))
  object.cell.data <- object@cell.data
  if (is.null(colorBy)) {
      colorBy <- colnames(object.cell.data)[-which(as.vector(sapply(object.cell.data,
          is.numeric)))]
  }
  for (variable in colorBy) {
      if (is.null(colVars[[variable]])) {
          colVars[[variable]] <- setNames(rainbow(length(unique(object.cell.data[,
              variable])),s=0.5), as.vector(sort(unique(object.cell.data[,
              variable]))))
          cellColor <- setNames(colVars[[variable]][object.cell.data[,
              variable]], object@cell.names)
      }
  }
  nmf.options(grid.patch = TRUE)
  NMF::aheatmap(topic.mat, scale = "none", revC = TRUE, main = "cisTopic contributions per cell",
      sub = "Column normalized topic contribution", Colv = dd.cells,
      annCol = object.cell.data[object@cell.names, colorBy,
          drop = FALSE], annColor = colVars, labCol = NA, color = colorPal(20))
}

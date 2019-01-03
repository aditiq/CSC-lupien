make_cicero_cds <- function(cds,
                            reduced_coordinates,
                            k=50,
                            summary_stats = NULL,
                            size_factor_normalize = TRUE,
                            silent = FALSE) {

  assertthat::assert_that(is(cds, "CellDataSet"))
  assertthat::assert_that(is.data.frame(reduced_coordinates) |
                            is.matrix(reduced_coordinates))
  assertthat::assert_that(assertthat::are_equal(nrow(reduced_coordinates),
                                                nrow(pData(cds))))
  assertthat::assert_that(setequal(row.names(reduced_coordinates),
                                   colnames(cds)))
  assertthat::assert_that(assertthat::is.count(k) & k > 1)
  assertthat::assert_that(is.character(summary_stats) | is.null(summary_stats))
  if(!is.null(summary_stats)) {
    assertthat::assert_that(all(summary_stats %in% names(pData(cds))),
                            msg = paste("One of your summary_stats is missing",
                                        "from your pData table. Either add a",
                                        "column with the name in",
                                        "summary_stats, or remove the name",
                                        "from the summary_stats parameter.",
                                        collapse = " "))
    assertthat::assert_that(sum(vapply(summary_stats, function(x) {
      !(is(pData(cds)[,x], "numeric") | is(pData(cds)[,x], "integer"))}, 1)) == 0,
                            msg = paste("All columns in summary_stats must be",
                                        "of class numeric or integer.",
                                        collapse = " "))
  }
  assertthat::assert_that(is.logical(size_factor_normalize))
  assertthat::assert_that(is.logical(silent))

  reduced_coordinates <- as.data.frame(reduced_coordinates)

  reduced_coordinates <- reduced_coordinates[colnames(cds),]

  # Create a k-nearest neighbors map
  nn_map <- as.data.frame(FNN::knn.index(reduced_coordinates, k=(k-1)))

  row.names(nn_map) <- row.names(reduced_coordinates)
  nn_map$agg_cell <- seq_len(nrow(nn_map))

  good_choices <- seq_len(nrow(nn_map))
  choice <- sample(seq_len(length(good_choices)), size = 1, replace = FALSE)
  chosen <- good_choices[choice]
  good_choices <- good_choices[good_choices != good_choices[choice]]
  it <- 0

  while (length(good_choices) > 0 & it < 5000) { # slow
    it <- it + 1
    choice <- sample(seq_len(length(good_choices)), size = 1, replace = FALSE)
    new_chosen <- c(chosen, good_choices[choice])
    good_choices <- good_choices[good_choices != good_choices[choice]]
    cell_sample <- nn_map[new_chosen,]
    combs <- data.frame(seq_len(nrow(cell_sample)-1), nrow(cell_sample))
    shared <- apply(combs, 1, function(x) {
      (k * 2) - length(unique(as.vector(as.matrix(cell_sample[x,]))))
    })

    if (max(shared) < .9 * k) {
      chosen <- new_chosen
    }
  }

  cell_sample <- nn_map[chosen,]
  combs <- t(combn(nrow(cell_sample), 2))
  shared <- apply(combs, 1, function(x) {  #slow
    (k * 2) - length(unique(as.vector(as.matrix(cell_sample[x,]))))
  })

  if(!silent) {
  message(paste0("Overlap QC metrics:\nCells per bin: ", k,
                 "\nMaximum shared cells bin-bin: ", max(shared),
                 "\nMean shared cells bin-bin: ", mean(shared),
                 "\nMedian shared cells bin-bin: ", median(shared)))

  if (mean(shared)/k > .1) warning("On average, more than 10% of cells are
                                   shared between paired bins.")
  }

  exprs_old <- exprs(cds)

  x <- lapply(seq_len(nrow(cell_sample)), function(x) {
    return(Matrix::rowSums(exprs_old %*%
                             Matrix::Diagonal(x=seq_len(ncol(exprs_old)) %in%
                                                cell_sample[x,])))})
  new_exprs <- do.call(rbind, x)

  pdata <- pData(cds)
  new_pdata <- plyr::adply(cell_sample, 1, function(x) {
    sub <- pdata[as.numeric(x[1,]),]
    df_l <- list()
    df_l["temp"] <- 1
    for (att in summary_stats) {
      df_l[paste0("mean_", att)] <- mean(sub[,att])
    }
    data.frame(df_l)
  })

  new_pdata <- new_pdata[,(k):ncol(new_pdata)]

  new_pdata$agg_cell <- paste("agg", new_pdata$agg_cell, sep="")

  row.names(new_pdata) <- new_pdata$agg_cell
  row.names(new_exprs) <- new_pdata$agg_cell
  new_exprs <- as.matrix(t(new_exprs))

  fdf <- fData(cds)
  new_pdata$temp <- NULL

  fd <- new("AnnotatedDataFrame", data = fdf)
  pd <- new("AnnotatedDataFrame", data = new_pdata)

  cicero_cds <-  suppressWarnings(newCellDataSet(new_exprs,
                                phenoData = pd,
                                featureData = fd,
                                expressionFamily=negbinomial.size(),
                                lowerDetectionLimit=0))

  cicero_cds <- monocle::detectGenes(cicero_cds, min_expr = .1)
  cicero_cds <- BiocGenerics::estimateSizeFactors(cicero_cds)
  cicero_cds <- suppressWarnings(BiocGenerics::estimateDispersions(cicero_cds))

  if (any(!c("chr", "bp1", "bp2") %in% names(fData(cicero_cds)))) {
    fData(cicero_cds)$chr <- NULL
    fData(cicero_cds)$bp1 <- NULL
    fData(cicero_cds)$bp2 <- NULL
    fData(cicero_cds) <- cbind(fData(cicero_cds),
                               df_for_coords(row.names(fData(cicero_cds))))
  }

  if (size_factor_normalize) {
    Biobase::exprs(cicero_cds) <-
      t(t(Biobase::exprs(cicero_cds))/Biobase::pData(cicero_cds)$Size_Factor)
  }

  cicero_cds
}

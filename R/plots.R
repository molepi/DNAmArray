utils::globalVariables(c(
  "qcProbes", "res", "value", "Row.names", 
  "gcscatterplot", "column_to_rownames", "rownames_to_column"
))

##' basic GC scatter plot
##'
##' @title quality control plot
##' @param data input data
##' @param x x
##' @param y y
##' @param col col
##' @param xlab xlab
##' @param ylab ylab
##' @param threshold_h threshold_h
##' @param threshold_v threshold_v
##' @param outliers outliers
##' @param xlim xlim
##' @param ylim ylim
##' @param main main
##' @return plot
##' @author tkuipers
##' @export
##' @import ggplot2
.gcscatterplot <- function(data,
                          x,
                          y,
                          col,
                          xlab = NULL,
                          ylab = NULL,
                          threshold_h = NULL,
                          threshold_v = NULL,
                          outliers = NULL,
                          xlim,
                          ylim,
                          main) {
  
  if (nrow(unique(data[col])) > 25) {
    if (all(!is.na(as.numeric(data[, col])))) {
      data[, col] <- as.numeric(data[, col])
    } else {
      col <- NULL
    }
  }
    
  out <- data[data$Row.names %in% outliers, ]
  data <- data[!(data$Row.names %in% outliers), ]
  ggplot(data = data, aes_string(x = x, y = y, color = col)) +
    geom_point(
      size = 1,
      shape = 16
    ) +
    {
      if (is.numeric(data[col]))
        scale_color_gradientn(
          colours = c("lightgrey", "darkblue"),
          guide = guide_colorbar(nbin = 10, barheight = 20)
        )
    } +
    geom_point(
      data = out,
      aes_string(x = x, y = y),
      shape = 4
    ) +
    geom_text_repel(
      data = out,
      aes(label = Row.names),
      show.legend = FALSE,
      box.padding = 0.5,
      size=3,
      segment.color = NA
    ) +
    geom_hline(
      yintercept = threshold_h,
      color = "grey60",
      size = 1,
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = threshold_v,
      color = "grey60",
      size = 1,
      linetype = "dashed"
    ) +
    theme_bw() +
    ggtitle(main) +
    xlab(xlab) +
    ylab(ylab)
}

##' rotate data
##'
##' @title rotate data
##' @param data input data
##' @param columns col
##' @return data
##' @author tkuipers
##' @export

.rotateData <- function(data, columns) {
  data[, columns] <-
    c(0.5 * (data[, columns[1]] + data[, columns[2]]),
      data[, columns[1]] - data[, columns[2]])
  data
}

##' Rotated MU plot
##'
##' @title Rotated MU plot
##' @description
##' Plots the median methylated and unmethylated log2 intensities per sample,
##' flagging low intensity outliers.
##' 
##' @param object input sData
##' @param threshold array-specific threshold
##' @param col variable in targets to colour by
##' @return plot
##' @author tkuipers, ljsinke
##' @export
##' @import ggplot2

plotMU <- function(object,
                   threshold,
                   col = NULL) {
  MU <- log2(t(object@MU))
  targets <- object@targets
  data <- merge(MU, targets, by = "row.names")
  data <- .rotateData(data, columns = c("Methylated", "Unmethylated"))
  
  ## Get outliers
  outliers <- data$Row.names[data$Methylated <= threshold]
  
  ## Plot
  .gcscatterplot(
    data,
    x = "Methylated",
    y = "Unmethylated",
    col = col,
    threshold_v = threshold,
    outliers = outliers,
    xlab = expression(paste(log[2], sqrt(M %*% U))),
    ylab = expression(paste(log[2], "(", M / U, ")")),
    main = "Rotated MU plot"
  )
}

##' Sample-dependent overall quality control
##'
##' @title Sample-dependent overall quality control
##' @description
##' NP control probes query each of the four nucleotides in a non-polymorphic region of the bisulfite genome. 
##' Signal intensity from these probes in the red (A and T) and green (C and G) channels can then be used to 
##' test overall performance of the assay, from amplification to detection. Intensity should be high in the 
##' red channel for NP probes querying A and T nucleotides, and high in the green channel for NP probes 
##' querying G and C nucleotides. The intensities for the relevant channel for NP probes are then combined 
##' and plotted per sample, where they should be above the specific threshold.
##' 
##' @param object input sData
##' @param threshold array-specific threshold
##' @param col variable in targets to colour by
##' @return plot
##' @author tkuipers, ljsinke
##' @export
##' @import ggplot2

plotOP <- function(object,
                   threshold,
                   col = NULL) {
  data <- object@plotdata
  
  d <- data[grepl(qcProbes["NP"], data$Type), ]
  dGrn <- d[d$ExtendedType %in% c("NP (C)", "NP (G)"), c(1:5, 7)]
  x <- tapply(dGrn$IntGrn, dGrn$Samples, mean)
  dRed <- d[d$ExtendedType %in% c("NP (A)", "NP (T)"), c(1:6)]
  y <- tapply(dRed$IntRed, dRed$Samples, mean)
  
  data <- data.frame(x, y)
  targets <- object@targets
  data <- merge(data, targets, by = "row.names", suffixes = c("", ".y"))
  
  data <- .rotateData(data, columns = c("x", "y"))
  
  ## Get outliers
  outliers <- data$Row.names[data$x <= threshold]
  
  .gcscatterplot(
    data,
    x = "x",
    y = "y",
    col = col,
    threshold_v = threshold,
    outliers = outliers,
    ylab = expression(paste(log[2], "(", R / G, ")")),
    xlab = expression(paste(log[2], sqrt(R %*% G))),
    main = "Sample-dependent overall quality control (NP)"
  )
}

##' Bisulfite Conversion I Quality Control
##'
##' @title BS quality control plot
##' @description
##' `plotBS` assesses bisulfite conversion efficiency using control probes, 
##' identifying poorly converted samples. Type I BC control probes monitor the 
##' efficiency of the bisulfite conversion. If the conversion reaction was 
##' successful, the 'C' (converted) probes will be extended, but if the sample 
##' has unconverted DNA, the 'U' (unconverted) probes are extended instead. 
##' 
##' Performance of the BC control probes C1 and C2 is monitored in the green 
##' channel, and C3 and C4 are monitored in the red channel. The intensities 
##' for the relevant channels are then combined for all BC control probes and 
##' plotted per sample to ensure they are above the specified threshold.
##' 
##' @param object input sData
##' @param threshold threshold value
##' @param col col
##' @return plot
##' @author tkuipers, ljsinke
##' @export
##' @import ggplot2

plotBS <- function(object,
                   threshold,
                   col = NULL) {
  data <- object@plotdata
  
  d <- data[grepl(qcProbes["BSI"], data$Type), ]
  dGrn <- d[grepl("C1|C2|C3", d$ExtendedType), c(1:5,7)]
  x <- tapply(dGrn$IntGrn, dGrn$Samples, mean)
  dRed <- d[grepl("C4|C5|C6", d$ExtendedType), c(1:6)]
  y <- tapply(dRed$IntRed, dRed$Samples, mean)
  
  data <- data.frame(x, y)
  targets <- object@targets
  data <- merge(data, targets, by = "row.names", suffixes = c("", ".y"))
  
  data <- .rotateData(data, columns = c("x", "y"))
  
  ## Get outliers
  outliers <- data$Row.names[data$x <= threshold]
  
  .gcscatterplot(
    data,
    x = "x",
    y = "y",
    col = col,
    threshold_v = threshold,
    outliers = outliers,
    ylab = expression(paste(log[2], "(", R / G, ")")),
    xlab = expression(paste(log[2], sqrt(R %*% G))),
    main = "Bisulfite Conversion I quality control"
  )
}

##' Sample-independent overall quality control (Hyb)
##'
##' @title Sample-independent overall quality control (Hyb)
##' @description
##' `plotHC` uses hybridization controls to assess the hybridization step 
##' using synthetic targets, which complement the array perfectly. These 
##' synthetic targets are present in the hybridization buffer (RA1) at three 
##' concentration levels, and their performance is only monitored in the green 
##' channel. The difference in green intensity between the high (H) and low (L) 
##' concentration is combined and can be plotted to ensure it is above the 
##' specified threshold.
##' 
##' @param object input sData
##' @param threshold threshold value
##' @param col col
##' @return plot
##' @author tkuipers, ljsinke
##' @export
##' @import ggplot2

plotHC <- function(object,
                   threshold,
                   col = NULL) {
  data <- object@plotdata
  
  d <- data[grepl(qcProbes["HYB"], data$Type), ]
  d <- d[order(d$Samples), ]
  x <- 0.5 * (d$IntGrn[grepl("High", d$ExtendedType)] + d$IntGrn[grepl("Low", d$ExtendedType)])
  y <- d$IntGrn[grepl("High", d$ExtendedType)] - d$IntGrn[grepl("Low", d$ExtendedType)]
  
  data <- data.frame(x, y, row.names = d$Samples[grepl("High", d$ExtendedType)])
  targets <- object@targets
  data <- merge(data, targets, by = "row.names", suffixes = c("", ".y"))
  
  ## Get outliers
  outliers <- data$Row.names[data$x <= threshold]
  
  .gcscatterplot(
    data,
    x = "x",
    y = "y",
    col = col,
    threshold_v = threshold,
    outliers = outliers,
    ylab = expression(paste(log[2], "(", H / L, ")")),
    xlab = expression(paste(log[2], sqrt(H %*% L))),
    main = "Sample-independent overall quality control (Hyb)"
  )
}

##' Detection p-value based on NP probes
##'
##' @title Detection p-value based on NP probes
##' @description
##' This plot uses negative control probes, which are randomly permuted 
##' sequences that should not hybridize to the DNA template. The mean signal 
##' of these probes defines the background signal. This plot therefore show 
##' the fraction of probe visually distinct from the background signal in each 
##' sample, coloured by array number.
##' 
##' @param object input RGset
##' @param detP detection p-value threshold
##' @param threshold proportion of probes required
##' @param col col
##' @return plot
##' @author tkuipers, ljsinke
##' @export
##' @import ggplot2

plotDP <- function(object,
                   detP,
                   threshold,
                   col = NULL) {
  
  detP_matrix <- detectionP(object)
  y <- colMeans(detP_matrix < detP)
  x <- 1:length(y)
  data <- data.frame(x, y, row.names = names(y))
  targets <- as.data.frame(colData(RGset))
  data <- merge(data, targets, by = "row.names", suffixes = c("", ".y"))
  
  ## Get outliers
  outliers <- data$Row.names[data$y <= threshold]
  
  .gcscatterplot(
    data,
    x = "x",
    y = "y",
    col = col,
    threshold_h = threshold,
    outliers = outliers,
    xlab = "Samples",
    ylab = paste0("Proportion of probes with p-value < ", detP),
    main = "Detection p-value based on negative control probes"
  )
}

##' Get outliers
##'
##' @title Get outliers
##' @description
##' A simple function to extract outliers from targets.
##' 
##' @param object input sData
##' @param thresholds threshold values
##' @return outliers
##' @author tkuipers, ljsinke
##' @export

get_outliers <- function(object, thresholds) {
  outliers <- list()
  
  ## MU
  MU <- log2(t(object@MU))
  targets <- object@targets
  data <- merge(MU, targets, by = "row.names")
  outliers$MU <- data$Row.names[data$Methylated <= thresholds[1]]
  
  ## OP
  data <- object@plotdata
  d <- data[grepl(qcProbes["NP"], data$Type), ]
  dGrn <- d[d$ExtendedType %in% c("NP (C)", "NP (G)"), c(1:5, 7)]
  x <- tapply(dGrn$IntGrn, dGrn$Samples, mean)
  dRed <- d[d$ExtendedType %in% c("NP (A)", "NP (T)"), c(1:6)]
  y <- tapply(dRed$IntRed, dRed$Samples, mean)
  data <- data.frame(x, y)
  targets <- object@targets
  data <- merge(data, targets, by = "row.names", suffixes = c("", ".y"))
  data <- .rotateData(data, columns = c("x", "y"))
  outliers$OP <- data$Row.names[data$x <= thresholds[2]]
  
  ## BS
  data <- object@plotdata
  d <- data[grepl(qcProbes["BSI"], data$Type), ]
  dGrn <- d[grepl("C1|C2|C3", d$ExtendedType), c(1:5,7)]
  x <- tapply(dGrn$IntGrn, dGrn$Samples, mean)
  dRed <- d[grepl("C4|C5|C6", d$ExtendedType), c(1:6)]
  y <- tapply(dRed$IntRed, dRed$Samples, mean)
  data <- data.frame(x, y)
  targets <- object@targets
  data <- merge(data, targets, by = "row.names", suffixes = c("", ".y"))
  data <- .rotateData(data, columns = c("x", "y"))
  outliers$BS <- data$Row.names[data$x <= thresholds[3]]
  
  ## HC
  data <- object@plotdata
  d <- data[grepl(qcProbes["HYB"], data$Type), ]
  d <- d[order(d$Samples), ]
  x <- 0.5 * (d$IntGrn[grepl("High", d$ExtendedType)] + d$IntGrn[grepl("Low", d$ExtendedType)])
  y <- d$IntGrn[grepl("High", d$ExtendedType)] - d$IntGrn[grepl("Low", d$ExtendedType)]
  data <- data.frame(x, y, row.names = d$Samples[grepl("High", d$ExtendedType)])
  targets <- object@targets
  data <- merge(data, targets, by = "row.names", suffixes = c("", ".y"))
  outliers$HC <- data$Row.names[data$x <= thresholds[4]]
  
  ## DP
  y <- object@DPfreq
  x <- 1:length(y)
  data <- data.frame(x, y, row.names = names(y))
  targets <- object@targets
  data <- merge(data, targets, by = "row.names", suffixes = c("", ".y"))
  outliers$DP <- data$Row.names[data$y <= thresholds[5]]
  
  ## Clean output
  temp <- na.omit(plyr::ldply(outliers, data.frame))
  temp <- unstack(temp)
  
  all_res <- c("MU", "OP", "BS", "HC", "DP")
  
  outliers <- temp %>%
    rownames_to_column(var="sample") %>% 
    mutate(value = TRUE) %>%
    pivot_wider(names_from = res, values_from = value, values_fill = FALSE) %>%
    column_to_rownames("sample")
  
  missing_cols <- setdiff(all_res, colnames(outliers))
  outliers[missing_cols] <- FALSE
  
  outliers <- outliers[, all_res]
}

#functions to handle package installation and loading####
#this function check if the needed packages are installed, if not try to install from CRAN and if that fails try to install from bioconductor
install_missing_packages <- function(missing_packages) {
  if (length(missing_packages) == 0) {
    message("All required packages are already installed!")
    return(invisible(NULL))
  }
  
  # Ensure BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in missing_packages) {
    message(paste("Installing", pkg, "..."))
    
    # Try installing from CRAN
    install.packages(pkg, dependencies = TRUE)
    
    # Check if installation was successful
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("CRAN installation failed for", pkg, "- Trying Bioconductor..."))
      
      # Try installing from Bioconductor
      BiocManager::install(pkg, ask = FALSE)
      
      # Final check if still missing
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Failed to install", pkg, "from both CRAN and Bioconductor."))
      } else {
        message(paste("Successfully installed", pkg, "from Bioconductor."))
      }
    } else {
      message(paste("Successfully installed", pkg, "from CRAN."))
    }
  }
  
  message("Installation process completed.")
}

#this function checks if required packages are installed and return an error if not. Used inside functions.
check_required_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but is not installed. Please install it."))
    }
  }
}

#install and load required packages####
required_packages <- c("BiocManager", "ggplot2", "dplyr", "ggalign", "Seurat",
                       "hues", "phateR", "scran", "SingleCellExperiment", 
                       "tradeSeq", "ComplexHeatmap", "circlize", "irlba", 
                       "purrr", "tidyr", "tibble", "rlang", "ggrastr")

#load the required packages
invisible(sapply(required_packages, library, character.only = TRUE))

#general purpose custom functions####

#this is a simple function that returns the statistical mode of a vector.
statistical_mode <- function(x, na.rm = FALSE) {
  # Handle NA values based on the na.rm parameter
  if (!na.rm && any(is.na(x))) {
    return(NA)  # Return NA if there are NA values and na.rm is FALSE
  }
  
  # Remove NA values if na.rm is TRUE
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  
  # If the vector is empty after removing NAs, return NA
  if (length(x) == 0) {
    return(NA)
  }
  
  # Calculate the frequency of each unique value
  freq_table <- table(x)
  
  # Find the maximum frequency
  max_freq <- max(freq_table)
  
  # Extract the value(s) with the maximum frequency
  mode_values <- names(freq_table[freq_table == max_freq])
  
  # Convert to the appropriate type (numeric, logical, etc.)
  if (is.numeric(x)) {
    mode_values <- as.numeric(mode_values)
  } else if (is.logical(x)) {
    mode_values <- as.logical(mode_values)
  }
  
  return(mode_values)
}

#this is a simple function that will take a vector of numerical values and return the peak of the density distribution
density_peak <- function(x, ...){
  x <- density(x, ...)
  x <- data.frame(x = x$x, y = x$y)
  x <- x$x[which(x$y == max(x$y))]
  return(x)
}

#create colors
#this function takes a vector as input and returns a named vector with colors as values and name as elements from the input vector
#it uses the hues package to assign colors
create_colors <-function(x, 
                         palette = c("default", "all colors", "colorblind friendly", "fancy light", 
                                      "fancy dark", "shades", "tarnish", "pastel", "pimp", 
                                      "intense", "fluo", "red roses", "ochre sand", "yellow lime",
                                      "green mint", "ice cube", "blue ocean", "indigo night", "purple wine"), 
                         seed.use = 123){
  #get required libraries  
  check_required_packages("hues")
  #test inputs
  palette <- match.arg(palette)
  #create pallete parameters
  palette_param <- list(default = c(hmin = 0, hmax = 360, cmin = 30, cmax = 80, lmin = 35, lmax = 80),
                        `all colors` = c(hmin = 0, hmax = 360, cmin = 0, cmax = 100, lmin = 0, lmax = 100),
                        `colorblind friendly` = c(hmin = 0, hmax = 360, cmin = 40, cmax = 70, lmin = 15, lmax = 85),
                        `fancy light` = c(hmin = 0, hmax = 360, cmin = 15, cmax = 40, lmin = 70, lmax = 100),
                        `fancy dark` = c(hmin = 0, hmax = 360, cmin = 8, cmax = 40, lmin = 7, lmax = 40),
                        shades = c(hmin = 0, hmax = 240, cmin = 0, cmax = 15, lmin = 0, lmax = 100),
                        tarnish = c(hmin = 0, hmax = 360, cmin = 0, cmax = 15, lmin = 30, lmax = 70),
                        pastel = c(hmin = 0, hmax = 360, cmin = 0, cmax = 30, lmin = 70, lmax = 100),
                        pimp = c(hmin = 0, hmax = 360, cmin = 30, cmax = 100, lmin = 25, lmax = 70),
                        intense = c(hmin = 0, hmax = 360, cmin = 20, cmax = 100, lmin = 15, lmax = 80),
                        fluo = c(hmin = 0, hmax = 300, cmin = 35, cmax = 100, lmin = 75, lmax = 100),
                        `red roses` = c(hmin = 330, hmax = 20, cmin = 10, cmax = 100, lmin = 35, lmax = 100),
                        `ochre sand` = c(hmin = 20, hmax = 60, cmin = 20, cmax = 50, lmin = 35, lmax = 100),
                        `yellow lime` = c(hmin = 60, hmax = 90, cmin = 10, cmax = 100, lmin = 35, lmax = 100),
                        `green mint` = c(hmin = 90, hmax = 150, cmin = 10, cmax = 100, lmin = 35, lmax = 100),
                        `ice cube` = c(hmin = 150, hmax = 200, cmin = 0, cmax = 100, lmin = 35, lmax = 100),
                        `blue ocean` = c(hmin = 220, hmax = 260, cmin = 8, cmax = 80, lmin = 0, lmax = 50),
                        `indigo night` = c(hmin = 260, hmax = 290, cmin = 40, cmax = 100, lmin = 35, lmax = 100),
                        `purple wine` = c(hmin = 290, hmax = 330, cmin = 0, cmax = 100, lmin = 0, lmax = 40))
  
  if(length(x) != length(unique(x))){
    warning("Input vector contain multiple identical values. Colapsing the vector to unique values only...")
    x <- unique(x)
  }
  
  set.seed(seed.use)
  colors <- hues::iwanthue(n = length(x), 
                     hmin = palette_param[[palette]]["hmin"], 
                     hmax = palette_param[[palette]]["hmax"], 
                     lmin = palette_param[[palette]]["lmin"], 
                     lmax = palette_param[[palette]]["lmax"], 
                     cmin = palette_param[[palette]]["cmin"], 
                     cmax = palette_param[[palette]]["cmax"])
  x <- setNames(colors, x)
  return(x)
}

#functions to plot scRNA data from seurat objects####
  
#this function plots features in seurat Objects using custom names
#currently it doesn't support split.by
FeaturePlot_with_names <- function(seuratObj, features, name_by, reduction = "umap", assay = NULL, layer = "data", na_cuttoff = 1e-09, na_color = "lightgray", colors_use = NULL, pt.size = NULL){
  required_packages <- c("Seurat", "tibble", "tidyr", "dplyr", "cowplot")
  check_required_packages(required_packages)
  
  if(is.null(colors_use)){
    colors_use <- c("#fff9cc", "#fde954", "#ea5436", "#445ea1", "#20294f")
  }
  
  embeddings <- Embeddings(seuratObj, reduction = reduction)
  colnames(embeddings)[1:2] <- c("dim1", "dim2")
  
  if(is.null(assay)){
    assay <- DefaultAssay(seuratObj)
  }
  
  if(is.null(pt.size)){
    pt.size <- ifelse(ncol(seuratObj) < 2000, 1, 0.5)
  }
  
  reduction_names <- colnames(Embeddings(seuratObj, reduction = reduction))[1:2]
  
  to_plot <- FetchData(seuratObj, slot = slot, layer = layer, assay = assay, vars = c(features, reduction_names))%>%
    rownames_to_column("cell") %>%
    pivot_longer(cols = any_of(unname(features)), names_to = "gene_id", values_to = "expression") %>%
    mutate(gene_name = seuratObj[["RNA"]]@meta.data[gene_id,name_by]) %>%
    mutate(gene_name = substr(gene_name, 1, 40)) %>%
    mutate(gene_name = ifelse(nchar(gene_name) == 40, paste0(gene_name, "..."), gene_name))
  
  #plot the figure
  plot_list <- list()
  for(gene in features){
    gene_name <- unique(filter(to_plot, gene_id == gene)$gene_name)
    plot_list[[gene]] <- to_plot %>%
      filter(gene_id == gene) %>%
      arrange(expression) %>%
      ggplot(aes(x = .data[[reduction_names[1]]], 
                                   y = .data[[reduction_names[2]]], 
                                   color = expression)) +
      geom_point(size = pt.size) +
      labs(color = NULL, title = gene_name, subtitle = gene, 
                    x = paste0(reduction, "_1"), y = paste0(reduction, "_2")) +
      scale_color_gradientn(colors = colors_use, limits = c(na_cuttoff, NA), na.value = na_color) +
      theme_cowplot()
    
  }
  
  #arrange by alphabetic order
  plot_list <- plot_list[names(sort(sapply(plot_list, function(x) x$labels$title)))]
  
  #calculate optimal plot width and height
  plot_width <- round(sqrt(length(plot_list)))
  plot_height <- floor(length(plot_list)/plot_width)
  
  cat(paste0("Recomended dimensions:\nwidth: ", plot_width*6, " inches\nheight: ", plot_height*5, " inches"))
  
  return(plot_grid(plotlist = plot_list, ncol = plot_width))
}

flag_outlier_cells <- function(seuratObj, qc_vars, out_method = c("quantile_up", "quantile_down", "boxplot"), quantile_threshold = 0.99){
  
  #validate inputs
  out_method <- match.arg(out_method)

  missing_vars <- qc_vars[!qc_vars %in% colnames(seuratObj@meta.data)]
  if(length(missing_vars) != 0){
    stop(paste0("The following features are not present in the object metadata: ", paste0(missing_vars, collapse = ", ")))
  }
  
  if(!is.numeric(quantile_threshold)){
    stop("`quantile_threshold` must be numeric.")
  }
  
  #prepare data for outlier flagging
  out_flags <- seuratObj@meta.data %>%
    rownames_to_column("cell_barcode") %>%
    select(any_of(qc_vars), cell_barcode) %>%
    pivot_longer(cols = qc_vars, names_to = "feature") %>%
    group_by(feature)
  
  #flag outliers based on out_method
  if(out_method == "quantile_up"){
    out_flags <- mutate(out_flags, is_outlier = value > quantile(value, quantile_threshold))
  }
  
  if(out_method == "quantile_down"){
    out_flags <- mutate(out_flags, is_outlier = value <= quantile(value, quantile_threshold))
  }
  
  if(out_method == "boxplot"){
    out_flags <- mutate(out_flags, is_outlier = value %in% boxplot.stats(value)$out)
  }
  
  
  #pivot it back to wide format
  out_flags <- out_flags %>%
    select(-value) %>%
    pivot_wider(names_from = feature, values_from = is_outlier, names_prefix = "is_outlier_")
  
  #append the flags to the metadata
  seuratObj@meta.data <- seuratObj@meta.data %>%
    rownames_to_column("cell_barcode") %>%
    select(-any_of(setdiff(colnames(out_flags), "cell_barcode"))) %>%
    left_join(out_flags, by = "cell_barcode") %>%
    column_to_rownames("cell_barcode")
  
  return(seuratObj)
}


#FeatLabels
#This function replace feature labels in seurat plot
#labels should be a named vector where names are the gene IDs and values are the labels to replace the IDs. 
FeatLabels <- function(plot, labels, labels_max_char = 40){
  required_packages <- c("dplyr", "forcats", "ggplot2")
  check_required_packages(required_packages)
    if(is.numeric(labels_max_char)){
      labels <- labels %>%
        substr(start = 1, stop = labels_max_char) %>%
        ifelse(nchar(.) == labels_max_char, paste0(., "..."), .)
    }
    
    if(inherits(plot, "ggplot")){
      plot <- list(plot)
    }
  
    for(i in c(1:length(plot))){
      plot_data <- plot[[i]]$data
      column <- grep("feature|gene", tolower(colnames(plot_data)))
      plot_data$feature <- plot_data[,column]
      plot_data <- plot_data %>%
        mutate(new_label = labels[as.character(feature)]) %>%
        mutate(new_label = fct_reorder(new_label, as.numeric(feature)))
      plot[[i]]$data[,column] <- plot_data$new_label
    }
  if(length(plot) == 1){
    plot <- plot[[1]]
  }
  return(plot)
}

#Functions to perform custom analysis in Seurat Objects####
##RunPHATE####
#this function run PHATE dimensionality reduction in a seurat object. 
#It basically wraps the phate() function from phater and stores the result in the correct slot of the seurat Object.
library(phateR)
library(Seurat)
library(dplyr)

RunPHATE <- function(seuratObj, dims = NULL, features = NULL, layer = "counts", assay = "RNA", reduction.key = "PHATE_", ...){
  required_packages <- c("Seurat", "phateR")
  check_required_packages(required_packages)
  if(is.null(dims)){
    if(is.null(features)){
    features <- Features(seuratObj, assay = assay, layer = layer)
    }
    phate_map <- t(GetAssayData(seuratObj, assay = assay, layer = layer)[features,])
    }else{
    phate_map <- Embeddings(seuratObj, reduction = "pca")[,dims]
  }
  
 
  phate_map <- phate(phate_map, ...)
  seuratObj@reductions$phate <- CreateDimReducObject(embeddings = 100*as.matrix(phate_map),
                                                     assay = assay,
                                                     key = reduction.key)
  return(seuratObj)
}

##scran normalize####
scran_normalize <- function(x, assay = "RNA"){
  required_packages <- c("SeuratObject", "SingleCellExperiment", "SummarizedExperiment", "scran", "scuttle")
  check_required_packages(required_packages)
  
  # Restrict to raw-count layers; excludes existing data.*, scale.data, etc.
  count_layers <- grep(
    "^counts(\\.|$)",
    SeuratObject::Layers(x[[assay]]),
    value = TRUE
  )
  
  if(!length(count_layers) > 0){
   stop("There are no count layers in the object") 
  }else{
    message(paste("Found", length(count_layers), "layers in the", assay, "assay. Normalizing each individually..."))
  }
  
  for (count_layer in count_layers) {
    counts <- SeuratObject::LayerData(
      object = x,
      assay = assay,
      layer = count_layer
    )
  
    # One SCE per sample/layer
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = counts)
    )
    
    sce <- scuttle::computePooledFactors(sce)
    sce <- scuttle::logNormCounts(sce)
    
    # counts.sampleA -> data.sampleA; counts -> data
    data_layer <- sub("^counts", "data", count_layer)
    
    SeuratObject::LayerData(
      object = x,
      assay = assay,
      layer = data_layer
    ) <- SummarizedExperiment::assay(sce, "logcounts")
  }
  return(x)
}

#plot_smoothers#####
#this function is an alternative to the PlotSmoothers() function from Tradeseq. It allows multiple genes as well as to split the plots by lineages or conditions.
#It also allow to display only specific lineages.
plot_smoothers <- function(models, counts, genes, lineages = NULL, gene_names = NULL, gene_names_max_char = 40, split.by = c("none", "condition", "lineage"), nPoints = 100, colors = NULL, show_knots = FALSE){
  required_packages <- c("tradeSeq", "ggplot2", "tibble", "dplyr", "hues", "rlang", "tidyr")
  check_required_packages(required_packages)
  
  split.by <- match.arg(split.by)
  #if counts are not provided, try to extract it from the SingleCellExperiment object
  if(missing(counts)){
    if(class(models) == "SingleCellExperiment"){
      counts <- models@assays@data@listData$counts[genes, , drop = FALSE]
    }else{
      stop("counts missing, and can't get it from models")
    }
  }
  
  #determine if there are conditions
  if("conditions" %in% names(models@colData$tradeSeq)){
    condition <- as.character(models@colData$tradeSeq$conditions)
  }else{
    condition <- rep(NA, times = ncol(counts))
  }
  
  #get smoothers
  smoothers <- predictSmooth(models, gene = genes, nPoints = nPoints, tidy = TRUE) %>%
    mutate(count = yhat)
  
  if(is.numeric(lineages)){
    smoothers <- filter(smoothers, lineage %in% lineages)
  }
  if(!"condition" %in% colnames(smoothers)){
    smoothers$condition <- NA
  }
  colnames(smoothers) <- gsub("time", "pseudotime", colnames(smoothers))
  
  #get knots
  knots <- models@metadata$tradeSeq$knots
  
  #get cell counts 
  cells <- counts %>%
    t() %>%
    as.data.frame() %>%
    cbind(., condition) %>%
    rownames_to_column("cell") %>%
    pivot_longer(cols = -c("cell", "condition"), values_to = "count", names_to = "gene")
  
  #get to which lineage each cell belongs
  lineage_data <- models@colData@listData$crv %>%
    as.data.frame() %>%
    select(contains("cellWeights")) 
  lineage_data$lineage <- apply(lineage_data, 1, function(x) x <- which(x == max(x))[1])
  lineage_pseudotime <- models@colData@listData$crv %>%
    as.data.frame() %>%
    select(contains("pseudotime.Lineage"))
  lineage_data$pseudotime <- NA
  for(i in c(1:nrow(lineage_data))){
    lineage_data$pseudotime[i] <- lineage_pseudotime[i,lineage_data$lineage[i]]
  }
  cells <- cbind(cells, lineage_data[cells$cell,])
  rownames(cells) <- NULL
  
  if(is.numeric(lineages)){
    cells <- filter(cells, lineage %in% lineages)
  }
  
  #set which variable will be used for coloring
  split_vs_color <- c(none = "group", condition = "lineage", lineage = "condition")[split.by]
  #add a column combining conditions and lineages and indicate which column will be used to split
  cells <- cells %>%
    mutate(lineage = paste0("Lineage ", lineage)) %>%
    mutate(group = ifelse(is.na(condition), lineage, paste0(lineage, "_", condition))) %>%
    mutate(none = NA) %>%
    mutate(split = !!data_sym(split.by)) %>%
    mutate(color = !!data_sym(split_vs_color))
  
  smoothers <- smoothers %>%
    mutate(lineage = paste0("Lineage ", lineage)) %>%
    mutate(group = ifelse(is.na(condition), lineage, paste0(lineage, "_", condition))) %>%
    mutate(none = NA) %>%
    mutate(split = !!data_sym(split.by)) %>%
    mutate(color = !!data_sym(split_vs_color))
  
  if(is.null(colors)){
    set.seed(123)
    colors <- iwanthue(length(unique(cells$color)), hmin = 0, hmax = 360, cmin = 40, cmax = 70, lmin = 15, lmax = 85)  
    names(colors) <- unique(cells$color)  
  }
  
  #change gene names if requested
  if(is.null(gene_names)){
    gene_names <- setNames(genes, genes)
    gene_labeller <- function(x){return(as.character(x))}
  }else{
    gene_labeller <- function(x){
      x <- as.character(x)
      label <- gene_names[x]
      label <- as.character(label)
      label <- paste(label,x , sep = "\n")
      return(label)
    }
  }
  #plot    
  plot <- cells %>%
    ggplot(aes(x = pseudotime, y = log1p(count), color = factor(color), group = factor(color)))+
    geom_point(size = 0.1)+
    geom_path(data = smoothers, linewidth = 1, color = "white")+
    geom_path(data = smoothers, linewidth = 0.5)+
    {if(show_knots)geom_vline(xintercept = knots, linetype = "dashed", color = "lightgrey")}+
    scale_color_manual(values = colors)+
    labs(color = split_vs_color, y = "Log(Count + 1)", x = "Pseudotime")+
    {if(split.by == "none")facet_wrap(vars(factor(gene, levels = genes)), scales = "free", labeller = as_labeller(gene_labeller))}+
    {if(split.by != "none")facet_grid(rows = vars(factor(gene, levels = genes)), cols = vars(split), scales = "free", labeller = as_labeller(gene_labeller))}+
    theme_classic()+
    theme(strip.text = element_text(hjust = 0, margin = margin(l = 3, r = 3, t = 2, b = 0), size = 7))
  return(plot)  
}

#meta_from_ref
#this function will project a query seurat object into a reference seurat object to extract metadata from the reference into the query.
#it allows the use of an `idconversion` named vector to convert gene ids into orthologs for instance. This way reference and query objects generated with different genomes can be used. 
meta_from_ref <- function(query_obj, ref_obj, normalization_method = c("LogNormalize", "SCT"), dims = 1:20, metadata_to_extract, ignore_values = NULL, id_conversion = NULL, subsample_ref = NULL, prefix = "predicted.", vars_to_regress = NULL){
  #libraries
  required_packages <- c("Seurat", "dplyr", "purrr", "tidyr", "tibble")
  check_required_packages(required_packages)
  
  #test if the required metadata are present in the reference dataset
  if(sum(!metadata_to_extract %in% colnames(ref_obj@meta.data)) != 0){
    missing_meta <- metadata_to_extract[which(!metadata_to_extract %in% colnames(ref_obj@meta.data))]
    stop(paste0("Could not find the \'", missing_meta, "\' metadata in the reference object."))
  }
  
  #remove cells with values to be ignored from the reference
  if(!is.null(ignore_values)){
    ncells_before <- ncol(ref_obj)
    for(i in c(1:length(ignore_values))){
      metadata <- names(ignore_values)[i]
      values_to_ignore <- ignore_values[[i]]
      if(!metadata %in% colnames(ref_obj@meta.data)){
        stop(paste0("could not find '", metadata, "' in the reference metadata"))
      }
      if(sum(!values_to_ignore %in% ref_obj@meta.data[,metadata]) > 0){
        missing_values <- values_to_ignore[!which(values_to_ignore %in% ref_obj@meta.data[,metadata])]
        stop(paste0("could not find the variable '", paste(missing_values, collapse = ", "), "' in the '", metadata, "' column in the reference metadata"))
      }
      print(paste0("Removing cells with '", metadata, "' equal to '", paste(values_to_ignore, collapse = ", "), "'"))
      cells_to_keep <- rownames(ref_obj@meta.data[which(!ref_obj@meta.data[,metadata] %in% values_to_ignore), , drop = FALSE])
      ref_obj <- subset(ref_obj, cells = cells_to_keep)
    }
    ncells_after <- ncol(ref_obj)
    print(paste0("Removed a total of ", ncells_before-ncells_after, " cells out of ", ncells_before, " cells from the reference"))
  }
  
  #subsample reference
  if(is.numeric(subsample_ref)){
    print(paste("reference object is being sampled for a total of", subsample_ref, "cells"))
    ref_obj <- subset(ref_obj, cells = sample(Cells(ref_obj), subsample_ref))
  }
  
  #determine normalization method
  normalization_method <- match.arg(normalization_method)
  
  #save the current cell metadata 
  query_meta <- query_obj@meta.data
  #create a query count matrix with reference IDs
  query_counts <- query_obj[["RNA"]]$counts #get count matrix
  
  #simplify the reference object to save memory
  ref_obj <- CreateSeuratObject(LayerData(ref_obj, assay = "RNA", layer = "counts"), meta.data = ref_obj@meta.data)
  
  if(!is.null(id_conversion)){
    print("converting IDs of the query object...")
    ref_total <- rownames(ref_obj[["RNA"]]$counts)
    ref_not_in_conversion <- sum(!ref_total %in% id_conversion)
    query_total <- rownames(query_counts)
    query_not_in_conversion <- sum(!query_total %in% names(id_conversion))
    if(ref_not_in_conversion > 0){
      warning(paste0(ref_not_in_conversion, " features out of ", length(ref_total), " reference features were not present in the `id_conversion` vector."))
    }
    if(query_not_in_conversion > 0){
      warning(paste0(query_not_in_conversion, " features out of ", length(query_total), " query features were not present in the `id_conversion` vector."))
    }
    query_counts <- query_counts[rownames(query_counts) %in% names(id_conversion),]
    rownames(query_counts) <- id_conversion[rownames(query_counts)]
  }
  #at this point we should have a query count matrix with feature IDs matching the IDs of the reference.
  #Now we create a new seurat object just with this count matrix
  new_query_obj <- CreateSeuratObject(query_counts, meta.data = query_meta)
  
  if(normalization_method == "SCT"){
    #run SCTransform on both objects
    print("Running SCT transformation on both objects...")
    new_query_obj <- SCTransform(new_query_obj, return.only.var.genes = FALSE, verbose = FALSE, vars.to.regress = vars_to_regress)
    ref_obj <- SCTransform(ref_obj, return.only.var.genes = FALSE, verbose = FALSE, vars.to.regress = vars_to_regress)
  }else{
    print("Running LogNormalize on both objects")
    new_query_obj <- new_query_obj %>%
      NormalizeData() %>%
      ScaleData(vars.to.regress = vars_to_regress) %>%
      FindVariableFeatures()
    ref_obj <- ref_obj %>%
      NormalizeData() %>%
      ScaleData(vars.to.regress = vars_to_regress) %>%
      FindVariableFeatures()
  }
  print("Running PCA on reference object...")
  ref_obj <- RunPCA(ref_obj, feature = VariableFeatures(ref_obj), verbose = FALSE)
  print("Running PCA on query object...")
  ref_obj <- RunPCA(ref_obj, feature = VariableFeatures(ref_obj), verbose = FALSE)
  
  #Finding Anchors
  print("Finding Transfer Anchor features...")
  Anchors <- FindTransferAnchors(query = new_query_obj, reference = ref_obj, normalization.method = normalization_method, reference.reduction = "pca", verbose = FALSE, dims = dims)
  
  #Map query to reference and retrieve metadata. 
  print("Mapping query to reference and retrieving metadata...")
  for(meta in metadata_to_extract){
    ref_meta <- MapQuery(anchorset = Anchors, reference = ref_obj, query = query_obj, refdata = meta, verbose = FALSE)@meta.data[,c("predicted.id.score", "predicted.id"), drop = FALSE]
    colnames(ref_meta) <- gsub("id", meta, colnames(ref_meta))
    colnames(ref_meta) <- gsub("predicted.", prefix, colnames(ref_meta))
    query_obj@meta.data[colnames(ref_meta)] <- ref_meta
  }
  return(query_obj)
}

#DoHeatmap_plus####
#function that makes more advanced heatmaps for seurat objects
DoHeatmap_plus <- function(seuratObj, 
                           features, 
                           cells = NULL, 
                           group.by = "ident", 
                           layer = "scale.data", 
                           assay = NULL, 
                           show_feature_labels = NULL, 
                           feature_labels = NULL, 
                           feature_labels_max_char = 40, 
                           feature_label_size = 11,
                           group.colors = NULL,  
                           hclust_method = "complete", 
                           expression_colors = NULL, 
                           raster = NULL, 
                           raster_resolution = 300){
  
  required_packages <- c("Seurat", "ggalign", "dplyr", "ggrastr")
  check_required_packages(required_packages)
  #set defaults
  if(is.null(cells)){
    cells <- colnames(seuratObj)
  }
  if(is.null(assay)){
    assay <- DefaultAssay(seuratObj)
  }
  if("ident" %in% group.by){
    seuratObj$ident <- Idents(seuratObj)
  }
  if(is.null(show_feature_labels)){
    if(!is.null(feature_labels)){
      show_feature_labels <- TRUE
    }else{
      show_feature_labels <- length(features) <= 100  
    }
  }
  
  #get the matrix
  heat_matrix <- as.matrix(GetAssayData(seuratObj, assay = assay, layer = layer))
  #test if the features are in the provided matrix
  missing_features <- features[which(!features %in% rownames(heat_matrix))]                         
  if(length(missing_features) != 0){
    warning(paste(length(missing_features), "features are not found in the matrix in assay:", assay, "and layer:", layer))
    features <- features[which(features %in% rownames(heat_matrix))]
  }
  #test if the cells are in the provided matrix
  missing_cells <- features[which(!cells %in% colnames(heat_matrix))]                         
  if(length(missing_cells) != 0){
    warning(paste(length(missing_cells), "cells are not found in the matrix in assay:", assay, "and layer:", layer))
    cells <- cells[which(cells %in% colnames(heat_matrix))]
  }  
  #reverse the order of features so heatmap looks correct
  features <- rev(features)
  
  #subset the matrix
  heat_matrix <- heat_matrix[features, cells, drop = FALSE]
  
  if(is.null(raster)){
    raster <- max(dim(heat_matrix)) >= 2000
  }
  
  #cap the expression matrix to quantiles to remove outliers
  limits <- quantile(heat_matrix, c(0.05, 0.95))
  heat_matrix[heat_matrix < limits[1]] <- limits[1]
  heat_matrix[heat_matrix > limits[2]] <- limits[2]
  
  #create color scale for the matrix
  if(is.null(expression_colors)){
    expression_colors <- c("#cc00ff", "black", "#fbff00")
  }
  intervals <- c(limits[[1]], 0, limits[[2]])
  breaks <- c(ceiling(limits[[1]]), floor(limits[[2]]))
  breaks <- breaks[1]:breaks[2]
  labels <- as.character(breaks)
  
  expression_colors <- scale_fill_gradientn(name = layer, 
                                            limits = limits, 
                                            colors = expression_colors, 
                                            breaks = breaks,
                                            labels = labels,
                                            values = scales::rescale(intervals))
  
  #cluster the cells
  if(!is.null(group.by)){
    group_data <- seuratObj@meta.data[colnames(heat_matrix), group.by, drop = FALSE]
    message(paste("clustering", ncol(heat_matrix), "cells based on", nrow(heat_matrix), "features using the method:", hclust_method))
    cells_cluster <- group_data %>%
      rownames_to_column("cell") %>%
      group_by(!!!syms(c(group.by))) %>%
      mutate(order = ifelse(n() >1, hclust(dist(t(heat_matrix[,cell])), method = hclust_method)$order, 1)) %>%
      arrange(!!!syms(c(group.by, "order"))) %>%
      mutate(new_order = row_number()) %>%
      ungroup()
    hclust_method <- NULL
    #determine cells order
    cells_order <- cells_cluster %>%
      select(cell, new_order) %>%
      deframe()  
    #determine index where columns should be break
    col_breaks <- cells_cluster %>%
      pull(group.by[1])
    #reorder the matrix
    heat_matrix <- heat_matrix[,names(cells_order)]
  }else{
    col_breaks <- NULL
    cells_annotation <- NULL
  }
  #replace gene names by labels
  if(show_feature_labels){
    if(!is.null(feature_labels)){
      feature_labels <- ifelse(nchar(feature_labels) > feature_labels_max_char, str_sub(feature_labels, 1, feature_labels_max_char), feature_labels)
      feature_labels <- ifelse(nchar(feature_labels) == feature_labels_max_char, paste0(feature_labels, "..."), feature_labels)
    }else{
      feature_labels <- setNames(rownames(heat_matrix), rownames(heat_matrix))
    }
    feature_labels <- scale_y_continuous(name = NULL, breaks = rownames(heat_matrix), labels = feature_labels[rownames(heat_matrix)])
  }else{
    feature_labels <- scale_y_continuous(name = NULL, breaks = NULL, labels = NULL)
  }
  
  #prepare geom tiles according to `raster` parameter
  if(raster){
    message(paste("Rasterizing heatmap with dpi:", raster_resolution))
    heat_tile <- ggrastr::rasterise(geom_tile(aes(fill = value)), dpi = raster_resolution)
    anno_tile <- ggrastr::rasterise(geom_tile(aes(y = 1, fill = factor(value))), dpi = raster_resolution)
  }else{
    message("Generating heatmap without rastering.")
    heat_tile <- geom_tile(aes(fill = value))
    anno_tile <- geom_tile(aes(y = 1, fill = factor(value)))
  }
  
  #make base plot
  plot <- ggheatmap(heat_matrix, filling = FALSE)+
    scale_x_continuous(labels = NULL, breaks = NULL)+
    expression_colors+
    feature_labels+
    theme(axis.text.y = element_text(size = feature_label_size))+
    heat_tile

  #add grouping annotations
    for(group_var in group.by){
      annotation <- setNames(group_data[[group_var]], rownames(group_data)) 
      if(is.factor(annotation)){
        annotation <- unfactor(annotation)
      }
      
      if(length(unique(annotation)) == 1){
        warning(paste0(group_var, " has only 1 value so it was skipped..."))
        group.by <- group.by[-which(group.by == group_var)] 
        next
      }
      #use the first group_var to also split the heatmap
      if(group_var == group.by[1]){
        plot <- plot + 
          anno_top() + 
          align_group(annotation[colnames(heat_matrix)]) 
      }
      plot <- plot + anno_top(size = 0.05) +
        ggalign(data = annotation[colnames(heat_matrix)])+
        anno_tile+
        scale_y_continuous(name = group_var)
      if(group_var %in% names(group.colors)){
        plot <- plot + scale_fill_manual(name = group_var, values = group.colors[[group_var]])
      }else{
        message(paste0("Colors not provided for ", group_var, ". Using standard ggplot colors."))
        plot <- plot + scale_fill_discrete(name = group_var)
      }
        plot <- plot +
        theme_void() +
        theme(axis.title.y = element_text(),
                strip.text = element_blank())
    }
  
    
  return(plot)
}

#annotate_cells
#this function annotate cells based on the expression of gene modules.
#gene modules should be provided as a named list in `feature_list` argument. 
annotate_cells <- function(seuratObj, feature_list, annotation_name = "annotation", unassigned_name = "not_annotated", method = c("ModuleScore", "MetaFeature"), threshold = 0, assay = NULL, layer = "data", ...){
  required_packages <- c("Seurat", "dplyr", "tidyr", "tibble")
  check_required_packages(required_packages)
  
  #check inputs
  method <- match.arg(method)
  if(is.null(assay)){
    assay <- DefaultAssay(seuratObj)
  }
  #check if there are missing genes in the data
  missing_genes <- unique(unlist(feature_list))
  missing_genes <- missing_genes[!which(missing_genes %in% rownames(seuratObj[[assay]]))]
  if(length(missing_genes) != 0){
    stop(paste("Could not find the following genes in the", assay, "assay:"), paste0(missing_genes, collapse = "; "))
  }
  
  #remove pre-existing columns from metadata
  cols_to_remove <- c(names(feature_list), annotation_name, paste0(annotation_name, "_score"))
  if(sum(cols_to_remove %in% colnames(seuratObj@meta.data)) > 0){
    warning(paste("Columns", paste(cols_to_remove, collapse = ", "), "were already present in the meta data and will be replaced by new values."))
    cols_to_remove <- paste0(cols_to_remove, collapse = "|")
    seuratObj@meta.data <- seuratObj@meta.data[,-grep(cols_to_remove, colnames(seuratObj@meta.data))]
  }
  rm(cols_to_remove)
  
  #assing gene module scores
  if(method == "ModuleScore"){
    #this assigns a score to each cell for each gene module using the AddModuleScore() function
    seuratObj <- AddModuleScore(seuratObj, features = feature_list, name = names(feature_list), assay = assay, slot = layer, ...)
    
    #the AddModuleScore function create weird names. Fix it.
    colnames(seuratObj@meta.data)[match(paste0(names(feature_list), c(1:length(feature_list))), colnames(seuratObj@meta.data))] <- names(feature_list)
  }
  
  if(method == "MetaFeature"){
    #this assigns a score to each cell for each gene module using the MetaFeature() function
    for(annotation in names(feature_list)){
      seuratObj <- MetaFeature(seuratObj, features = feature_list[[annotation]], assay = assay, slot = layer, meta.name = annotation, ...)  
    }
    seuratObj@meta.data[,names(feature_list)] <-  seuratObj@meta.data[,names(feature_list)] / colMeans(seuratObj@meta.data[,names(feature_list)])
  }
  
  annotation_data <- seuratObj@meta.data[,names(feature_list)] %>%
    rownames_to_column("cell") %>%
    pivot_longer(cols = -cell, names_to = "annotation", values_to = "score") %>%
    group_by(cell) %>%
    summarise(max_score = max(score), annotation = annotation[which(score == max_score)]) %>%
    mutate(annotation = ifelse(max_score < threshold, unassigned_name, annotation)) %>%
    mutate(annotation = factor(annotation, levels = c(names(feature_list), unassigned_name))) %>%
    select(cell, annotation, max_score) %>%
    dplyr::rename(!!annotation_name := annotation, !!paste0(annotation_name, "_score") := max_score) %>%
    column_to_rownames("cell") 
  
  seuratObj@meta.data <- cbind(seuratObj@meta.data, annotation_data[rownames(seuratObj@meta.data),]) 
  return(seuratObj)
}


#RunPCA_plus
#this function runs PCA on a seurat object but allow more flexibility in the PCA parameters, such as wether features are scaled or not for instance
RunPCA_plus <- function(seuratObj, assay = NULL, features = NULL, layer = "scale.data", scale = TRUE, npcs = 50, reduction.name = "pca", reduction.key = "PC_", weight.by.var = TRUE, seed.use = 42){
  required_packages <- c("Seurat", "irlba")
  check_required_packages(required_packages)
  
  #set defaults
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if(is.null(assay)){
    assay <- DefaultAssay(seuratObj)
  }
  if(is.null(features)){
    features <- VariableFeatures(seuratObj, assay = assay)
  }
  
  #get the counts
  data_matrix <- as.matrix(FetchData(seuratObj, features, assay = assay, layer = layer))
  #scale it
  if(scale & !grepl("scale", layer)){
    print("scaling the matrix...")
    data_matrix <- scale(data_matrix)
  }
  pca.results <- irlba(data_matrix, nv = npcs)
  
  feature.loadings <- pca.results$v
  sdev <- pca.results$d/sqrt(max(1, nrow(data_matrix) - 1))
  if (weight.by.var) {
    cell.embeddings <- pca.results$u %*% diag(pca.results$d)
  }
  else {
    cell.embeddings <- pca.results$u
  }
  rownames(feature.loadings) <- colnames(data_matrix)
  colnames(feature.loadings) <- paste0(reduction.key, 1:npcs)
  rownames(cell.embeddings) <- rownames(data_matrix)
  colnames(cell.embeddings) <- colnames(feature.loadings)
  
  seuratObj@reductions[[reduction.name]] <- CreateDimReducObject(embeddings = cell.embeddings, 
                       loadings = feature.loadings, assay = assay, stdev = sdev, 
                       key = reduction.key)
   
  return(seuratObj)
}

VariableFeaturePlot_SCT <- function(seuratObj, assay = "SCT", nFeat_label = 20, label_info = "feature_id", label_trim = 30){
  
  #test if assay is present
  if(!assay %in% names(seuratObj@assays)){
    stop(paste0("Cannot find assay `", assay, "` in the provided seurat object."))
  }
  
  #prepare data for plotting
  plot_list <- list()
  for(i in c(1:length(seuratObj[[assay]]@SCTModel.list))){
    to_plot <- seuratObj[[assay]]@SCTModel.list[[i]]@feature.attributes %>%
      as.data.frame() %>%
      rownames_to_column("feature") %>%
      left_join(
        seuratObj[["RNA"]]@meta.data %>% #get metadata from "RNA" assay
          rownames_to_column("feature") %>%
          mutate(feature_id = feature) %>%
          select(feature, any_of(label_info))
        , by = "feature"
      ) %>%
      mutate(is_variable = feature %in% VariableFeatures(seuratObj)) %>%
      mutate(label = ifelse(feature %in% VariableFeatures(seuratObj)[1:nFeat_label], cellranger_name, NA)) %>% #get labels for top genes
      mutate(label = ifelse(nchar(label) > label_trim, paste0(strtrim(label, label_trim), "..."), label))
    
    color_labels <- to_plot %>%
      group_by(is_variable) %>%
      summarise(count = n()) %>%
      deframe()
    
    color_labels["FALSE"] <- paste("Non-variable count:", color_labels["FALSE"])
    color_labels["TRUE"] <- paste("Variable count:", color_labels["TRUE"])
    
    plot_list[[i]] <- ggplot(to_plot, aes(x = gmean, y = residual_variance, color = is_variable))+
      geom_point(size = 0.8)+
      geom_text_repel(data = filter(to_plot, !is.na(label)), aes(label = label), show.legend = FALSE, min.segment.length = 0, color = "black")+
      scale_x_continuous(transform = "log10")+
      theme_cowplot()+
      ggtitle(paste0("SCT model for library: ", levels(factor(seuratObj$orig.ident))[i]))+
      scale_color_manual(values = c("black", "red"), labels = color_labels, name = NULL)+
      theme(legend.position = "top")
  }
  final_plot <- ggalign::align_plots(!!!plot_list)
  return(final_plot)
}

#plot_table
## this function plots a table using the geom_tile, returning a ggplot object that can be alligned with other plots with ggaling
plot_table <- function(x){
  #check packages
  required_packages <- c("dplyr", "tibble", "ggplot2")
  check_required_packages(required_packages)
  
  #check if table is a dataframe
  if(!is.data.frame(x)){
    stop("object must be a data frame")
  }
  
  #prepare table header
  header <- tibble(
    row = ".header",
    column = factor(colnames(x), levels = colnames(x)),
    value = colnames(x),
    header = TRUE
  )
  
  x <- x %>%
    rownames_to_column("row") %>%
    mutate(row = factor(row, levels = rownames(x))) %>%
    mutate(across(-row, as.character)) %>%
    pivot_longer(-row, names_to = "column", values_to = "value") %>%
    mutate(
      column = factor(column, levels = levels(header$column)),
      header = FALSE
    ) %>%
    bind_rows(header) %>%
    mutate(
      row = factor(row, levels = c(rev(unique(row[row != ".header"])), ".header"))
    ) %>%
    group_by(column) %>%
    mutate(col_width = nchar(as.character(column))) %>%
    mutate(col_width = max(col_width)) %>%
    ggplot(aes(x = column, y = row, label = value)) +
    geom_tile(fill = "white", colour = "black") +
    geom_text(aes(fontface = if_else(header, "bold", "plain")), size = 3) +
    theme_void()
  
  return(x)
}

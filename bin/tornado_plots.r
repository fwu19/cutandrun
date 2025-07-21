#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)

## functions ####
make_reference <- function(
    ref_bed,
    ref_type = c('center', 'start', 'end', 'region'), 
    ignore_strand = T,
    split_col = NULL,
    order_col = NULL
){
    require(dplyr)
    require(GenomicRanges)
    require(plyranges)
    
    ## process reference ####
    if (!file.exists(ref_bed)){ return(NULL)}
    ref <- read.delim(ref_bed, header = F)
    
    if(!is.null(split_col)){
        ref$split_by <- ref[, split_col]
    }else{
        ref$split_by <- 'no_split'
    }
    
    if(!is.null(order_col)){
        ref$order_by <- ref[, order.col]
    }
    
    ref_gr <- makeGRangesFromDataFrame(df = ref, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', strand.field = 'V6', keep.extra.columns = T, starts.in.df.are.0based = T, ignore.strand = ignore_strand)
    
    if (ref_type[1] == 'center'){
        ref_regions <- plyranges::mutate(
            plyranges::anchor_center(ref_gr), width=1)
    }else if (ref_type[1] == 'start'){
        ref_regions <- plyranges::mutate(
            plyranges::anchor_5p(ref_gr), width=1)
    }else if (ref_type[1] == 'end'){
        ref_regions <- plyranges::mutate(
            plyranges::anchor_3p(ref_gr), width=1)
    }else if (ref_type[1] == 'region'){
        ref_regions <- ref_gr
    }
 
    return(ref_regions)   
}

compute_matrix <- function(ref_regions, bw_files, extend = 1000, w = 50, target_ratio = 0.5, k = 40){
    require(GenomicRanges)
    require(plyranges)
    require(EnrichedHeatmap)
    
    bw_list <- lapply(bw_files, read_bigwig)

    ## Compute normalized matrix
    if (width(ref_regions[1]) == 1){
        mtx_list <- lapply(
            bw_list, 
            normalizeToMatrix,
            value_column = "score", 
            target = ref_regions, extend = extend,
            mean_mode = "w0", w = w, keep = c(0,1)
        )
        
    }else{
        mtx_list <- lapply(
            bw_list, 
            normalizeToMatrix,
            value_column = "score", 
            target = ref_regions, extend = extend,
            target_ratio = target_ratio, k = k,
            mean_mode = "w0", w = w, keep = c(0,1)
        )
        
    }
    

    ## return
    return(mtx_list)
}

make_one_heatmap <- function(mtx, name, row_title=NULL, nsplit=1, axis = T, show_legend = T, legend.title = 'CPM', col_fun = circlize::colorRamp2(c(0, ceiling(max(colMeans(mtx)))), c("white", "red")), axis_name = '', ylim=range(colMeans(mtx))){
    require(EnrichedHeatmap)
    require(plyranges)
    
    EnrichedHeatmap(
        mtx, 
        name = name,
        col = col_fun,
        top_annotation = HeatmapAnnotation(
            lines = anno_enriched(
                ylim = ylim,
                axis = axis, axis_param = list(side = 'right', facing = 'inside'),
                show_legend = F,
                show_annotation_name = F,
                gp = gpar(fontsize = 6, col = 1:nsplit)
            )
        ),
        axis_name = axis_name,
        axis_name_rot = 90, axis_name_gp = gpar(fontsize = 6),
        column_title = name, column_title_gp = gpar(fontsize = 6),
        row_title = row_title, row_title_gp = gpar(fontsize = 6),
        show_row_names = F, 
        show_heatmap_legend = show_legend,
        heatmap_legend_param = list(
            title = legend.title, 
            title_gp = gpar(fontsize = 6, fontface = "bold"),
            labels_gp = gpar(fontsize = 6)),
        cluster_row_slices = FALSE
    )
    
}

heatmaps_nosplit_rows <- function(
    ref_regions,
    mtx_list,  
    order_by = NULL,
    out_base = 'tornado_plots', 
    cols = c("steelblue", "yellow", "red"), 
    axis_name = NULL, 
    fig_width = NULL,
    fig_height = 6,
    hlim = NULL
){
    require(dplyr)
    require(plyranges)
    require(EnrichedHeatmap)
    
    ## compute ranges for heatmaps and line plots ####
    ylim <- range(unlist(lapply(mtx_list, colMeans)))
    ylim <- c(ifelse(ylim[1]>0, ylim[1]*0.9, ylim[1]*1.1), ifelse(ylim[2]>0, ylim[2]*1.1, ylim[2]*0.9))
    
    if (is.null(hlim)){
        hlim <- ylim[2]
    }
    
    ## assembly heatmaps ####
    ngroups <- length(mtx_list)
    ht_list = make_one_heatmap(mtx = mtx_list[[1]], nsplit = 1, name = names(mtx_list)[1], legend.title = names(mtx_list)[1], axis_name = axis_name, ylim = ylim, col_fun = circlize::colorRamp2(seq(0, hlim, length.out = length(cols)), cols))
    for (i in names(mtx_list)[2:ngroups]){
        ht_list <- ht_list+
            make_one_heatmap(mtx = mtx_list[[i]], nsplit = 1, name = i, legend.title = i, axis_name = axis_name, ylim = ylim, col_fun = circlize::colorRamp2(seq(0, hlim, length.out = length(cols)), cols))
    }
    
    ## order rows ####
    if ('order_by' %in% colnames(ref_regions)){
        row_order <- order(ref_regions$order_by, decreasing = T)
    }else if (is.null(order_by)){
        row_order <- order(rowMeans(as.matrix(sapply(mtx_list, enriched_score))), decreasing = T)
    }else{
        row_order <- order(rowMeans(as.matrix(sapply(mtx_list[order_by], enriched_score))), decreasing = T)
    }
    
    ## add a side bar if a custom order is used ####
    fwidth <- ngroups + 1
    if ('order_by' %in% colnames(ref_regions)){
        ht_list <-  Heatmap(
            ref_regions$order_by, name = ".", row_title = NULL, 
            show_row_names = FALSE, width = unit(3, "mm"), 
            heatmap_legend_param = list(
                title = order_col, 
                title_gp = gpar(fontsize = 6, fontface = "bold"),
                labels_gp = gpar(fontsize = 6)),
            cluster_row_slices = FALSE
        )+
            ht_list
        fwidth <- fwidth + 0.3  
    }
    
    ## save to a file #### 
    out_dir <- dirname(out_base)
    if(!dir.exists(out_dir)){dir.create(out_dir, recursive = T)}
    
    pdf(paste(out_base, 'pdf', sep = '.'), width = ifelse(is.null(fig_width), fwidth, fig_width), height = fig_height)
    draw(ht_list, row_order=row_order, heatmap_legend_side = "right", ht_gap = unit(2, "mm"))
    dev.off()
    
}

heatmaps_split_rows <- function(
    ref_regions,
    mtx_list,  
    order_by = NULL,
    out_base = 'tornado_plots', 
    cols = c("steelblue", "yellow", "red"), 
    axis_name = NULL, 
    fig_width = NULL,
    fig_height = 6,
    hlim = NULL
){
    require(dplyr)
    require(plyranges)
    require(EnrichedHeatmap)
    
    ## split rows and compute ranges ####
    nsplit <- length(unique(ref_regions$split_by))
    ylim <- range(
        sapply(
            mtx_list, 
            function(mtx){
                range(sapply(split(1:nrow(mtx), ref_regions$split_by), function(i){colMeans(mtx[i,,drop=F])}))
            })
    )
    ylim <- c(ifelse(ylim[1]>0, ylim[1]*0.9, ylim[1]*1.1), ifelse(ylim[2]>0, ylim[2]*1.1, ylim[2]*0.9))
    if (is.null(hlim)){
        hlim <- ylim[2]
    }
    
    ## assembly heatmaps ####
    ngroups <- length(mtx_list)
    ht_list <- make_one_heatmap(mtx = mtx_list[[1]], nsplit = nsplit, name = names(mtx_list)[1], legend.title = names(mtx_list)[1], axis_name = axis_name, ylim = ylim, col_fun = circlize::colorRamp2(seq(0, hlim, length.out = length(cols)), cols))
    for (i in names(mtx_list)[2:ngroups]){
        ht_list <- ht_list+
            make_one_heatmap(mtx = mtx_list[[i]], nsplit = nsplit, name = i, legend.title = i, axis_name = axis_name, ylim = ylim, col_fun = circlize::colorRamp2(seq(0, hlim, length.out = length(cols)), cols))
    }
    
    ## order rows ####
    if ('order_by' %in% colnames(ref_regions)){
        row_order <- order(ref_regions$order_by, decreasing = T)
    }else if (is.null(order_by)){
        row_order <- order(rowMeans(as.matrix(sapply(mtx_list, enriched_score))), decreasing = T)
    }else{
        row_order <- order(rowMeans(as.matrix(sapply(mtx_list[order_by], enriched_score))), decreasing = T)
    }
    
    ## add a side bar if a custom order is used ####
    fwidth <- ngroups + 1
    if ('order_by' %in% colnames(ref_regions)){
        ht_list <-  Heatmap(
            ref_regions$order_by, name = ".", row_title = NULL, 
            show_row_names = FALSE, width = unit(3, "mm"), 
            heatmap_legend_param = list(
                title = order_col, 
                title_gp = gpar(fontsize = 6, fontface = "bold"),
                labels_gp = gpar(fontsize = 6)),
            cluster_row_slices = FALSE
        )+
            ht_list
        fwidth <- fwidth + 0.3  
    }
    
    ## save to a file #### 
    ht_list <- ht_list +
        Heatmap(
            ref_regions$split_by, col = 1:nsplit, name = ".", row_title = NULL, 
            show_row_names = FALSE, width = unit(3, "mm"), 
            heatmap_legend_param = list(
                title = NULL, 
                title_gp = gpar(fontsize = 0),
                labels_gp = gpar(fontsize = 6)),
            cluster_row_slices = FALSE
        )
    fwidth <- fwidth + 0.3

    out_dir <- dirname(out_base)
    if(!dir.exists(out_dir)){dir.create(out_dir, recursive = T)}
    pdf(paste(out_base, 'pdf', sep = '.'), width = ifelse(is.null(fig_width), fwidth, fig_width), height = fig_height)
    draw(ht_list, row_order=row_order, row_split = ref_regions$split_by, heatmap_legend_side = "right", ht_gap = unit(2, "mm"))
    dev.off()
    
    
}

## make plots ####
args <- as.vector(commandArgs(T)) # 
for (arg in args){
    v <- unlist(strsplit(arg, split = '='))
    assign(v[1], v[2])
}                                                                        
if(!exists('ref_type')){ref_type <- 'center'}
if(!exists('ignore_strand')){ignore_strand <- T}
ignore_strand <- as.logical(ignore_strand)
if(!exists('axis_name')){ 
    axis_name <- c('-1Kb', 'Center', '+1Kb')
}else{
    axis_name <- unlist(strsplit(axis_name, split = ','))
}
ref_beds <- file.path('bed', list.files('bed/'))
bw_files <- file.path('bigwig', list.files('bigwig/'))

for (ref_bed in ref_beds){
    if (file.size(ref_bed) < 1){ next }
    
    ref_regions <- make_reference(
        ref_bed, ref_type, ignore_strand, 
        order_col = NULL, split_col = NULL)
    
    mtx_list <- compute_matrix(ref_regions, bw_files)
    names(mtx_list) <- tools::file_path_sans_ext(basename(bw_files))

    out_base <- tools::file_path_sans_ext(basename(ref_bed))
    if(!grepl(target, out_base)){out_base <- paste(out_base, target, sep = '_')}
    
    heatmaps_nosplit_rows(
        ref_regions = ref_regions,
        mtx_list = mtx_list,
        out_base = out_base,
        hlim = 4,
        axis_name = axis_name
    )
}


## make heatmaps around TSS
gene_beds <- file.path('gene/', list.files('gene/'))
for (ref_bed in gene_beds){
    if (file.size(ref_bed) < 1){ next }
    
    ref_type <- 'start'
    ignore_strand <- F
    out_base <- 'gene_TSS'
    if(!grepl(target, out_base)){out_base <- paste(out_base, target, sep = '_')}
    axis_name <- c('-1Kb', 'TSS', '+1Kb')
    
    ref_regions <- make_reference(
        ref_bed, ref_type = 'start', ignore_strand = F, 
        order_col = NULL, split_col = NULL)
    
    mtx_list <- compute_matrix(ref_regions, bw_files)
    names(mtx_list) <- tools::file_path_sans_ext(basename(bw_files))
    
    heatmaps_nosplit_rows(
        ref_regions = ref_regions,
        mtx_list = mtx_list,
        out_base = out_base,
        hlim = 4,
        axis_name = axis_name
    )
}






# heatmaps_split_rows(
#     ref_regions = ref_regions,
#     mtx_list = mtx_list,
#     out_prefix = out_prefix, 
#     hlim = 4,
#     axis_name = axis_name
# )



#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(ggplot2)
library(patchwork)

## functions ####
count2dgelist <- function(counts.tsv=NULL, return.counts = F, pattern2remove="^X|.bam$", counts=NULL, out.dir=NULL, feature.cols=1:7, samples = NULL){
    options(stringsAsFactors = F)
    require(edgeR)
    
    if(is.null(counts)){
        counts <- read.delim(counts.tsv)
    }
    if(!is.null(pattern2remove)){
        colnames(counts) <- gsub(pattern2remove, '', colnames(counts))
    }
    
    y0 <- DGEList(counts=counts[-(feature.cols)], genes=counts[feature.cols], remove.zeros = T, samples = samples)
    y0 <- calcNormFactors(y0)
    
    if(return.counts){
        if(is.null(out.dir)){out.dir <- dirname(counts.tsv)}
        if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
        write.table(cbind(y0$genes, y0$counts), paste(out.dir, 'rawCounts.txt', sep = '/'), quote = F, row.names = F)
    }
    return(y0)
}

run_da <- function(
    y0, out.prefix,
    control.group, test.group, group=NULL,
    fdr=0.01, fc=1, fdr2=NULL, fc2=NULL,
    report.cpm=F, report.rpkm=T,
    TMM=T, method = 'QL',
    rename.feature = NULL, feature.length = 'length',
    design.object = ~0+group,
    target = NULL
){
    require(edgeR)

    ## retrieve and process data ####
    if(!is.null(group)){y0$samples$group <- group}
    if (sum(y0$samples$group %in% control.group) < 2 & sum(y0$samples$group %in% test.group) < 2){
        return(NULL)
    }else if (sum(y0$samples$group %in% control.group) == 0 | sum(y0$samples$group %in% test.group) == 0){
        return(NULL)
    }
    j <- y0$samples$group %in% c(control.group, test.group)

    y <- y0[,j]
    y$samples$group <- ifelse(y$samples$group %in% control.group, 'control', 'test')

    if(TMM){
        keep <- filterByExpr(y, group = y$samples$group, min.count=10, min.total.count = 15)

        y <- y[keep,,keep.lib.sizes=F]
        y<-calcNormFactors(y)
    }


    ## Create design matrix ####
    design <- model.matrix(design.object,data = y$samples)
    colnames(design) <- gsub('^group','',colnames(design))

    ## Make contrasts ####
    contrasts <- makeContrasts(
        cmp = test - control,
        levels = design
    )

    ## Run DP test ####
    out.dir <- dirname(out.prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
    
    y <- estimateDisp(y, design, robust = T)

    pdf(paste(out.dir,'bcv.pdf',sep = '/'),width = 4, height = 4); plotBCV(y); dev.off()

    if (method == 'QL'){
        fit<-glmQLFit(y, design=design, dispersion = y$trended.dispersion, robust = T)

        # pdf(paste(out.dir,'qldisp.pdf',sep = '/'),width = 4, height = 4); plotQLDisp(fit); dev.off()

        test <- glmQLFTest(fit, contrast = contrasts)
    }else{
        fit <- glmFit(y, design=design, dispersion = y$trended.dispersion, robust = T)
        test <- glmLRT(fit, contrast = contrasts)
    }

    ## prepare for plots
    df <- test$table
    df$FDR <- p.adjust(df$PValue, method = 'BH')
    df$is.sig <- (df$FDR < fdr) * sign(df$logFC) * (abs(df$logFC) > log2(fc))

    ## convert y$group back ####
    y$samples$group <- ifelse(y$samples$group %in% 'control', control.group, test.group)

    ## write out results ####
    
    if(!is.null(rename.feature)){colnames(test$genes)[1] <- rename.feature}
    df <- cbind(test$genes,df[c('logFC','logCPM','PValue','FDR','is.sig')])
    if(!is.null(fdr2) & !is.null(fc2)){
        df$is.sig2 <- (df$FDR < fdr2) * sign(df$logFC) * (abs(df$logFC) > log2(fc2))
    }
    if(report.cpm){
        cpm <- cpm(y, normalized.lib.sizes = T, log = F)
        colnames(cpm) <- paste('CPM.TMMnormalized', colnames(cpm), sep = '.')
        df <- cbind(df, cpm)
    }
    if(report.rpkm){
        rpkm <- rpkm(y, gene.length = feature.length, normalized.lib.sizes = T, log = F)
        colnames(rpkm) <- paste('FPKM.TMMnormalized', colnames(rpkm), sep = '.')
        df <- cbind(df, rpkm)
    }
    write.table(df,paste(out.prefix,'txt',sep = '.'), sep = '\t',quote = F,row.names = F)

    ## return results ####
    df_sum <- data.frame(
        output.folder = gsub('.*differential_peaks/', '', out.dir),
        control.group = control.group,
        test.group = test.group,
        control.samples = sum(y$samples$group %in% control.group),
        test.samples = sum(y$samples$group %in% test.group),
        features.tested = nrow(y),
        features.up = sum(df$is.sig %in% 1),
        features.down = sum(df$is.sig %in% -1),
        FDR.cutoff = fdr,
        FC.cutoff = fc
    )
    if(!is.null(target)){
        df_sum$target <- target
    }

    if('is.sig2' %in% colnames(df)){
        df_sum <- cbind(
            df_sum,
            data.frame(
                features.up2 = sum(df$is.sig2 %in% 1),
                features.down2 = sum(df$is.sig2 %in% -1),
                FDR.cutoff2 = fdr2,
                FC.cutoff2 = fc2
            )
        )
    }

    return(list(summary = df_sum, y = y, df = df, design=design, test=test))

}


## PCA
plot_pca <- function(y, out.prefix, var.genes = NULL, color = NULL, plot.title = '', sample.label = T, feature.length = 'gene_length'){
    options(stringsAsFactors = F)
    require(ggrepel)
    require(edgeR)

    ## Create outdir if needed
    out.dir <- dirname(out.prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir,recursive = T)}


    log2rpkm <- rpkm(y, gene.length = feature.length, normalized.lib.sizes = T, log = T)

    if (!is.null(var.genes)){
        keep <- rank(-apply(log2rpkm, 1, var)) < var.genes
        log2rpkm <- log2rpkm[keep,]
    }

    ## Run PCA
    pca <- prcomp(t(log2rpkm), center = T, scale = T)

    ## Scree plot
    pca.variance.prop <- (pca$sdev^2)/sum(pca$sdev^2)*100

    pdf(paste(out.dir,'pca.scree.plot.pdf',sep = '/'))
    barplot(
        pca.variance.prop[1:50],
        cex.names = 1,
        xlab = 'Principal component (PC), 1-50',
        ylab = 'Proportion of variance (%)',
        main = 'Scree plot',
        ylim = c(0,80)
    )

    points(
        cumsum(pca.variance.prop)[1:50], col = 'red', type = 'l'
    )
    dev.off()

    ## PC1 vs PC2
    df <- cbind(pca$x[,1:2],data.frame(label=rownames(pca$x)))
    if(is.null(color)){
        df$color <- y$samples$group
    }else{
        df$color <- color
    }

    if(length(unique(df$color))>1){
        p <- ggplot(df,aes(x=PC1,y=PC2,label=label,color=color))
    }else{
        p <- ggplot(df,aes(x=PC1,y=PC2,label=label))
    }
    p <- p +
        geom_point(shape = 1)

    if(sample.label){
        p <- p +
            geom_text_repel(size = 2.4, color = 'black', position = 'jitter',max.overlaps = 80)
    }

    p <- p +
        labs(
            title = plot.title,
            color = '',
            x = paste0('PC1 (',round(pca.variance.prop[1],1),'%)'),
            y = paste0('PC2 (',round(pca.variance.prop[2],1),'%)')
        )+
        theme_bw()+
        theme(
            legend.position = 'top'
        )
    ggsave(paste(out.prefix,'PCA.pdf',sep='.'),width = 5,height = 5)

    ## return data
    return(p)

}

## Plot MD
plot_MD <- function(df, out.prefix, plot.title = ""){
    require(ggplot2)

    ## Create outdir if needed
    out.dir <- dirname(out.prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir,recursive = T)}

    p <- ggplot(df, aes(x = logCPM, y = logFC, color=factor(is.sig)))+
        geom_hline(yintercept = 0)+
        geom_point(size = 0.2)+
        scale_color_manual(
            '',
            values = c("-1" = "blue", "0" = "gray", "1" = "red"),
            breaks = c(-1,0,1), labels = c('Down', 'No Sig.', 'Up')
        )+
        labs(
            x = "Average log CPM",
            y = "log-fold-of-change",
            title = plot.title
        )+
        theme_bw()+
        theme(
            text = element_text(size = 8),
            legend.position = 'top'
        )

    ggsave(paste(out.prefix,'MD.pdf',sep = '.'),width = 4,height = 5)
    return(p)
}

## Plot volcano
plot_volcano <- function(df, out.prefix, plot.title = ""){
    require(ggplot2)

    ## Create outdir if needed
    out.dir <- dirname(out.prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir,recursive = T)}

    p <- ggplot(df,aes(x=logFC,y=-log10(FDR),color=factor(is.sig)))+
        geom_point(size = 0.2)+
        scale_color_manual(
            values = c('-1'='blue','0'='gray','1'='red'),
            breaks = c('-1','0','1'),
            labels = c('Down','No Sig.','Up'),
            drop = T
        )+
        labs(
            x='log2(fold change)',
            y='-log10(FDR)',
            color='',
            title = plot.title
        )+
        theme_bw()+
        theme(
            text = element_text(size = 8),
            legend.position = 'top'
        )

    ggsave(paste(out.prefix,'Volcano.pdf',sep = '.'),width = 4,height = 5)
    return(p)
}

## recompute is.sig2
recal_sig <- function(txt, col.sig, fdr, fc){
    de <- read.delim(txt)
    de[,col.sig] <- sign(de$logFC) * (abs(de$logFC) > log2(fc)) * (de$FDR < fdr)
    write.table(de, txt, sep = '\t', quote = F, row.names = F)
    return(
        data.frame(
            source.file = basename(txt),
            modify.col = col.sig,
            features.test = nrow(de),
            features.up = sum(de[,col.sig] %in% 1),
            features.down = sum(de[,col.sig] %in% -1),
            FDR.cutoff = fdr,
            FC.cutoff = fc
        )
    )
}

## wrapper
run_one_comparison <- function(y0, control.group, test.group, out.dir, prefix, plot.title, fdr, fc, fdr2, fc2, tgt){
    k <- sapply(
        strsplit(y0$genes$sample.groups, split = ','),
        function(v){sum(c(control.group, test.group) %in% v) > 0 }) > 0 # filter peaks present in either control or test group
    if(sum(k) == 0){ return(NULL)}

    out.prefix <- paste(out.dir, prefix, prefix, sep = '/')
    lst <- run_da(
        y0[k,],
        out.prefix,
        control.group = gsub('-', '_', control.group),
        test.group = gsub('-', '_', test.group),
        group = gsub('-', '_', y0$samples$sample_group),
        feature.length = 'length',
        fdr = fdr, fc = fc, fdr2 = fdr2, fc2 = fc2,
        target = tgt
    )

    if(is.null(lst)){return(NULL)}

    y <- lst$y
    df <- lst$df
    lst$plots <- list(
        PCA = plot_pca(
            y,
            out.prefix,
            color = y$samples$sample_group,
            sample.label = T,
            plot.title = "",
            var.genes = 500,
            feature.length = "length"),
        MD = plot_MD(
            df, out.prefix = out.prefix,
            plot.title = plot.title
        ),
        volcano = plot_volcano(
            df, out.prefix = out.prefix,
            plot.title = plot.title
        )
    )


    return(lst)


}

wrapper_one_conp <- function(ss, cmp, tgt, cts.file, out.base, fdr, fc, fdr2, fc2){
    
    cts <- read.delim(cts.file)
    if (nrow(cts) < 10){ return (NULL)} # do not test if less than 10 peaks
    
    ## create DGElist ####
    ssi <- ss %>%
        filter(id %in% colnames(cts)) %>%
        dplyr::select(id, target, sample_group, sample_replicate) %>%
        arrange(factor(id, levels = colnames(cts)[8:ncol(cts)]))
    
    y0 <- count2dgelist(
        counts = cts,
        feature.cols = 1:7,
        samples = ssi
    )
    saveRDS(y0, paste0(out.base, '.y0.rds'))
    
    
    ## run DGE ####
    dp <- mapply(
        run_one_comparison,
        MoreArgs = list(y0 = y0, out.dir = out.base, fdr = fdr, fc = fc, fdr2 = fdr2, fc2 = fc2, tgt = tgt),
        cmp$control.group,
        cmp$test.group,
        cmp$out.prefix,
        cmp$plot.title,
        SIMPLIFY = F)
    names(dp) <- basename(cmp$out.prefix)

    if(length(dp) > 0){
        dp[sapply(dp, is.null)] <- NULL
    }

    if (length(dp) > 0){
        return(dp)
    }else{
        return(NULL)
    }

}


## read arguments ####
args <- as.vector(commandArgs(T)) # ss, cmp

lst <- strsplit(args, split = '=')
for (x in lst){
    assign(x[1],x[2])
} # read arguments: ss, comparison, cts, fdr, fc, fdr2, fc2
rm(lst)

ss <- read.csv(ss_csv)

if (grepl('dummy_file', cmp_file)){
    error_message <- paste(cmp_file, "is a dummy file! Provide --comparison path/to/comparison_file (a comparison table in csv, txt, tsv or rds format)!")
    write.table(
        error_message, paste0(tgt, '.README.txt'), sep = '\n', quote=F, row.names = F, col.names = F
    )
    quit()
}else if (file.size(cmp_file) == 0){
   error_message <- paste( cmp_file, "is empty!")
   write.table(
       error_message, paste0(tgt, '.README.txt'), sep = '\n', quote=F, row.names = F, col.names = F
   )
   quit()
}else if (grepl('.csv$', cmp_file)){
    cmp <- read.csv(cmp_file)
}else if (grepl('.rds$', cmp_file)){
    cmp <- readRDS(cmp_file)
}else if (grepl('.txt$|.tsv$', cmp_file)){
    cmp <- read.delim(cmp_file)
}else{
    stop(paste(cmp_file, "should be .csv, .txt, .tsv or .rds!"))
}
if(!'out.prefix' %in% colnames(cmp)){
    cmp$out.prefix <- paste(cmp$test.group, cmp$control.group, sep = '_vs_')
}
if(!'plot.title' %in% colnames(cmp)){
    cmp$plot.title <- paste(cmp$test.group, cmp$control.group, sep = ' vs ')
}

if (exists('fdr')){ fdr <- as.numeric(fdr) }else{ fdr <- 0.05 }
if (exists('fc')){ fc <- as.numeric(fc) }else{ fc <- 1.5 }
if (exists('fdr2')){ fdr2 <- as.numeric(fdr2) }else{ fdr2 <- 0.01 }
if (exists('fc2')){ fc2 <- as.numeric(fc2) }else{ fc2 <- 2 }

cts.files <- list.files('counts/', full.names = T)

## detect differential peaks ####
dp.list <- list()
for (cts.file in cts.files){
    out.base <- gsub('\\.raw_counts.txt$', '', basename(cts.file))
    dp.list[[out.base]] <- wrapper_one_conp(ss, cmp, tgt, cts.file, out.base, fdr, fc, fdr2, fc2)

}
if (length(dp.list) > 0){
    dp.list[sapply(dp.list,is.null)] <- NULL
}

if (length(dp.list) == 0){
    error_message <- c(
        "No comparison was done for this target.",
        "Check the following: ",
        paste("Test and control groups in", args[2], "should match sample_groups in", args[1], "."),
        paste("At least one row in", args[2], "should be a valid comparison, i.e. both test and control groups are present in", args[1], ".")
        )
    write.table(
        error_message, paste0(tgt, '.README.txt'), sep = '\n', quote=F, row.names = F, col.names = F
    )
    system('rm -r *_peaks/')
}else{
    saveRDS(dp.list, paste(tgt, 'dp.rds', sep = '.'))
}


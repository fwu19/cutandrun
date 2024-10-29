#!/usr/bin/env Rscript

options(stringsAsFactors = F)
options(scipen = 99)
library(dplyr)
library(ggplot2)
library(patchwork)

## functions ####
generate_count_matrix <- function(conp.bed, infiles){

    conp <- read.delim(conp.bed, header = F, col.names = c('chrom', 'start', 'end', 'conp.id', 'sample.groups', 'npeak'))
    
    cts <- cbind(
        conp[1:5],
        read.delim(infiles[1], header = T, comment.char = '#')[6],
        lapply(
        infiles, function(fname){
            df <- read.delim(fname, header = T, comment.char = '#')
            stopifnot(identical(conp$conp.id, df$Geneid))
            df[7]
        }
    ))
    colnames(cts)[7:ncol(cts)] <- gsub('.fragmentCounts.txt', '', basename(infiles))
    cts
}

count2dgelist <- function(counts.tsv=NULL, return.counts = T, pattern2remove="^X|.bam$", counts=NULL, out.dir=NULL, feature.cols=1:7, samples = NULL){
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
    
    if(is.null(out.dir)){out.dir <- dirname(counts.tsv)}
    if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
    saveRDS(y0, paste(out.dir, 'y0.rds', sep = '/'))
    
    if(return.counts){
        write.table(cbind(y0$genes, y0$counts), paste(out.dir, 'rawCounts.txt', sep = '/'), quote = F, row.names = F)
    }
    return(y0)
}

## for all
pca_log2rpkm <- function(y, out.dir, prefix, var.genes = NULL, color = NULL, plot.title = '', sample.label = T, feature.length = 'gene_length'){
    options(stringsAsFactors = F)
    require(ggrepel)
    
    log2rpkm <- rpkm(y, gene.length = feature.length, normalized.lib.sizes = T, log = T)
    
    if (!is.null(var.genes)){
        keep <- rank(-apply(log2rpkm, 1, var)) < var.genes
        log2rpkm <- log2rpkm[keep,]
    }
    
    ## Run PCA
    pca <- prcomp(t(log2rpkm), center = T, scale = T)
    
    ## Create outdir if needed
    if(!dir.exists(out.dir)){dir.create(out.dir,recursive = T)}
    
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
    
    p +
        labs(
            title = plot.title,
            color = '', 
            x = paste0('PC1 (',round(pca.variance.prop[1],1),'%)'),
            y = paste0('PC2 (',round(pca.variance.prop[2],1),'%)')
        )+
        theme_bw()
    ggsave(paste(out.dir,paste(prefix,'PCA.pdf',sep='.'),sep = '/'),width = 6,height = 5)
    
    ## return data
    return(pca)
    
}

## for all
run_da <- function(
    y0, out.prefix, 
    control.group, test.group, group=NULL, 
    fdr=0.01, lfc=1, fdr2=NULL, lfc2=NULL, 
    report.cpm=F, report.rpkm=T, 
    TMM=T, method = 'QL', 
    rename.feature = NULL, feature.length = 'length', 
    design.object = ~0+group,
    target = NULL
){
    require(edgeR)
    
    ## create output directory ####
    out.dir <- dirname(out.prefix)
    if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
    
    ## retrieve and process data ####
    if(!is.null(group)){y0$samples$group <- group}
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
    y <- estimateDisp(y, design, robust = T)
    
    pdf(paste(out.dir,'bcv.pdf',sep = '/'),width = 4, height = 4); plotBCV(y); dev.off()
    
    if (method == 'QL'){
        fit<-glmQLFit(y, design=design, dispersion = y$trended.dispersion, robust = T)
        
        pdf(paste(out.dir,'qldisp.pdf',sep = '/'),width = 4, height = 4); plotQLDisp(fit); dev.off()
        
        test <- glmQLFTest(fit, contrast = contrasts)
    }else{
        fit <- glmFit(y, design=design, dispersion = y$trended.dispersion, robust = T)
        test <- glmLRT(fit, contrast = contrasts)
    }
    
    ## prepare for plots
    is.sig <- decideTests(object = test,adjust.method = 'BH',p.value = fdr,lfc = lfc)
    df <- test$table
    df$FDR <- p.adjust(df$PValue, method = 'BH')
    df$is.sig <- as.vector(is.sig)
    
    ## convert y$group back ####
    y$samples$group <- ifelse(y$samples$group %in% 'control', control.group, test.group)
    
    ## write out results ####
    saveRDS(list(y=y, design=design, fit=fit, test=test), paste(out.dir,'da.rds',sep = '/'))
    if(!is.null(rename.feature)){colnames(test$genes)[1] <- rename.feature}
    df <- cbind(test$genes,df[c('logFC','logCPM','PValue','FDR','is.sig')])
    if(!is.null(fdr2) & !is.null(lfc2)){
        df$is.sig2 <- (df$FDR < fdr2) * sign(df$logFC) * (abs(df$logFC) > lfc2)
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
        features.up = sum(is.sig %in% 1),
        features.down = sum(is.sig %in% -1),
        FDR.cutoff = fdr,
        FC.cutoff = round(2^lfc,1)
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
                FC.cutoff2 = round(2^lfc2,1)
            )
        )
    }
    
    return(list(summary = df_sum, y = y, df = df))
    
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
recal_sig <- function(txt, col.sig, fdr, lfc){
    de <- read.delim(txt)
    de[,col.sig] <- sign(de$logFC) * (abs(de$logFC) > lfc) * (de$FDR < fdr)
    write.table(de, txt, sep = '\t', quote = F, row.names = F)
    return(
        data.frame(
            source.file = basename(txt),
            modify.col = col.sig,
            features.test = nrow(de),
            features.up = sum(de[,col.sig] %in% 1),
            features.down = sum(de[,col.sig] %in% -1),
            FDR.cutoff = fdr,
            FC.cutoff = round(2^lfc,1)
        )
    )
}

## wrapper
wrapper_one_conp <- function(ss, cmp, cnp, in.dir, out.dir, fdr = 0.05, lfc = log2(1.5), fdr2 = 0.01, lfc2 = 1){
    out.dir <- paste(out.dir, gsub('\\.bed$', '', basename(cnp)), sep = '/')
    if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
    
    ## generate count matrix ####
    ssi <- ss %>%
        filter(conp == cnp) 
    tgt <- ssi$target[1]
    cts <- generate_count_matrix(paste(in.dir, '07_consensus_peaks', cnp, sep = '/'), paste(in.dir, '07_consensus_peaks/featureCounts', sub('.bed$', '', cnp), paste(ssi$sample_id, 'fragmentCounts.txt', sep = '.'), sep = '/'))
    
    ## create DGElist ####
    y0 <- count2dgelist(
        counts = cts, 
        out.dir = out.dir, 
        feature.cols = 1:6, 
        samples = ssi %>% 
            dplyr::select(sample_group, sample_replicate, sample_id) %>% 
            arrange(factor(sample_id, levels = colnames(cts)[7:ncol(cts)]))
    )
    
    ## subset comparison table ####
    icmp <- cmp %>% filter(target == tgt)
    
    ## run DGE ####
    run_one_comp <- function(control, test){
        k <- grepl(paste(c(control, test), collapse = '|'), y0$genes$sample.groups) # filter peaks present in samples for comparison
        if(sum(k) > 0){
            out.prefix <- paste(out.dir, paste(test,control,sep = '_vs_'), paste(test,control,sep = '_vs_'), sep = '/')
            lst <- run_da(
                y0[k,], 
                out.prefix, 
                control.group = gsub('-', '_', control), 
                test.group = gsub('-', '_', test), 
                group = gsub('-', '_', y0$samples$sample_group), 
                feature.length = 'Length',
                fdr = fdr, lfc = lfc, fdr2 = fdr2, lfc2 = lfc2,
                target = tgt
            )
            
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
                    feature.length = "Length"),
                MD = plot_MD(
                    df, out.prefix = out.prefix, 
                    plot.title = paste0(tgt, ': ', test, " vs ", control, ' (control)')
                ),
                volcano = plot_volcano(
                    df, out.prefix = out.prefix, 
                    plot.title = paste0(tgt, ': ', test, " vs ", control, ' (control)')
                )
            )
            
            
            return(lst)
        }
        
    }
    
    dp <- mapply(run_one_comp, icmp$control, icmp$test, SIMPLIFY = F)
    names(dp) <- paste(icmp$test, icmp$control, sep = '_vs_')
    
    return(dp)
}


## read arguments ####
args <- as.vector(commandArgs(T)) # ss, ssdp, cmp
ss <- read.csv(args[1]) %>% 
    filter(!grepl('IgG', target)) %>%
    mutate(
        sample_id = paste(group, replicate, sep = '_R')
    ) %>% 
    right_join(
        read.csv(args[2]), by = 'sample_id'
    ) 
cmp <- read.csv(args[3])
in.dir <- normalizePath(args[4])
out.dir <- 'differential_peaks'
if (!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

## detect differential peaks ####
dp.list <- lapply(unique(ss$conp), wrapper_one_conp, ss=ss, cmp=cmp, in.dir = in.dir, out.dir = out.dir)
names(dp.list) <- unique(ss$conp)
saveRDS(dp.list, paste(out.dir, 'dp.rds', sep = '/'))

dp_sum <- bind_rows(lapply(dp.list, function(dp){
    bind_rows(lapply(dp, function(x){x$summary}))
})) %>% 
    dplyr::relocate(target, output.folder, control.group, test.group, features.tested:FC.cutoff, control.samples:test.samples)
dp_sum %>% 
    write.table(paste(out.dir, 'DP_summary.txt', sep = '/'), sep = '\t', quote = F, row.names = F)

suppressWarnings(suppressMessages(library(balance)))
suppressWarnings(suppressMessages(library(zCompositions)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(ggpubr)))
pound2inch = 1/72
Linesp=1.5

sbp.subset <- function (sbp, ternary = TRUE, ratios = TRUE){
    if (!ternary & !ratios) {
        message("Alert: Skipping balance subset.")
        return(sbp)
    }
    if (ternary & !ratios) {
        message("Alert: Using 'ternary' enables 'ratios' too.")
        ratios <- TRUE
    }
    b.size <- apply(sbp, 2, function(x) sum(abs(x)))
    keep <- rep(FALSE, ncol(sbp))

                    keep <- keep | (b.size == 3)
    if (ratios) 
        keep <- keep | (b.size == 2)
    sbp <- sbp[, keep, drop = FALSE]
    joinNames <- function(bp, x) make.names( 
        paste0(
        c(paste0(x[bp==1],collapse = "_"), #分子 
          paste0(x[bp==-1], collapse = "_")),  #分母
        collapse="__") #分子__分母1_分母2 or 分子1_分子2__分母
    ) 
    colnames(sbp) <- apply(sbp, 2, function(x) joinNames(x[x !=0], rownames(sbp)[x !=0]))
    sbp
}  
                           
balance_plot_yang <- function (x, y, d.group, n.group, boxplot.split = TRUE, weigh.var = FALSE, 
    size.text = 20, size.pt = 4) 
{
    Amber=c('#ff5b47', '#8cfd94','#ffd700', '#ffa001','#8aceea', '#d4a1d1','#fd6ab3','#9775de',
        '#b9c206','#f5b5c9','#a3fe41','#ddb689','#437bb4','#dd6b28','#47d8c3','#c56dc5')
    cols <- rep(Amber,2)
    if (ncol(x) != nrow(y)) {
        stop("Please check that ncol(x) = nrow(y) = D.")
    }
    if (any(x == 0)) {
        stop("Please remove zeros before analysis.")
    }
    if (is.null(colnames(x))) {
        colnames(x) <- paste0("Component", 1:ncol(x))
    }
    if (is.null(colnames(y))) {
        colnames(y) <- paste0("Balance", 1:ncol(y))
    }
    if (is.null(rownames(x))) {
        rownames(x) <- as.character(1:nrow(x))
    }
    x <- as.data.frame(x)
    y <- as.data.frame(y)
    rownames(y) <- colnames(x)
    b.weight <- apply(y, 2, function(i) sum(abs(i)))
    b.order <- order(b.weight, decreasing = TRUE)
    d.weight <- apply(y[, b.order], 1, function(i) sum(i * 1/1:length(i)))
    d.order <- order(d.weight, decreasing = TRUE)
    pt <- y
    pt <- wide2long(pt)
    pt <- pt[pt$value != 0, ]
    colnames(pt) <- c("Part", "BalanceID", "Component")
    pt$Component <- factor(pt$Component, colnames(x)[d.order])
    pt$BalanceID <- factor(pt$BalanceID, rev(colnames(y)[b.order]))
    pt$Part <- factor(pt$Part, levels = c(1, -1))
    pt$Group <- paste0(pt$BalanceID, pt$Part)
    B <- balance.fromSBP(x, y)
    colnames(B) <- colnames(y)
    rownames(B) <- rownames(x)
    dt <- wide2long(B)
    colnames(dt) <- c("SampleValue", "BalanceID", "Index")
    dt$BalanceID <- factor(dt$BalanceID, rev(colnames(y)[b.order]))
    if (weigh.var) {
        vars <- apply(B, 2, stats::var)
        vars <- vars/sum(vars)
        linewidth <- data.frame(BalanceID = colnames(B), LineWidth = vars)
        dt <- merge(dt, linewidth)
    }
    if (!all(pt$BalanceID %in% dt$BalanceID)) 
        stop("Unexpected Error: try renaming balances.")
    if (!all(dt$BalanceID %in% pt$BalanceID)) 
        stop("Unexpected Error: try renaming balances.")
    if (!missing(d.group)) {
        if (length(unique(d.group)) > 32) 
            stop("Only 32 groups for 'd.group' supported.")
        d.group <- as.character(d.group)
        names(d.group) <- colnames(x)
        pt$d.group <- d.group[as.character(pt$Component)]
        d.cols <- cols[1:length(unique(pt$d.group))]
    }
    else {
        pt$d.group <- "1"
        d.cols <- "black"
    }
    if (!missing(n.group)) {
        if (length(unique(n.group)) > 8) 
            stop("Only 8 groups for 'n.group' supported.")
        n.group <- as.character(n.group)
        names(n.group) <- rownames(x)
        dt$n.group <- n.group[as.character(dt$Index)]
        n.cols <- cols[1:length(unique(dt$n.group))]
    }
    else {
        dt$n.group <- "1"
        n.cols <- "black"
    }
    if (boxplot.split) {
        dt$Group <- paste(dt$BalanceID, dt$n.group)
    }
    else {
        dt$Group <- dt$BalanceID
    }
    balance.partition <- ggplot2::ggplot(pt, ggplot2::aes_string(x = "BalanceID", 
        y = "Component", shape = "Part", group = "Group")) + 
        ggplot2::geom_line() + ggplot2::geom_point(#ggplot2::aes_string(col = "d.cols"), 
        size = size.pt) + #ggplot2::scale_colour_manual(values = d.cols) + 
        ggplot2::xlab("Balance ID") + ggplot2::ylab("Component ID") + 
        ggplot2::labs(col = "Component Group", shape = "Partition") + 
        ggplot2::coord_flip() + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
        hjust = 1, vjust = 0.5)) + ggplot2::theme(text = ggplot2::element_text(size = size.text)) + 
        ggplot2::theme(legend.position = "top")
    balance.distribution <- ggplot2::ggplot(dt, ggplot2::aes_string(x = "BalanceID", 
        y = "SampleValue", group = "BalanceID"), col = "black") + 
        ggplot2::scale_colour_manual(values = n.cols) + ggplot2::xlab("") + 
        ggplot2::ylab("Sample-wise Distribution of Balance") + 
        ggplot2::ylim(-1.1 * max(abs(dt$SampleValue)), 1.1 * 
            max(abs(dt$SampleValue))) + ggplot2::labs(col = "Sample Group") + 
        ggplot2::coord_flip() + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
        hjust = 1, vjust = 0.5)) + ggplot2::theme(text = ggplot2::element_text(size = size.text)) + 
        ggplot2::theme(legend.position = "top")
    if (weigh.var) {
        balance.distribution <- balance.distribution + ggplot2::geom_line(ggplot2::aes_string(size = "LineWidth"), 
            alpha = 0.5) + ggplot2::guides(size = FALSE)
    }
    else {
        balance.distribution <- balance.distribution + ggplot2::geom_line(alpha = 0.5) + 
            ggplot2::guides(size = FALSE)
    }
    if (boxplot.split) {
        balance.distribution <- balance.distribution + ggplot2::geom_boxplot(ggplot2::aes_string(group = "Group", 
            col = "n.group"), outlier.size = size.pt/2, position = ggplot2::position_dodge(0.8), 
            alpha = 0.5)
    }
    else {
        balance.distribution <- balance.distribution + ggplot2::geom_point(ggplot2::aes_string(col = "n.group"), 
            size = size.pt/2, position = ggplot2::position_dodge(0.8), 
            alpha = 0.5)
    }
    grid::grid.newpage()
    grob <- cbind(ggplot2::ggplotGrob(balance.partition), ggplot2::ggplotGrob(balance.distribution), 
        size = "first")
    grid::grid.draw(grob)
    res <- methods::new("bplot")
    res@balance.partition <- balance.partition
    res@balance.distribution <- balance.distribution
    return(res)
}                         

                        

#做一个table来存储otu_file, sample_file, tax_file,
distal_DBA <- function(otu_file, sample_file, #taxtable,
                    #    phynotype="Disease.MESH.ID", 
                    #    tax_level  = 'Phylum', # Domain  Phylum  Class   Order   Family  Genus   Species 
                       limitsize = FALSE, #for ggsave
                       dlimits=0.0001
                       ){
    sep <- ifelse(endsWith(otu_file,'csv'),',','\t')
    otu_table <- read.table(otu_file, header=T, row.names=1,sep=sep) #%>% t()    
    sep <- ifelse(endsWith(sample_file,'csv'),',','\t')
    sample_data_withIndex <- read.table(sample_file, header=T, row.names=1,comment.char='', sep=sep) #%>%
                #    .[eval(phynotype)] 
    sample_data <- sample_data_withIndex[,1]
    # tax_group <- taxtable[rownames(otu_table),] %>% .[, eval(tax_level)]
    otu_table <- t(otu_table) 

#     tax_group <- read.table(tax_file, header=T, row.names=1,comment.char='', sep=sep) %>%
#                    .[, eval(tax_level)] 
    # write.csv(sample_data_withIndex,file=paste0(sample_file,'.diseaseId.csv'))
    
                    
    dl <- rep(dlimits, ncol(otu_table)) 
    #Numeric vector or matrix of detection limits/thresholds 0被替换为小于dl的值
    if (any(otu_table==0)) {
        otu_table <- multRepl(otu_table, label=0, dl=dl, closure=100) #百万分之1
        #closure=10^6 if ppm or closure=100 if percentages
        # imp.missing = TRUE 会把0当成缺失值，补充为可能较大的值
    }

    sbp <- sbp.fromADBA(otu_table, sample_data) #all balances
    sbp <- sbp.subset(sbp) # subset D-2,D-3 part balances  
    #分子__分母1_分母2 or 分子1_分子2__分母
    #输出sbp矩阵要保存下来，用于物种相对丰度的balance特征提取（balance.fromSBP），
    # 然后用模型的最优特征集进行预测（z特征取一部分z）。
    colnames(sbp)<-str_replace(colnames(sbp),'^X','') #去掉数字前面的X
    otu_file <- str_sub(otu_file, end=-5) #去掉后缀
    write.csv(sbp,file=paste0(otu_file,'.sbp.csv')) 
    
    # 可视化balance
    balance_plot <- balance_plot_yang(otu_table, sbp, 
                                #  d.group=tax_group, #phylum 取消
                                 n.group=sample_data, #物种分组(只能小于8个，之分析门水平)和样品分组 ,颜色不够用
                                 boxplot.split = TRUE,
                                 weigh.var = TRUE, size.text = 10, size.pt = 2
                                )
    
    # save(balance_plot, file = paste0(otu_file,'.balance_plot.RData')) ##后续优化绘图算法 
    balance_ggar <- ggarrange(balance_plot@balance.partition + 
                              rremove("x.text") + rremove("x.ticks") +
                              theme(legend.direction ="vertical",legend.position = "top",
                                    axis.text.x = element_text(face = "italic") #genus species italic
                                   ) + #horizontal 
                              guides(col = guide_legend(ncol = 2, byrow = F,keywidth =0., keyheight = 0.6)),
                              
                              balance_plot@balance.distribution + rremove("y.text")+
                              theme(legend.direction ="vertical",legend.position = "top") +
                              guides(col = guide_legend(reverse = T)), 
                              
                              nrow = 1, align = "h",
                              widths = c(5, 2)
                             )
    
    lt10 <- 10*pound2inch*Linesp*length(balance_plot@balance.partition$data$Part)/2
    # lt10 <- 10*pound2inch*Linesp* dim(top20_sbp)
    height <- ifelse(lt10>10, lt10, 10)
    # ggsave(plot=balance_ggar, filename=paste0(otu_file,'.balance.partition_distribution.svg'), 
        #    width = 10, height = height, units = "in", limitsize = limitsize)limitsize = limitsize
    ggsave(plot=balance_ggar, filename=paste0(otu_file,'.balance.partition_distribution.png'), 
           width = 10, height = height, units = "in")
    ggsave(plot=balance_ggar, filename=paste0(otu_file,'.balance.partition_distribution.pdf'), 
           width = 10, height = height, units = "in" )
    
    write.csv(balance_plot@balance.partition$data, file=paste0(otu_file,'.balance.partition.csv')) #存储绘图数据用于后续图标优化
    write.csv(balance_plot@balance.distribution$data, file=paste0(otu_file,'.balance.distribution.csv'))

    z <- balance.fromSBP( #Compute Balances from an SBP Matrix
        x = otu_table,
        y = sbp
    )
    rownames(z) <- rownames(sample_data_withIndex)# for sklearn logistic regression
    write.csv(z,file=paste0(otu_file,'.balance.csv')) 
    
    # rownames(z) <- sample_data #for weka transform
    # write.csv(z,file=paste0(otu_file,'.balance.weka.csv')) 
    log <- 'run distal_DBA finished'
    return(log)
}

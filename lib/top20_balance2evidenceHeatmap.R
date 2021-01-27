suppressWarnings(suppressMessages(library(ComplexHeatmap)))
suppressWarnings(suppressMessages(library(circlize)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(RColorBrewer)))

Top20balanceMicroPhenoDB_Heatmap <- function (
  # otuf="otu_table.xls", 
  # groupmff='sample_data.xls.diseaseId.csv',
  lrcoefile='balance.csv.LogisticRegression.hypertuned.coef.csv',
  sbpf='otu_table.xls.balance.sbp.csv',
  MicroPhenoDBfile=NA, #判断taxid是数字的话用ncbi版本，物种名就用老版本MicroPhenoDB.xls
  Mesh2diseasef=NA,
  topn=20,
  outdir='./'
){
#   otuf='/opt/services/djangoapp/src/mAML/balance/dashApps/utils/testoutdir/otu_table.xls'
#   groupmff='/opt/services/djangoapp/src/mAML/balance/dashApps/utils/testoutdir/sample_data.xls.diseaseId.csv'
#   sbpf='/opt/services/djangoapp/src/mAML/balance/dashApps/utils/testoutdir/otu_table.xls.balance.sbp.csv'
#   lrcoefile='/opt/services/djangoapp/src/mAML/balance/dashApps/utils/testoutdir/otu_table.xls.balance.csv.LogisticRegression.hypertuned.coef.newids.csv'
#   MicroPhenoDBfile='/opt/services/djangoapp/src/mAML/balance/static/MicroPhenoDB/MicroPhenoDB.GI.ncbitaxid.xls'
#   Mesh2diseasef='/opt/services/djangoapp/src/mAML/balance/dashApps/utils/GMrepoMeshIDdiseaseUniq.xls'
  pound2inch = 1/72 #字体大小到长度的转换# 1.5倍行距的图片，  
  Linesp=1.5 #行间距
  
  # Microbe Disease Score
  # 543     Ankylosing spondylitis  -0.12776538018793834
  # 216816  Non-alcoholic fatty liver disease       -0.2148012081746108
  # 433284  Type i diabetes mellitus        -0.3252255354848429

  ## readin data
  # sep <- ifelse(endsWith(otuf,'csv'),',','\t')
  # otu <- read.csv(otuf, row.names = 1,sep=sep,check.names=F) # rows are taxids 
  # rownames(otu) <- str_replace_all(row.names(otu), "[\\[\\]]", "")

  # sep <- ifelse(endsWith(groupmff,'csv'),',','\t')
  # groupmf <- read.csv(groupmff, row.names = 1, sep=sep,check.names=F)

  sep <- ifelse(endsWith(lrcoefile,'csv'),',','\t')
  D015212_lrcoefall <- read.csv(lrcoefile,header=T,sep=sep,check.names=F) #lrcoef  
    #   dataset77.Collected_in_the_outer__Homogenised.xls.balance.csv.LogisticRegression.hypertuned.coef.csv
    #   ,homogenised_(LH)
    # Desulfotomaculum.ruminis__Prevotella.stercorea,0.0
    # Faecalibacterium.prausnitzii__Coprococcus.eutactus,0.0

  

  topn <- ifelse(nrow(D015212_lrcoefall)>topn, topn, nrow(D015212_lrcoefall))  
  D015212_lrcoef <- D015212_lrcoefall[order(-abs(D015212_lrcoefall[,2])),][1:topn,]
  rownames(D015212_lrcoef)<- D015212_lrcoef[,1]
  otufname <- str_replace(basename(sbpf),'.balance.sbp.csv','')
  otuf <- paste0(outdir,otufname) #output dir .balance.sbp.csv otu_table.xls.balance.top32_relAbundance.pdf
  write.csv(D015212_lrcoef,file=paste0(otuf,'.balance.top',topn,'.LRcoefficient.csv'), row.names = F)
  diseasename <- colnames(D015212_lrcoef)[2]
  colnames(D015212_lrcoef)<-c("BalanceID", "LRcoef")
  # colnames(D015212_lrcoef)[1] <- "BalanceID"
#   evidencefile =   paste0(otuf,'.balance.partition.csv') #
#   ncbi2spe_pnylumn_evidence <- read.csv(evidencefile,row.names=1) # .balance.partition.csv
# # "","Part","BalanceID","Component","Group","d.group"
# # "139","1","Arcanobacterium__Dethiosulfovibrio_Gramella","Arcanobacterium","Arcanobacterium__Dethiosulfovibrio_Gramella1","p__Actinobacteria"
#   colnames(ncbi2spe_pnylumn_evidence) <- c("Part","BalanceID","Species","Group","Phylum")
#   ncbi2spe_pnylumn_evidence$Species <- str_replace_all(ncbi2spe_pnylumn_evidence$Species, "[\\[\\]]", "")
#   D015212_ncbi_phylumn_evidence <- semi_join(ncbi2spe_pnylumn_evidence, D015212_lrcoef, by = "BalanceID")  
#     #balance degree changeg2 points %>%
#   D015212_ncbi_phylumn_evidence$Species <- str_replace_all(D015212_ncbi_phylumn_evidence$Species, "[\\[\\]]", "")
#   D015212_ncbi_phylumn_evidence[,1:3] %>%  
#     spread(Species,Part)  %>% # d[,c(2,3,4)] %>% spread(b,fac) 
#     column_to_rownames(var = 'BalanceID') %>%  t() ->  top20_sbp #以这个为顺序依据 #读取sbp 矩阵
  # top20_sbp[is.na(top20_sbp)] <- 0
  #与top20_sbp顺序一致 [] 都是 "."


  # ncbi2spe_pnylumn_evidence[,c(3,5)] %>% distinct() -> uniqEvid #species contain "[]"
  # rownames(uniqEvid) <- uniqEvid$Species 
  # uniqEvid<-uniqEvid[rownames(top20_sbp),] #species 都带'[]'
  
  

  sep <- ifelse(endsWith(sbpf,'csv'),',','\t')
  sbp <- read.csv(sbpf,header=T,sep=sep,check.names=F,row.names=1)
  sbp[,rownames(D015212_lrcoef)] %>% .[rowSums(abs(.))>0,] ->top20_sbp
  D015212_lrcoef <- D015212_lrcoef[colnames(top20_sbp),]


  map(rownames(top20_sbp), str_length) %>% as.numeric() %>%  max() -> longestSpeName
  longestSpeName_len <- 10*pound2inch*longestSpeName

  map(colnames(top20_sbp), str_length) %>% as.numeric() %>%  max() -> longestBalanceName
  longestBalanceName_len <- 10*pound2inch*longestBalanceName 

  LRcoef_col_fun = colorRamp2(c(0,range(summary(D015212_lrcoef$LRcoef))),
                              c("white","#00BFC4", "red")) 
  topha = HeatmapAnnotation(
      LR_coefficient = anno_simple(D015212_lrcoef$LRcoef, col = LRcoef_col_fun),
      annotation_name_side = "right")

  lgd_LRcoef  <-  Legend(title = "LR coefficient", col = LRcoef_col_fun, 
                        at = c(0.0, range(D015212_lrcoef$LRcoef)),
                        labels = c("0.0", round(range(D015212_lrcoef$LRcoef),1)))
  
  colnames(top20_sbp) %>% str_replace_all("^\\.+|\\.+$",'') %>% #trim space
    str_replace_all("_\\.+|\\.+_$",'_') %>% #trim _ space
    str_replace_all("\\.+",' ') -> colnames(top20_sbp)

  ht1 <- Heatmap(as.matrix(top20_sbp), col=c("white","white"),  na_col = "white",
        rect_gp = gpar(col= "grey"),
        layer_fun = function(j, i, x, y, width, height, fill) {
          v = pindex(as.matrix(top20_sbp), i, j)
          l = v > 0
          p = v < 0
          grid.points(x[l], y[l], pch = 1, size = unit(2, "mm"))
          grid.points(x[p], y[p], pch = 2, size = unit(2, "mm"))
          },
        
        
        cluster_rows = T,show_row_dend = F,
        cluster_columns = F,show_column_dend  = F,
        bottom_annotation = topha,
        # left_annotation = leftha,
        # right_annotation = rightha,
        # row_order = order(uniqEvid$Phylum),
        column_order = order(D015212_lrcoef$LRcoef,decreasing = T),
        row_names_side = "left",
        # row_labels = uniqEvid$Species, #[替换没了
        row_names_max_width = unit(longestSpeName_len,"inch"), #unit(6, "cm"),
        column_names_max_height = unit(0.5*longestBalanceName,"inch"), 
        row_names_gp = gpar(fontface = "italic"),#fontsize = 12,
        column_names_side = "top",
        show_column_names = TRUE, #balance name too long to show
        column_names_rot = 45,
        column_names_gp = gpar(fontface = "italic"),
        
        show_heatmap_legend = F,
        width = dim(top20_sbp)[2]
        # heatmap_width = dim(top20_sbp)[2] + 2*longestSpeName_len fontsize = 14, fontface = "bold"
        # height = 12,
        )
  # width and height only control the width/height of the heamtap body
  lgd_partition  <-  Legend(title = "Partition", 
                            type = "points", pch=1:2,
                            labels = c(1,-1))

  ## save heatmap pdf and svg

  lt10 <- 10*pound2inch*Linesp*dim(top20_sbp)[1] +  0.1*longestBalanceName # 1/2 longest colnames
  height <- ifelse(lt10>10, lt10, 10) #不能太小
  wLinesp=1.5 #column行间距


  if(!is.na(MicroPhenoDBfile)){ #upload evidence or use microPhenoDB
     
  
    sep <- ifelse(endsWith(MicroPhenoDBfile,'csv'),',','\t')
    MicroPhenoDB <- read.csv(MicroPhenoDBfile,sep=sep,check.names=F) #core_table.GI.txt
    colnames(MicroPhenoDB) <- c('Microbe', 'Disease', 'Score')
    
    MicroPhenoDB  %>% 
      spread(Disease,Score) %>% 
      column_to_rownames("Microbe")  -> MicroPhenoDB_mat


    ## MicroPhenoDB evidence for D015212
      # Marinobacter__Aquifex_Candidatus.Carsonella1 R会把空格变成"."
    # uniqEvid$Species <- str_replace_all(uniqEvid$Species, "[\\[\\]]", "")
    MicroPhenoDB_mat[as.character(rownames(top20_sbp)),] -> D015212_MicroPhenoDB_mat
    #NA NA NA <NA> <NA> 没有证据的id填充NA
    D015212_MicroPhenoDB_mat[is.na(D015212_MicroPhenoDB_mat)] <- 0
    rownames(D015212_MicroPhenoDB_mat) <- rownames(top20_sbp)
      
      
    #过滤都是0的列
    D015212_MicroPhenoDB_mat[,map(abs(D015212_MicroPhenoDB_mat), max)>0] -> D015212_MicroPhenoDB_mat_del0 
    
    write.csv(D015212_MicroPhenoDB_mat_del0,file=paste0(otuf,'.balance.top',topn,'.MicroPhenoDB.csv') )

    if (str_starts(diseasename,'D0')){
      sep <- ifelse(endsWith(Mesh2diseasef,'csv'),',','\t')  
      Mesh2disease  <- read.csv(Mesh2diseasef,sep=sep,check.names=F,row.names=1,header=F)
      if (diseasename %in% rownames(Mesh2disease)){ 
        diseasename=Mesh2disease[diseasename,1] 
      }
    }

    # evidence
    dcolnames <- tolower(str_replace(colnames(D015212_MicroPhenoDB_mat_del0),"\'s",""))
    dname <- tolower(diseasename)  
    diseasename=colnames(D015212_MicroPhenoDB_mat_del0)[which(dcolnames==dname)]
    diseaseScore <- D015212_MicroPhenoDB_mat_del0[,diseasename]
    diseaseScores <- round(diseaseScore[diseaseScore!=0],2)
    rownames <- rownames(D015212_MicroPhenoDB_mat_del0)[diseaseScore!=0]
    speciesEvidence <- data.frame(Species=rownames,diseaseCorr=diseaseScores)
    filename <- paste0(otuf,'.balance.top',topn,'.evidence.csv')
    write.csv(speciesEvidence,filename, row.names = F,quote=F)

    numers <- str_split(D015212_lrcoef[,1],'_+',simplify = T)
    
    spe2evidScore <- function(spe){
      score=ifelse(spe %in% speciesEvidence$Species, 
            speciesEvidence[speciesEvidence$Species==spe,]$diseaseCorr,
            NA
      )
    }
    numers[,1] %>% map(spe2evidScore) %>% unlist() -> D015212_lrcoef$NumeratorE
    numers[,2] %>% map(spe2evidScore) %>% unlist() -> D015212_lrcoef$Denominator1E
    numers[,3] %>% map(spe2evidScore) %>% unlist() -> D015212_lrcoef$Denominator2E
    D015212_lrcoef[is.na(D015212_lrcoef)]<-0 #NA to 0
    speciesInfer <- function(BalanceID, LRcoef,NumeratorE,Denominator1E,Denominator2E){
      nums <- str_split(BalanceID,'_+',simplify = T)
      #      [,1]   [,2]   
      # [1,] "1531" "45851"
      infer='-'
      if(dim(nums)[2]==2){ #2part balance
        if(LRcoef>0){
          if(all(NumeratorE<0, Denominator1E==0)){
            infer=str_c(nums[1,2],"_-")
          }else if(all(Denominator1E>0, NumeratorE==0)){
            infer=str_c(nums[1,1],"_+")
          }
        }else if(LRcoef<0){
          if(all(NumeratorE>0, Denominator1E==0)){
            infer=str_c(nums[1,2],"_+")
          }else if(all(Denominator1E<0, NumeratorE==0)){
            infer=str_c(nums[1,1],"_-")
          }
        }
      }else if(dim(nums)[2]>2){ #3part balance
        ds<-c(Denominator1E, Denominator2E)
        if(LRcoef>0){         
          if(NumeratorE<0 & all(any(ds>0), 0 %in% ds)){
            infer=str_c(nums[1,ds==0],"_-")
          }else if(all(ds>0)){
            infer=str_c(nums[1,1],"_+")
          }
        }else if(LRcoef<0){
          if(NumeratorE>0 & all(any(ds<0), 0 %in% ds) ){
            infer=str_c(nums[1,inferds==0],"_+")
          }else if(all(ds<0)){
            infer=str_c(nums[1,1],"_-")
          }
        }
      }
      infer
    }
    D015212_lrcoef %>% pmap(speciesInfer) %>% unlist() -> D015212_lrcoef$Species_infer 
    # D015212_lrcoef$Species_infer[is.na(D015212_lrcoef$Species_infer)]<-'-'
    D015212_lrcoef[D015212_lrcoef==0]<-'-'
    filename <- paste0(otuf,'.balance.top',topn,'.species_infer.csv')
    write.csv(D015212_lrcoef,filename, row.names = F,quote=F)


    D015212_MicroPhenoDB_mat_del0 %>% 
      round(1) -> D015212_MicroPhenoDB #abs(score) > 0.3 证据可靠，这个疾病需要保留
    # D015212_MicroPhenoDB_mat %>% round(1) -> D015212_MicroPhenoDB

    # "Pervasive developmental disorder - not otherwise specified" #too long name 
    # colnames(D015212_MicroPhenoDB)[7]  <- "Pervasive developmental disorder" ## 编号可能会变
    # colnames(D015212_MicroPhenoDB)[map(colnames(D015212_MicroPhenoDB),str_length) >35] <- "Pervasive developmental disorder"

    #cor_diseases <- c("Crohn's disease", "Ulcerative colitis")
    coln <- colnames(D015212_MicroPhenoDB)
    colnames(D015212_MicroPhenoDB)[coln =="Pervasive developmental disorder - not otherwise specified"] <- "Pervasive developmental disorder"

    #reindex <- c(which(coln %in% cor_diseases), which(! coln %in% cor_diseases)) #cor_diseases in first columns
    #D015212_MicroPhenoDB <- D015212_MicroPhenoDB[,reindex]



    ## MicroPhenoDB evidence heatmap for D015212

    MicroPhenoDB_col_fun = colorRamp2(c(0,range(D015212_MicroPhenoDB)),
                                      c("white","#00BFC4", "red")) 
    rangeHas3 <- c(-0.3,0.3)[range(abs(D015212_MicroPhenoDB))>0.3]
    lgd_MicroPhenoDB  <-  Legend(title = "Association score", col = MicroPhenoDB_col_fun, 
                          at = c(0.0, range(D015212_MicroPhenoDB), rangeHas3),
                          labels = c("0.0", round(range(D015212_MicroPhenoDB),1), rangeHas3))
    

    map(as.character(colnames(D015212_MicroPhenoDB)),str_length) %>% as.numeric() %>%  max() -> longestName
    longestName_len <- 10*pound2inch*longestName  
    
    

    ht2 <- Heatmap(as.matrix(D015212_MicroPhenoDB), name="Association score",
                na_col = "white",
                col=MicroPhenoDB_col_fun,
          rect_gp = gpar(col= "grey"),
          
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(D015212_MicroPhenoDB[i,j]!= 0) {
              grid.text(             
                D015212_MicroPhenoDB[i,j],
                x = x, y = y,
                gp = gpar(fontsize = 8),
                just = c("center","center"),check.overlap = T
                )
              }
            },                 
          
          cluster_rows = F,show_row_dend = F,
          cluster_columns = F,show_column_dend  = F,

          row_names_side = "right",
          show_row_names = F,
          column_names_side = "top",
          column_names_rot = 45,
          column_names_max_height = unit(longestName_len, "inch"), 
          show_heatmap_legend = F,
          width = dim(D015212_MicroPhenoDB)[2]
          )

    width <-  10*pound2inch*wLinesp*(dim(top20_sbp)[2] + dim(D015212_MicroPhenoDB)[2]) +  longestName_len + longestSpeName_len
      
      
    ht_list = ht1 + ht2 #设置heatmap 方块放个的大小，就能统一list里热土的比例
    pd_list <- packLegend(lgd_partition, lgd_MicroPhenoDB, lgd_LRcoef,
                          # max_height = unit(10, "cm"), lgd_Phylum, 
                          row_gap = unit(1, "cm"))


    filename <- paste0(otuf,'.balance.top',topn,'.evidence')
    pdf(str_c(filename,'.pdf'), width=width, height=height)
    draw(ht_list, #row_title = "top20Balance species", 
        row_title_gp = gpar(col = "black"),
        annotation_legend_list = pd_list,
        merge_legend = TRUE,
        # heatmap_legend_side = "bottom", 
        # annotation_legend_side = "bottom",
        # column_title = "Inflammatory Bowel Diseases (D015212)", 
        column_title_gp = gpar(fontsize = 16))

    dev.off()

    png(str_c(filename,'.png'), width=width, height=height, units = "in", res = 100)
    # par(mar=c(5.1, 28.1, 4.1, 12), xpd=F) 
    draw(ht_list, #row_title = "top20Balance species", 
        row_title_gp = gpar(col = "black"),
        annotation_legend_list = pd_list,
        # list(lgd_partition,  lgd_LRcoef, lgd_BalanceDegree, lgd_MicroPhenoDB),
        merge_legend = TRUE,
        # heatmap_legend_side = "bottom", 
        # annotation_legend_side = "bottom",
        # column_title = "Inflammatory Bowel Diseases (D015212)", 
        column_title_gp = gpar(fontsize = 16))

    dev.off()

  }else{
    
    width <-  10*pound2inch*wLinesp*(dim(top20_sbp)[2]) + longestSpeName_len
      
      
    ht_list = ht1
    pd_list <- packLegend(lgd_partition, lgd_LRcoef,
                          # max_height = unit(10, "cm"), lgd_Phylum, 
                          row_gap = unit(1, "cm"))

    D015212_lrcoef[D015212_lrcoef==0]<-'-'
    filename <- paste0(otuf,'.balance.top',topn,'.csv')
    write.csv(D015212_lrcoef,filename, row.names = F,quote=F)


    filename <- paste0(otuf,'.balance.top',topn)
    pdf(str_c(filename,'.pdf'), width=width, height=height)
    draw(ht_list, #row_title = "top20Balance species", 
        row_title_gp = gpar(col = "black"),
        annotation_legend_list = pd_list,
        merge_legend = TRUE,
        # heatmap_legend_side = "bottom", 
        # annotation_legend_side = "bottom",
        # column_title = "Inflammatory Bowel Diseases (D015212)", 
        column_title_gp = gpar(fontsize = 16))

    dev.off()

    svg(str_c(filename,'.svg'), width=width, height=height)
    # par(mar=c(5.1, 28.1, 4.1, 12), xpd=F) 
    draw(ht_list, #row_title = "top20Balance species", 
        row_title_gp = gpar(col = "black"),
        annotation_legend_list = pd_list,
        # list(lgd_partition,  lgd_LRcoef, lgd_BalanceDegree, lgd_MicroPhenoDB),
        merge_legend = TRUE,
        # heatmap_legend_side = "bottom", 
        # annotation_legend_side = "bottom",
        # column_title = "Inflammatory Bowel Diseases (D015212)", 
        column_title_gp = gpar(fontsize = 16))

    dev.off()


  }

  log <- paste0("finish biomarker discovery: ",lrcoefile)
  return(log)


}

library(devtools)
install_github("Coolgenome/iTALK", build_vignettes = TRUE)
install.packages("ps")
install.packages("scater")
BiocManager::install("scater")

devtools::install_local("G://R-3.6.2/library/iTALK-master.zip")
devtools::install_github('arc85/celltalker')
install.packages("G://google Downloads/iTALK-master/",repos=NULL,type="source")

library(iTALK)
library(Seurat)

#不愉快的开始，iTALK分析

#清除环境
rm(list = ls()) 
options(stringsAsFactors = F) 
options(digits = 2)




#设置路径
setwd("G:/Project/NTC_CD8/肿瘤组织样本CD8/")

genecount <- read.csv("Genecount_TTC.csv",header = T ,row.names = 1)



genemeta <- read.csv("TTC_monocleMeta.csv", header = T, row.names = 1)


TTC <- CreateSeuratObject(genecount,
                          meta.data = genemeta,
                          project = "TTC",
                          min.cells = 5,
                          min.features  = 200)


##准备iTALK输入文件
# 输入文件格式为数据框，行为细胞名称，列为基因名和细胞注释信息，基因名命名的列值为表达数据。
# 细胞注释信息必选“cell_type”，值为细胞类型注释信息；可选“compare_group”，值为细胞的样本分组信息。

cell.meta <- subset(TTC@meta.data, select=c("sampleType","pseudotime"))

names(cell.meta) <- c("compare_group", "cell_type")

cell.expr <- data.frame(t(as.matrix(TTC@assays$RNA@counts)), check.names=F)
data.italk <- merge(cell.expr, cell.meta, by=0)
rownames(data.italk) <- data.italk$Row.names
data.italk <- data.italk[,-1]

#data.italk[1:5,c(1:2,18027:18028)]

##设置绘图颜色和其他变量
mycolor <- c("#92D0E6","#F5949B","#E11E24","#FBB96F","#007AB8","#A2D184","#00A04E","#BC5627","#0080BD","#EBC379","#A74D9D")
cell_type <- unique(data.italk$cell_type)
cell_col <- structure(mycolor[1:length(cell_type)], names=cell_type)
#通讯类型变量
comm_list<-c('growth factor','other','cytokine','checkpoint')
#圈图展示配体-受体对的数量
PN=20

##==细胞通讯关系概览==##
dir.create('iTALK/Overview')
setwd("iTALK/Overview")

#实际上一个样本之内的细胞才能通讯，这里提取一组样本分析是为了克服单细胞数据稀疏性造成的误差
data1 <- subset(data.italk, subset=data.italk$compare_group=="TTC")
##寻找高表达的配体受体基因,top_genes=50代表提取平均表达量前50%的基因
highly_exprs_genes <- rawParse(data1, top_genes=50, stats="mean")
res<-NULL
for(comm_type in comm_list){
  #comm_type = 'cytokine'   #测试循环代码的临时变量
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  res_cat <- res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs, decreasing=T),]
  write.csv(res_cat, paste0('LRpairs_Overview_',comm_type,'.csv'))
  pdf(paste0('LRpairs_Overview_',comm_type,'.pdf'), width=6, height=7)
  #png(paste0('LRpairs_Overview_',comm_type,'_net.png'), width=600, height=650)
  #绘制细胞通讯关系网络图
  NetView(res_cat, col=cell_col, vertex.label.cex=1.2, edge.label.cex=0.9, 
          vertex.size=30, arrow.width=3, edge.max.width=10, margin = 0.2)
  dev.off()
  #绘制topPN的配体-受体圈图
  #png(paste0('LRpairs_Overview_',comm_type,'_circ.png'), width=600, height=650)
  pdf(paste0('LRpairs_Overview_',comm_type,'.pdf'), width=6, height=7)
  LRPlot(res_cat[1:PN,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:PN],link.arr.width=res_cat$cell_to_mean_exprs[1:PN], text.vjust = "0.35cm")
  dev.off()
  res<-rbind(res,res_cat)
}

res <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs, decreasing=T),]

index_other <- res$comm_type == "other"
res_no_other <- res[!index_other,]

write.csv(res, 'LRpairs_Overview.csv')
res.top <- res_no_other[1:100,]

res.top <- res.top[order(res.top$cell_from_mean_exprs*res.top$cell_to_mean_exprs, decreasing=F),]

png('LRpairs_Overview_net.png', width=600, height=650)
NetView(res, col=cell_col, vertex.label.cex=1.2, edge.label.cex=0.9, 
        vertex.size=30, arrow.width=3, edge.max.width=10, margin = 0.2)
dev.off()
png('LRpairs_Overview_circ.png', width=600, height=650)

LRPlot(res.top, datatype='mean count', link.arr.lwd=res.top$cell_from_mean_exprs,
       cell_col=cell_col, link.arr.width=res.top$cell_to_mean_exprs)
dev.off()
setwd("~/project/2020/2007_10xDemo2")

res.top <- 
LRPlot(res.top, datatype='mean count', link.arr.lwd=res.top$cell_from_mean_exprs,
       cell_col=cell_col, link.arr.width=res.top$cell_to_mean_exprs)


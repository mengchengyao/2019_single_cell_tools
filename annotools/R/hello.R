# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}


#' Title go_listfile
#'
#' @param species the ref gene annotation from annorad loacal database(Go Term)
#' @param indir  the database indir
#'
#' @return go annotation file name
#' @export no export
#'
#' @examples go_list<-go_listfile("Mus_musculus.GRCm38.90.chr",indir='/annoroad/bioinfo/PMO/database/config_database')
go_listfile<-function(species,indir="/annoroad/bioinfo/PMO/database/config_database"){
  annofile=paste(indir,paste(species,".txt",sep=""),sep="/")
  print(paste("annotion file is :",annofile,sep=" "))
  ini.list <- fread(annofile,header=F,sep="=",stringsAsFactors = FALSE,quote = "",data.table = F,fill=TRUE)
  go_gile<-ini.list[ini.list$V1=="GO_annotate",]$V2
  return(go_gile)
}
#' Title go_gene
#'
#' @param infile  go annotation file name
#' @param outfile  Gene GO txt file,just only two col, and the first col is gene(ensemble), the second col is GO term name
#' @param python  absulately python path to deal go annotation file
#' @param script script to deal go annotation file
#'
#' @return data.fram, the first col is Go term names and the second col is Gene name
#' @export term2gene.xls
#'
#' @examples term2gene<-go_gene(infile,paste(outdir,"term2gene.xls",sep="/"),script=script,python=python)[,c("GO","Gene")]
go_gene<-function(infile,outfile,python="/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/python",script="/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RNA/sRNA/ngs_bioinfo/TET-1487/yaomengcheng/Analysis-test/clusterProfiler/gene.go.py"){
  CMD<-paste(python,script,infile,outfile,sep=" ")
  print(paste("开始GENE2GO文件转换，时间为：",Sys.time(),sep=''))
  system(CMD,intern = TRUE,ignore.stderr = TRUE)
  print(paste("开始GENE2GO文件读取，时间为：",Sys.time(),sep=''))
  tmp<-fread(outfile,sep="\t",header=T,stringsAsFactors = FALSE,quote = "",data.table = F)
  return(tmp)
}
#' Title seurat_marker
#'
#' @param markerfile  a markerfile from seurat diff for cluster
#' @param sep sep for marker file ,default=","
#' @param colname which col for group,defaul="cluster"
#'
#' @return list,
#' @export
#'
#' @examples
seurat_marker<-function(markerfile,sep=",",colname="cluster"){
  marker_list<-list()
  print(paste("开始读取markerfile文件，时间为：",Sys.time(),sep=" "))
  marker<-fread(markerfile,header=T,sep=sep,stringsAsFactors = FALSE,quote = "",data.table = F)
  for (i in unique(marker[[colname]])){
    marker_list[[as.character(i)]]<-as.vector(marker[marker[[colname]]==i,]$gene)
  }
  #names(marker_list)<-unique(marker[[colname]])
  print(paste("完成markerfile转换成marker_list，时间为：",Sys.time(),sep=" "))
  return (marker_list)
}
#' Title get gene symbol from ensemble name
#'
#' @param genefile genes.tsv from cellranger output
#' @param term2gene  data.fram, the first col is Go term names and the second col is Gene name
#' @param sep genefile sep,default sep="\t"
#'
#' @return term2gene  data.fram, the first col is Go term names and the second col is Gene name(symbol)
#' @export
#'
#' @examples term2gene<-ensemble2symbol(genefile,term2gene)
ensemble2symbol<-function(genefile,term2gene,sep="\t"){
  print(paste("开始读取gens.tsv文件，时间为：",Sys.time(),sep=" "))
  genes<-read.table(genefile,sep=sep)
  colnames(genes)<-c("Gene","Gene_symbol")
  print(paste("开始读取匹配ensemble和symbol文件，时间为：",Sys.time(),sep=" "))
  term2gene2<-merge(term2gene,genes,by="Gene",all=TRUE)
  print(paste("删掉go或者symbol为na的行，时间为：",Sys.time(),sep=" "))
  term2gene2<-term2gene2[!is.na(term2gene2$GO)&!is.na(term2gene2$Gene_symbol),c("GO","Gene_symbol")]
  print(paste("完成ensemble和symbol的转换，时间为：",Sys.time(),sep=" "))
  return (term2gene2)
}

#' Title
#'
#' @param markerfile       markerfile from seurat diff output
#' @param outdir           output dir for output file
#' @param genefile         ensemble and symbol file from cellrange output
#' @param prefix           output's prefix
#' @param species          gene go annotation for species,must in list for annorad annotation list
#' @param indir            database indir
#' @param python           which python for dealing with go annodataion file
#' @param script           which script for dealing with go annodataion file
#' @param go.class         Go Term and Go Term namse file
#' @param sep              marker file sep
#' @param colname          marker for group colnames
#' @param pvalueCutoff     compareCluster analysis for param,pvalueCutoff
#' @param pAdjustMethod    compareCluster analysis for param,pAdjustMethod
#' @param qvalueCutoff     compareCluster analysis for param,qvalueCutoff
#'
#' @return                 S4 object for result, from clusterprofiler
#' @export
#'
#' @examples               x<-compare(markerfile,outdir,genefile,prefix=prefix)
compare<-function(markerfile,outdir,genefile,prefix="compare_cluster",species="Mus_musculus.GRCm38.90.chr",indir="/annoroad/bioinfo/PMO/database/config_database",python="/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/python",script="/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RNA/sRNA/ngs_bioinfo/TET-1487/yaomengcheng/Analysis-test/clusterProfiler/gene.go.py",go.class="/annoroad/bioinfo/PMO/database/GO/current/go.class",sep=",",colname="cluster",pvalueCutoff = 0.1, pAdjustMethod = "BH", qvalueCutoff = 1){
  go_list<-go_listfile(species,indir=indir)
  term2gene<-go_gene(go_list,paste(outdir,"term2gene.xls",sep="/"),script=script,python=python)[,c("GO","Gene")]
  term2name <- fread(go.class,header=T,sep="\t",stringsAsFactors = FALSE,quote = "",data.table = F)[,c("Accession","Term_name")]
  gene_list<-seurat_marker(markerfile,sep=sep,colname=colname)
  term2gene<-ensemble2symbol(genefile,term2gene)
  print(paste("开始进行compareCluster富集分析，时间为：",Sys.time(),sep=" "))
  x3 <- compareCluster(gene_list, fun='enricher',TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)
  print(paste("完成compareCluster富集分析，时间为：",Sys.time(),sep=" "))
  write.table(x3@compareClusterResult,paste(outdir,paste(prefix,"compareClusterResult.xls",sep="_"),sep="/"),quote=F,row.names=F,sep="\t")
  return(x3)

}

#' Title
#'
#' @param x                 S4 object for result  from clusterprofiler
#' @param outdir            output dir
#' @param prefix            outputfile prefix
#' @param w                 the width for pdf file
#' @param h                 the heigh for pdf file
#' @param showCategory      show top numbers of analysis pdf file
#'
#' @return
#' @export
#'
#' @examples              compare_plot(x,outdir,prefix,w=12,h=8,showCategory=20)
compare_plot<-function(x,outdir,prefix,w=12,h=8,showCategory=10){
  print(paste("开始画图，时间为：",Sys.time(),sep=" "))
  pdf(paste(outdir,paste(prefix,"compareClusterResult.pdf",sep="_"),sep="/"),w=w,h=h)
  p0<-dotplot(x, showCategory=showCategory,includeAll=TRUE)
  p1 <- emapplot(x)
  p2 <- emapplot(x,legend_n=2)
  p3 <- emapplot(x,pie="count")
  p4 <- emapplot(x,pie="count", pie_scale=1.5, layout="kk")
  p5<-cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
  print(p0)
  print(p3)
  dev.off()
  print(paste("输出图片文件为：",paste(outdir,paste(prefix,"compareClusterResult.pdf",sep="_"),sep="/"),sep=" "))
  print(paste("完成画图，时间为：",Sys.time(),sep=" "))
}

#' Title
#'
#' @param marker_list    marker list from seurat diff output
#' @param outdir         outdir
#' @param prefix         outfile prefix
#' @param species        species for analysis, must be in single,Mouse,Human,all
#' @param showCategory   show top numbers for output pdf file
#' @param databese       cellmarker annotation file dir
#' @param w              The width of output pdf
#' @param h              The height of output pdf
#'
#' @return               The compareCluster result of output
#' @export               pdf file for result
#'
#' @examples             gene_list<-seurat_marker(markerfile,sep=",",colname="cluster")
#' @examples             cellmarker(gene_list,outdir,species="Mouse",prefix="cell_marker")
cellmarker<-function(marker_list,outdir,prefix="cellmarker",species="Mouse",showCategory=20,sep="\t",databese="/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RNA/sRNA/ngs_bioinfo/TET-1487/yaomengcheng/Analysis-test/clusterProfiler/database",w=12,h=12){
  cellmarkerfile<-paste(databese,paste(species,"cell_markers.txt",sep="_"),sep="/")
  cell_markers <- vroom::vroom(cellmarkerfile) %>%
    tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
    dplyr::select(cellMarker, geneSymbol) %>%
    dplyr::mutate(geneSymbol = strsplit(geneSymbol, ', '))
  print(cell_markers)
  y <- compareCluster(marker_list, fun='enricher',TERM2GENE=cell_markers,minGSSize=1,pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1)
  pdf(paste(outdir,paste(prefix,"compareCluster_cell_markers.pdf",sep="_"),sep="/"),w=w,h=h)
  p1<-dotplot(y, showCategory=20,includeAll=TRUE)
  print(p1)
  dev.off()
  write.table(y@compareClusterResult,paste(outdir,paste(prefix,"compareCluster_cell_markers.xls",sep="_"),sep="/"),quote=F,sep="\t")
  return (y)
}


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
#' @param indir   the database indir
#' @param type    the type for file[GO_annotate or KEGG_annotate]
#'
#' @return go annotation file name
#' @importFrom data.table fread
#' @examples go_list<-go_listfile("Mus_musculus.GRCm38.90.chr",indir='/annoroad/bioinfo/PMO/database/config_database',type="GO_annotate")
go_listfile<-function(species,indir="/annoroad/bioinfo/PMO/database/config_database",type="GO_annotate"){
  annofile=paste(indir,paste(species,".txt",sep=""),sep="/")
  print(paste("annotion file is :",annofile,sep=" "))
  ini.list <- fread(annofile,header=F,sep="=",stringsAsFactors = FALSE,quote = "",data.table = F,fill=TRUE)
  go_gile<-ini.list[ini.list$V1==type,]$V2
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
#' @importFrom data.table fread
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
#' @importFrom data.table fread
#' @examples  seurat_marker(markerfile,sep=",",colname="cluster")
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
#' Title go.list or ko.list to list
#'
#' @param gmt.file
#'
#' @return list,
#' @export
#' @examples  ko2map<-go_Pathways('/annoroad/data1/bioinfo/PMO/Public/database/Public/Quarter/KEGG/current/data/ko2map/ko2map.xls')
go_Pathways<-function(gmt.file)
{
  pathwayLines <- strsplit(readLines(gmt.file), "\t|\\|")
  pathways <- lapply(pathwayLines, tail, -1)
  names(pathways) <- sapply(pathwayLines, head, 1)
  pathways
}

#' Title mkdirs funtion(if dir exist, do not del the dir)
#' @param outdir outdir for mkdirs
#' @param fp dirnames for mkdirs
#'
#' @return mkdirs and retunr none
#' @export
#' @examples  mkdirs(outdir,'tmp/')
mkdirs <- function(outdir,fp) {
  if(!file.exists(file.path(outdir,fp))) {
    #		mkdirs(dirname(fp))
    dir.create(file.path(outdir,fp))
  }else{
    print(paste(fp,"Dir already exists!",sep="     "))
    #unlink(file.path(outdir,fp), recursive=TRUE)
    #dir.create(file.path(outdir,fp))
  }
}

#' Title go_updata(not use config.ini annotation file, using result of blastx and updata by python script with the same of annorad methods)
#'
#' @param outdir which mkdir tmp dirs for output
#' @param go_gile which come from  config.ini annotation file, but not using this annotation file and provided dirnames for blastx results
#' @param python_YMC which python to use
#' @param anno which python script anno.py to use
#' @param update_go which update_go script to use
#' @param sprot.anno which sprot.anno file to use
#' @param go.class which go.class file to use
#' @param go.update which go.update file to use
#' @return go_updata.list  a file go_updata.list
#'
#' @examples go_updata(outdir,go_gile)
go_updata<-function(outdir,go_gile,python_YMC="/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/python",anno='/annoroad/data1/bioinfo/PMO/Public/Pipeline/Stable/RNA/Public_Module/current/Function_GO//anno.py',update_go="/annoroad/data1/bioinfo/PMO/Public/Pipeline/Stable/RNA/Public_Module/current/Function_GO///update_go.py",sprot.anno='/annoroad/data1/bioinfo/PMO/Public/database/Public/Quarter/Uniport/current/data/sprot.anno',go.class='/annoroad/data1/bioinfo/PMO/Public/database/Public/Quarter/GO/current/data/go.class',go.update='/annoroad/data1/bioinfo/PMO/Public/database/Public/Quarter/GO/current/data/go.update'){
  ###判断输出文件是否存在
  mkdirs(outdir,'tmp/')
  setwd(paste(outdir,'tmp/',sep='/'))
  CMD1<-paste(python_YMC ,' ',anno,' -i ',paste(dirname(go_gile),"/../gtf/uniprot/gene.uniq.blastx",sep=''),' -a ',sprot.anno,' -o ','go.list',' -c 4 -t GO ',sep="")
  CMD2<-paste(python_YMC,' ',update_go,' -c ',go.class,' -u ',go.update,' -i ','go.list',' -o ','go_updata.list',' > ','replaced.go',sep="")
  print(paste("开始go list生成，时间为：",Sys.time(),sep=''))
  system(CMD1,intern = TRUE,ignore.stderr = TRUE)
  print(paste("开始go updata，时间为：",Sys.time(),sep=''))
  system(CMD2,intern = TRUE,ignore.stderr = TRUE)
  print(paste("完成go updata，时间为：",Sys.time(),sep=''))
  #tmp<-fread(outfile,sep="\t",header=T,stringsAsFactors = FALSE,quote = "",data.table = F)
  return(paste(outdir,'/tmp//go_updata.list',sep=''))
}##

#' Title gene_ko_update(not use config.ini annotation file, using result of blastx and updata by python script with the same of annorad methods)
#'
#' @param outdir which mkdir tmp dirs for output
#' @param go_gile which come from  config.ini annotation file, but not using this annotation file and provided dirnames for blastx results
#' @param python_YMC which python to use
#' @param anno which python script anno.py to use
#' @param perl which perl to use
#' @param extract_ko which extract_ko script to use
#' @param sprot.anno which sprot.anno file to use
#' @return ko_updata.list a data.frame(the first colum is Gene and the second colum is ko)
#'
#' @examples gene_ko_update(outdir,ko_gile)
gene_ko_update<-function(outdir,go_gile,python_YMC="/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/python",anno='/annoroad/data1/bioinfo/PMO/Public/Pipeline/Stable/RNA/Public_Module/current/KEGG//anno.py',perl='/annoroad/share/software/install//perl-5.16.2/bin/perl',extract_ko="/annoroad/data1/bioinfo/PMO/Public/Pipeline/Stable/RNA/Public_Module/current/KEGG/extract_ko_from_bgi_result.pl",sprot.anno='/annoroad/data1/bioinfo/PMO/Public/database/Public/Quarter/Uniport/current/data/sprot.anno'){
  ###判断输出文件是否存在
  mkdirs(outdir,'tmp/')
  setwd(paste(outdir,'tmp/',sep='/'))
  CMD1<-paste(python_YMC ,' ',anno,' -i ',paste(dirname(go_gile),"/../gtf/uniprot/gene.uniq.blastx",sep=''),' -a ',sprot.anno,' -o ','ko.list',' -c 3 -t KO ',sep="")
  #CMD2<-paste(python,' ',update_go,' -c ',go.class,' -u ',go.update,' -i ','go.list',' -o ','go_updata.list',' > ','replaced.go',sep="")
  print(paste("开始gene ko 生成，时间为：",Sys.time(),sep=''))
  system(CMD1,intern = TRUE,ignore.stderr = TRUE)
  #print(paste("开始go updata，时间为：",Sys.time(),sep=''))
  #system(CMD2,intern = TRUE,ignore.stderr = TRUE)
  print(paste("完成gene ko，时间为：",Sys.time(),sep=''))
  tmp<-fread(paste(outdir,'/tmp//ko.list',sep=''),sep="\t",header=F,stringsAsFactors = FALSE,quote = "",data.table = F)
  colnames(tmp)<-c('Gene','KO')
  return(tmp)
}##


#' Title id2symbol(just for apply funciton)
#'
#' @param X A vector of gene list
#' @param gene_id and gene_symbol data.frame
#' @return a gene list with Gene_symbol vector
#'
#' @examples id2symbol(gene_list)
id2symbol<-function(x,genes){
  return(as.vector(genes[genes$Gene %in%x,]$Gene_symbol))
}


#' Title kegg_enrich_local_compareCluster
#'
#' @param gene_list              a gene_list for enrichment, the gene list of length range from 1-max
#' @param outdir                 output dir for output file
#' @param ensemble_symbol_file   ensemble and symbol file from cellrange output
#' @param ko_color               output geneKO with color type,default red
#' @param species                gene go annotation for species,must in list for annorad annotation list
#' @param indir                  database indir
#' @param ko2map                 ko2map file for ko and mapID
#' @param term2name_pre          term2name for kegg
#' @param sep                    marker file sep
#' @param category               category for species[animal or plant... ]
#' @param pvalueCutoff           compareCluster analysis for param,pvalueCutoff
#' @param pAdjustMethod          compareCluster analysis for param,pAdjustMethod
#' @param qvalueCutoff           compareCluster analysis for param,qvalueCutoff
#'
#' @return                       S4 object for result, from clusterprofiler
#' @export
#' @importFrom clusterProfiler compareCluster
#' @importFrom clusterProfiler enricher
#' @importFrom data.table fread
#' @importFrom reshape2 melt
#' @examples   x<-kegg_enrich_local_compareCluster(species,gene_list[[1]],outdir)
kegg_enrich_local_compareCluster<-function(species,gene_list,outdir,ko_color='red',ensemble_symbol_file=NULL,type='KEGG_annotate',indir="/annoroad/bioinfo/PMO/database/config_database",category='animal',term2name_pre='/annoroad/data1/bioinfo/PMO/Public/database/Public/Quarter/KEGG/current/data/pathway_',ko2map='/annoroad/data1/bioinfo/PMO/Public/database/Public/Quarter/KEGG/current/data/ko2map/ko2map.xls',pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff =0.05){
  ko2map<-go_Pathways(ko2map)
  ko_gile<-go_listfile(species,indir=indir,type=type)  #KEGG_annotate type="GO_annotate"
  gene_ko<-unique(gene_ko_update(outdir,ko_gile))
  gene_ko<-gene_ko[gene_ko$KO!='.',]
  ko2map=unique(reshape2::melt(ko2map))
  colnames(ko2map)<-c('kegg_term','KO')
  kegg_ko2gene<-unique(na.omit(merge(ko2map,gene_ko,by="KO",all=TRUE)))
  print(gene_list[[1]])
  gene_list2<-list()
  if (is.null(ensemble_symbol_file)){
    print("不进行转换基因名称转换")
    gene_list2<-gene_list
  }else{
    kegg2gene_symbol<-ensemble2symbol(ensemble_symbol_file,kegg_ko2gene,type='kegg') ##c("GO","Gene_symbol")
    colnames(kegg2gene_symbol)[3]<-'Gene'
    genes<-read.table(ensemble_symbol_file,sep='\t')
    colnames(genes)<-c('Gene','Gene_symbol')
    gene_list2<-lapply(gene_list,id2symbol,genes=genes)#function(x){as.vector(genes[genes$Gene %in%x,]$Gene_symbol)})
    kegg_ko2gene<-kegg2gene_symbol
    #lapply(gene_list,function(x){})
  }
  kegg_term2name <- fread(paste(term2name_pre,category,'.list',sep=''),header=F,sep="\t",stringsAsFactors = FALSE,quote = "",data.table = F)
  kegg_term2name$V1<-gsub('path:','',kegg_term2name$V1)
  colnames(kegg_term2name)<-c('kegg_term','name')
  kegg_ko2gene_tmp<-kegg_ko2gene[,2:3]
  #print(head(kegg_ko2gene_tmp))
  #print(head(gene_list[[1]]))
  #print(head(gene_list2[[1]]))
  if (length(names(gene_list))>1){
    x <- compareCluster(gene_list2, fun='enricher',TERM2GENE=kegg_ko2gene_tmp,TERM2NAME=kegg_term2name,pvalueCutoff = pvalueCutoff,  pAdjustMethod = pAdjustMethod ,qvalueCutoff =qvalueCutoff)
    x@compareClusterResult$geneKO<-as.vector(unlist(lapply(x@compareClusterResult$geneID, function(x){ paste(as.vector(kegg_ko2gene[kegg_ko2gene$Gene %in% unlist(strsplit(x,split = "/",fixed=T)),]$KO),paste(' ',ko_color,sep=''),sep='',collapse='|')} )))
  }else{
    x <- enricher(gene=gene_list2,pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff =qvalueCutoff,TERM2GENE = kegg_ko2gene_tmp,TERM2NAME =kegg_term2name )
    x@result$geneKO<-as.vector(unlist(lapply(x@result$geneID, function(x){ paste(as.vector(kegg_ko2gene[kegg_ko2gene$Gene %in% unlist(strsplit(x,split = "/",fixed=T)),]$KO),paste(' ',ko_color,sep=''),sep='',collapse='|')} )))
  }
  return (x)
}

#' Title kegg_enrich_local_compareCluster
#'
#' @param gene_list              a gene_list for enrichment, the gene list of length range from 1-max
#' @param outdir                 output dir for output file
#' @param ensemble_symbol_file   ensemble and symbol file from cellrange output
#' @param go.class               go.class file for go annotation
#' @param species                gene go annotation for species,must in list for annorad annotation list
#' @param indir                  database indir
#' @param ont                    Ontology for go ananlysis[BP CC MF]
#' @param pvalueCutoff           compareCluster analysis for param,pvalueCutoff
#' @param pAdjustMethod          compareCluster analysis for param,pAdjustMethod
#' @param qvalueCutoff           compareCluster analysis for param,qvalueCutoff
#'
#' @return                       S4 object for result, from clusterprofiler
#' @export
#' @importFrom clusterProfiler compareCluster
#' @importFrom clusterProfiler enricher
#' @importFrom data.table fread
#' @importFrom reshape2 melt
#' @examples   x<-kegg_enrich_local_compareCluster(species,gene_list[[1]],outdir)


go_enrich_local_compareCluster<-function(species,gene_list,outdir,indir="/annoroad/bioinfo/PMO/database/config_database",go.class='/annoroad/bioinfo/PMO/database/GO/current/go.class',ensemble_symbol_file=NULL,ont='all',pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff =0.05){
  go_gile<-go_listfile(species,indir=indir)
  go_gile2=go_updata(outdir,go_gile)
  ##gene10272 有问题，这是在建库的时候有问题
  print(paste("1 go_term2gene 数据准备，时间为：",Sys.time(),sep=''))
  go_list<-go_Pathways(go_gile2)
  go_term2gene=unique(reshape2::melt(go_list))
  colnames(go_term2gene)<-c("GO","Gene")
  print(paste("1 go_term2gene 数据完成准备，时间为：",Sys.time(),sep=''))
  gene_list2<-list()
  if (is.null(ensemble_symbol_file)){
    print("不进行转换基因名称转换")
    gene_list2<-gene_list
  }else{
    go_term2gene_symbol<-ensemble2symbol(ensemble_symbol_file,go_term2gene) ##c("GO","Gene_symbol")
    colnames(go_term2gene_symbol)<-c("GO","Gene")
    genes<-read.table(ensemble_symbol_file,sep='\t')
    colnames(genes)<-c('Gene','Gene_symbol')
    gene_list2<-lapply(gene_list,id2symbol,genes=genes)#function(x){as.vector(genes[genes$Gene %in%x,]$Gene_symbol)})
    go_term2gene<-go_term2gene_symbol
  }
  ###
  ###  term2name  数据准备
  #go.class<-'/annoroad/bioinfo/PMO/database/GO/current/go.class'
  term2name <- fread(go.class,header=T,sep="\t",stringsAsFactors = FALSE,quote = "",data.table = F)[,c('Ontology',"Accession","Term_name")]
  colnames(term2name)<-c('Ontology','go_term','Term_name')
  rownames(term2name)<-term2name$go_term
  if (ont=='all'){
    print("enrich all ont")
  }else if( ont=='CC'){
    term2name<-term2name[term2name$Ontology=='cellular_component',]
  }else if( ont=='BP'){
    term2name<-term2name[term2name$Ontology=='biological_process',]
  }else if( ont=='MF'){
    term2name<-term2name[term2name$Ontology=='molecular_function',]
  }else{
    print(paste("输入的ont有问题，不在允许范围内，其可以输入为：CC/BP/MF,请注意核对参数，小伙你输入的ONT为：",ont,sep=''))
  }
  #gene_list<-list(cluster1=sample(go_term2gene$Gene,20),cluster2=sample(go_term2gene$Gene,200),cluster3=sample(go_term2gene$Gene_name,100))
  #go_enrich <- enricher(gene=gene_list[[1]],pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE = go_term2gene,TERM2NAME = term2name[,c('go_term','name')])
  #x <- compareCluster(gene_list, fun='enricher',TERM2GENE=go_term2gene,TERM2NAME=term2name[,c('go_term','Description')],pvalueCutoff =0.05, pAdjustMethod ='BH')##, qvalueCutoff = 0.1
  x<-list()
  a<-term2name[,c(2,3)]
  #print(head(gene_list2[1]))
  if (length(names(gene_list))>1){
    x <- compareCluster(gene_list2, fun='enricher',TERM2GENE=go_term2gene,TERM2NAME=a,pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff =qvalueCutoff)
    #print(head(x@compareClusterResult))
    x@compareClusterResult<-as.data.frame(x@compareClusterResult)
    #x@compareClusterResult$Description<-term2name[as.vector(x@compareClusterResult$ID),]$name
    #print(head(x@compareClusterResult)) #Description
    x@compareClusterResult$Ontology<-term2name[as.vector(x@compareClusterResult$ID),]$Ontology
    x@compareClusterResult<-na.omit(x@compareClusterResult)
  }else{
    x <- enricher(gene=gene_list2,pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff =qvalueCutoff,TERM2GENE = go_term2gene,TERM2NAME = a)
    #print(head(x@result))
    x@result<-as.data.frame(x@result)
    #x@result$Description<-term2name[as.vector(x@result$ID),]$name
    x@result$Ontology<-term2name[as.vector(x@result$ID),]$Ontology
    x@result<-na.omit(x@result)
  }
  return (x)
}





#' Title get gene symbol from ensemble name
#'
#' @param genefile genes.tsv from cellranger output
#' @param term2gene  data.fram, the first col is Go term names and the second col is Gene name
#' @param sep genefile sep, default \t
#' @return term2gene  data.fram, the first col is Go term names and the second col is Gene name(symbol)
#'
#' @examples go_term2gene_symbol<-ensemble2symbol(genefile,term2gene,sep="\t",type='go')
ensemble2symbol<-function(genefile,term2gene,sep="\t",type='go'){
  print(paste("开始读取gens.tsv文件，时间为：",Sys.time(),sep=" "))
  genes<-read.table(genefile,sep=sep)
  colnames(genes)<-c("Gene","Gene_symbol")
  print(paste("开始读取匹配ensemble和symbol文件，时间为：",Sys.time(),sep=" "))
  term2gene2<-merge(term2gene,genes,by="Gene",all=TRUE)
  print(paste("删掉go或者symbol为na的行，时间为：",Sys.time(),sep=" "))
  if (type=='go'){
    term2gene2<-term2gene2[!is.na(term2gene2$GO)&!is.na(term2gene2$Gene_symbol),c('GO','Gene_symbol')]
  }else{
    term2gene2<-term2gene2[!is.na(term2gene2$kegg_term)&!is.na(term2gene2$Gene_symbol),c('KO','kegg_term','Gene_symbol')]
  }
  print(paste("完成ensemble和symbol的转换，时间为：",Sys.time(),sep=" "))
  return (term2gene2)
}#go_term2gene_symbol<-ensemble2symbol(ensemble_symbol_file,go_term2gene)

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
#' @importFrom clusterProfiler compareCluster
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
#' @importFrom clusterProfiler dotplot emapplot
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
#' @export
#' @importFrom clusterProfiler compareCluster dotplot emapplot
#' @importFrom dplyr select mutate
#' @importFrom vroom vroom
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


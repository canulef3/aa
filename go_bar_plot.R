###################################
############GO ANALYSIS############
###################################
library("biomaRt")
library("org.Hs.eg.db")
library("GO.db")
library("GOstats")
library("goProfiles")
library("topGO")
library("clusterProfiler")
library("Rgraphviz")
library("ggplot2")
library("tidyverse")
setwd("/Users/francisca/Desktop/")

#1. GENERAR UNA LISTA EN ENTREZ DE LOS GENES QUE QUIERO ESTUDIAR
up_genes <- read.table("~/Desktop/DATA_TESIS/trabajos/common_gene_list_up_regulated_edge_deseq.txt", quote="\"", comment.char="")
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
human_list <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = up_genes, mart = ensembl )
write.table(human_list, file = "my_genes_up.txt", quote = FALSE, sep = "\t", append = FALSE)
my_genes_up <- read.delim("~/Desktop/DATA_TESIS/trabajos/my_genes_up.txt")
my_genes_up[,1] -> aa
as.vector(as.character(aa)) -> aa2
is.vector(aa2)
#2. Análisis GO
ego<-enrichGO(aa2,'org.Hs.eg.db', ont="BP",pvalueCutoff = 0.1 ,pAdjustMethod = "BH",readable = TRUE)
x2=enrichGO(aa2, OrgDb='org.Hs.eg.db', ont='ALL', pvalueCutoff = 0.1, pAdjustMethod = "BH",readable = TRUE)
x3<-x2@result #extraer la tabla de result a partir de los datos de X
x3 <- mutate(x3, -log10(x3$p.adjust)) #calcular el -log10(padj)
colnames(x3) <- c( "ONTOLOGY","ID","Description","GeneRatio", "BgRatio","pvalue","p.adjust","qvalue", "geneID","Count","EnrichmentScore")
#le cambié el nombre a la nueva columna por uno mas facil de describir


x4 <- x3 %>% 
  group_by(ONTOLOGY) %>% #agrupar los datos por ontologia genica
  arrange() %>% #ordenar los datos en orden creciente
  top_n(n=10) #seleccionar los primeros n de cada grupo

x4$Description <- factor(x4$Description, levels = x4$Description[order(x4$EnrichmentScore, decreasing = TRUE)])

"Gene Ontology" -> legend_title
g2 <- ggplot(x4, aes(Description, EnrichmentScore, fill=ONTOLOGY, split=ONTOLOGY)) 
g2 + geom_bar(stat = "identity") +
  scale_fill_manual(legend_title, labels = c("Biological Process", "Cellular Component", "Molecular Function"), values=c("#E495A5","#ABB065","#39BEB1"))+
  xlab("GO category") +
  ylab("Enrichment score of -Log10(pvalue)") +
  labs(title = "GO Enrichment in hnpcs",
       subtitle = "All Ontologies") +
  theme(axis.text.x=element_text(angle=-40, hjust=0, vjust=0.6),
        legend.text = element_text(size = 14),
        legend.title = element_text(face="bold", size = 16)) + 
  facet_wrap(.~ONTOLOGY, scales = "free_x") + #cambie el grid por wrap
  theme(axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10),
        legend.key.size =  unit(0.1, "in"), plot.title = element_text(size = 18))

setwd("C:/Users/stephan123/Desktop/DESeq2")

# load packages

library(DESeq2)

library(data.table)

library(readr)

library(biomaRt)

library(RColorBrewer)

library(pheatmap)

library(ggplot2)

library(devtools)

library(enrichR)

library(ggrepel)

# Set directory

dir<- "C:/Users/stephan123/Desktop/DESeq2"


#Construct count matrix

HTseqFiles<- grep('ARM',list.files(dir), value = TRUE)

HTseqAGE= c("48","54","48","52","49","55","58","46","64","29","33","31","34")

HTseqGender= c("F","F","F","F","F","F","F","F","M","F","M","M","M")

HTseqConditions<- c(rep("diseased",9),rep("healthy",4))

HTseqNAMES <- c("P1","P2","P3","P4","P5","P6","P7","P8","P9","C1","C2","C6","C7")

HTseqTable<- data.frame(sampleNames=HTseqNAMES,fileName=HTseqFiles,condition=HTseqConditions,age=HTseqAGE,gender=HTseqGender)

View(HTseqTable)

# subset data set to remove one male in cases and one female in controls/ to not skew analysis.
HTseqTable <- subset(HTseqTable, ! HTseqTable$fileName %in% c("CASE_09_ARM.txt","CTRL_01_ARM.txt","CASE_01_ARM.txt"))
View(HTseqTable)

#create DESeq dataset.
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = HTseqTable,directory = dir,design = ~ condition)

#biomart
listMarts()
genes <- rownames(ddsHTSeq)
list(genes)

ensembl=useMart("ensembl")

ensembl = useDataset("hsapiens_gene_ensembl",mart = ensembl)

#set filters for biotype.
geneBiotype <- getBM(attributes = c('ensembl_gene_id','gene_biotype','description','hgnc_symbol','chromosome_name'),filters = 'ensembl_gene_id',values = genes,mart = ensembl,verbose = FALSE)
geneBiotype$description <- lapply(geneBiotype$description, function(x) gsub(" \\[.*","",x))
head(geneBiotype)

geneBiotype_sub <- subset(geneBiotype, geneBiotype$gene_biotype == 'lincRNA')#c("miRNA","scRNA","sRNA","snRNA","lincRNA","macro_lncRNA","misc_RNA","scaRNA","snoRNA"))
head(geneBiotype_sub)
unique(geneBiotype_sub$gene_biotype)

ddsHTSeq <- subset(ddsHTSeq, rownames(ddsHTSeq) %in% geneBiotype_sub$ensembl_gene_id)

#simple filter counts below 1
ddsHTSeq<- ddsHTSeq[rowSums(counts(ddsHTSeq))>1,]

#set reference for comparison
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref='healthy') 

#standard analysis
ddsHTSeq_expres <- DESeq(ddsHTSeq)

res<- results(ddsHTSeq_expres, alpha = 0.05)
res
summary(res)

# stats analysis RLog
Rlog_HTSeq<- rlog(ddsHTSeq_expres, blind = FALSE)
head(assay(Rlog_HTSeq))


# distance estimates/ heatmap. shows similarity between samples.
HTSeq_Distance <- dist(t(assay(Rlog_HTSeq)))
HTSeq_distmatrix <- as.matrix(HTSeq_Distance) 
rownames(HTSeq_distmatrix) <- paste(Rlog_HTSeq$sampleNames,Rlog_HTSeq$gender,  Rlog_HTSeq$age, Rlog_HTSeq$condition, sep = " | ") 
colnames(HTSeq_distmatrix) <- NULL
colours <- colorRampPalette(rev(brewer.pal(9,"Blues")) )(255)

plot(HTSeq_Distance)
# save graph as PDF in file DESeq2.
pdf('Dist.Hest.pdf', width = 10, height = 7, onefile = FALSE)


HTseq_heatmap <- pheatmap(HTSeq_distmatrix,clustering_distance_rows = HTSeq_Distance,clustering_distance_cols = HTSeq_Distance,
                          col=colours,border_color = FALSE,main = "",fontface= "bold",fontsize = 18,fontsize_col = 18)
HTseq_heatmap
dev.off()



# DRAWING PCA plots.
HTSeq_PCA <- plotPCA(Rlog_HTSeq, intgroup= c("condition","gender","age"), returnData=TRUE)
percentVar <- round(100 * attr(HTSeq_PCA, "percentVar"))

# PCA plots by ggplot2
pdf('PCA_diseased_2.pdf',width = 10,height = 7, onefile = FALSE)

gg<- ggplot(HTSeq_PCA, aes(PC1,PC2,color=condition, shape=gender))+
  geom_point(size=3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  coord_fixed() +
  geom_text_repel(aes(label=paste(HTSeq_PCA$name, sep = '-')),size=4)+
  scale_x_continuous(limits = c(-25,+60))+
  scale_y_continuous(limits = c(-35,+20))+
  theme_classic()+
  ggtitle(label = 'PCA plot of gender vs condition')+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
gg

dev.off()

#heatmap RES expression profiles.

ENR_bio <- subset(geneBiotype, geneBiotype$ensembl_gene_id %in% selected, select = c('ensembl_gene_id','hgnc_symbol'))

# heatmap count matrix.
res <- res[order(res$padj),]
res <- subset(res, res$padj < 0.05)
selected <- rownames(res)

sig <- subset(Rlog_HTSeq, rownames(Rlog_HTSeq) %in% selected)

pdf('heatmap_RES.pdf',width = 10,height = 17, onefile = FALSE)

log2.norm.count <- assay(Rlog_HTSeq)[selected,]
log2.norm.count <- log2.norm.count - rowMeans(log2.norm.count)
rownames(log2.norm.count) <- paste(rownames(log2.norm.count),
                                   geneBiotype$description[match(rownames(log2.norm.count),
                                                                 geneBiotype$ensembl_gene_id)], sep = " - ")
df <- as.data.frame(colData(Rlog_HTSeq)[,c("condition", "gender")])
heat<-pheatmap(log2.norm.count,cluster_rows = TRUE, clustering_distance_rows = 'euclidean', 
         show_rownames = TRUE,cluster_cols = FALSE, cellwidth = 10, 
         cellheight = 10, annotation_col = df, fontsize_row = 8,fontsize_col = 8,
         main = "Heatmap of expressed genes",fontface= "bold",fontsize = 13,border_color = FALSE)
heat
dev.off()



#enriched genes

ENR_bio <- subset(geneBiotype, geneBiotype$ensembl_gene_id %in% selected, select = c('ensembl_gene_id','hgnc_symbol'))

for_enrichr <- subset(ENR_bio$hgnc_symbol,  ENR_bio$hgnc_symbol != "")

dbs <- listEnrichrDbs()$libraryName

enriched_genes <- enrichr(for_enrichr,dbs)

rich <- subset(enriched_genes, lapply(enriched_genes, function(x) nrow(x)) > 0)

View(enriched_genes)

# volcanoplot

new <- as.data.frame(subset(res, select = c(log2FoldChange, pvalue)))
# merge data
new$HGNC <- ENR_bio$hgnc_symbol[match(rownames(new), ENR_bio$ensembl_gene_id)]

pdf('Volcano_Log2.pdf',width = 10,height = 7 ,onefile = FALSE)
vovolcanoplot <- ggplot(data = new, aes(x = log2FoldChange ,y = -log10(pvalue)))+ 
  geom_point(data = new,size=2, color = ifelse(new$HGNC %in% lb, "red", "black"))+
  theme_gray()+
  geom_label_repel(inherit.aes = TRUE,fill='blue',segment.color = 'Blue',label = ifelse(new$HGNC %in% lb,new$HGNC, ''),color = ifelse(new$HGNC %in% lb, "white", "black")) +
  ggtitle("Gene expression change Log2")+
  theme(title = element_text(size = 12,face = "bold",hjust = 0.5),
        axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10))+
  ylab("-log10 (pvalue)")
vovolcanoplot
dev.off()

lb <- c("TTTY14","TTTY15","XIST","LINC00265","LINC00310","TSIX","MIR155HG","LINC00302","C20orf197","FAM30A")


# gene enrichr
they_affect <- lapply(rich$TRANSFAC_and_JASPAR_PWMs$Term, function(x) gsub(" \\(.*", "", x))
c <- enrichr(they_affect, dbs)

they_affect_2 <- lapply(rich$Human_Gene_Atlas$Term, function(x) gsub("_.*","", x))
c2 <- enrichr(they_affect_2, dbs)

they_affect_3 <- lapply(rich$ENCODE_Histone_Modifications_2015$Term, function(x) gsub("_.*","", x))
c3 <- enrichr(they_affect_3, dbs)

they_affect_5 <- lapply(rich$LINCS_L1000_Ligand_Perturbations_up$Term, function(x) gsub("-.*","",x))
c5 <- enrichr(they_affect_5, dbs)

# Results tabel

HGNC_Symbol <- c(new$HGNC)


results_table <- new
results_table$Description <- geneBiotype$description[match(row.names(results_table), geneBiotype$ensembl_gene_id)]
 #
Result_tab <- as.data.frame(res[,c('log2FoldChange','pvalue','padj')])
Result_tab$description <- geneBiotype_sub$description[match(rownames(Result_tab),geneBiotype_sub$ensembl_gene_id)]
Result_tab$hgnc <- geneBiotype_sub$hgnc_symbol[match(rownames(Result_tab),geneBiotype_sub$ensembl_gene_id)]
Result_tab$chromosome_name <- geneBiotype_sub$chromosome_name[match(rownames(Result_tab),geneBiotype_sub$ensembl_gene_id)]
Results <- subset(Result_tab,rownames(Result_tab) %in% c("ENSG00000176659","ENSG00000188185","ENSG00000234883","ENSG00000226777","ENSG00000176075","ENSG00000227456","ENSG00000270641","ENSG00000229807","ENSG00000176728","ENSG00000233864"))
head(Results)
# THE END.

# CITATIONS
citation("DESeq2")

citation("biomaRt") 

citation("pheatmap")

citation("ggplot2")

citation("enrichR")



test <- function(z) {
  lapply(rich, function(x) subset(x, x$Genes == z, select=c('Term','Genes')))
  }

test('LINC01094')

they_affect <- lapply(rich$TRANSFAC_and_JASPAR_PWMs$Term, function(x) gsub(" \\(.*", "", x))

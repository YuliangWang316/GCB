library(future)
plan("multiprocess", workers = 20)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#seurat 4.1.0
options(warn=-1)
set.seed(1)
pbmc1.data <- Read10X(data.dir="D:/BOAO single-cell sequencing/GCB1/outs/filtered_feature_bc_matrix/")
pbmc2.data <- Read10X(data.dir="D:/BOAO single-cell sequencing/GCB2/outs/filtered_feature_bc_matrix/")
pbmc1.data <- as.data.frame(pbmc1.data)
pbmc2.data <- as.data.frame(pbmc2.data)
for (i in 1:length(colnames(pbmc1.data))) {
  colnames(pbmc1.data)[i] <- paste(colnames(pbmc1.data)[i],"pbmc1",i,sep = "-")  
}

for (i in 1:length(colnames(pbmc2.data))) {
  colnames(pbmc2.data)[i] <- paste(colnames(pbmc2.data)[i],"pbmc2",i,sep = "-")  
}

pbmc.data <-cbind(pbmc1.data,pbmc2.data)
pbmc.data <- as.data.frame(pbmc.data)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,Ighd <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB_filtered_noIgd", min.cells = 3, min.features = 200)
rm(pbmc.data,pbmc.data_t,pbmc1.data,pbmc2.data)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
counts<-quantile(pbmc@meta.data$nCount_RNA,c(0.025,0.975))
Features<-quantile(pbmc@meta.data$nFeature_RNA,c(0.025,0.975))
percentmt<-quantile(pbmc@meta.data$percent.mt,c(0.025,0.975))
pbmc <- subset(pbmc, subset = nFeature_RNA > 754 & nFeature_RNA < 5169 & percent.mt < 41.25 & nCount_RNA < 32574 & nCount_RNA >1656 & percent.mt > 4.13)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)

all.genes <- rownames(pbmc)
pbmc1 <- ScaleData(pbmc,features = all.genes)

remove(counts,Features,i,percentmt)

pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(object = pbmc1),npcs = 1000)

#pbmc1_1 <- RunPCA(pbmc1, features = all.genes,npcs = 1000)

S1<-data.frame(0,0,0,0)
colnames(S1)<-c("j","p","q","r")

S1_1<-data.frame(0,0,0,0)
colnames(S1_1)<-c("j","p","q","r")




for (j in 5:1000) {
  for (p in c(50,100,200,300,400,500,600,700,800,900,1000)) {
    for (q in c(20,30,40,50,60,70,80,90,100)) {
      for (r in seq(0.5,10,by=0.1)) {
        
        pbmc1 <- FindNeighbors(pbmc1, dims = 1:j,n.trees = p,k.param = q)
        
        pbmc1_1 <- FindNeighbors(pbmc1_1, dims = 1:j,n.trees = p,k.param = q)
        
        pbmc1 <- FindClusters(pbmc1, resolution = r)
        
        pbmc1_1 <- FindClusters(pbmc1_1, resolution = r)
        
        
        pbmc1.markers <- FindAllMarkers(pbmc1, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        
        pbmc1_1.markers <- FindAllMarkers(pbmc1_1, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        if ("Myc" %in% pbmc1.markers$gene) {
          Myc1<-filter(pbmc1.markers,gene == "Myc")
          if (length(rownames(Myc1)) <=2){
            a<-unique(Myc1$cluster)
            Myc1_seurat<-subset(pbmc1,idents = a)
            Myc1_data<-FetchData(Myc1_seurat,vars = "Myc")
            b<-mean(Myc1_data$Myc)
            Myc1_data_hi<-as.data.frame(Myc1_data[Myc1_data > b,])
            c<-length(rownames(Myc1_data_hi))/length(rownames(Myc1_data))
            if (c >= 0.5){
              P1<-data.frame(j,p,q,r)
              colnames(P1)<-c("j","p","q","r")
              S1<-rbind(S1,P1)
            }
            
            
            
          }
        }
        
        
        if ("Myc" %in% pbmc1_1.markers$gene) {
          Myc1_1<-filter(pbmc1_1.markers,gene == "Myc")
          if (length(rownames(Myc1_1)) <=2){
            a<-unique(Myc1_1$cluster)
            Myc1_1_seurat<-subset(pbmc1_1,idents = a)
            Myc1_1_data<-FetchData(Myc1_1_seurat,vars = "Myc")
            b<-mean(Myc1_1_data$Myc)
            Myc1_1_data_hi<-as.data.frame(Myc1_1_data[Myc1_1_data > b,])
            c<-length(rownames(Myc1_1_data_hi))/length(rownames(Myc1_1_data))
            if (c >= 0.5){
              P1_1<-data.frame(i,j,p,q,r)
              colnames(P1_1)<-c("j","p","q","r")
              S1_1<-rbind(S1_1,P1_1)
            }
            
            
            
          }
        }
        
        
      }
    }
  }
}

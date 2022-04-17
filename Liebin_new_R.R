library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#seurat 4.1.0
options(warn=-1)
set.seed(1)
pbmc.data <- read.table("D:/Liebing_DZLZRT_Rawdata/LieBingData_counts and normalize count/RTLZ.counts.tsv",sep = "\t",header = TRUE,row.names = 1)

pbmc.data <- as.data.frame(pbmc.data)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,Ighd <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB_filtered_noIgd", min.cells = 3, min.features = 200)
rm(pbmc.data,pbmc.data_t)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
counts<-quantile(pbmc@meta.data$nCount_RNA,c(0.025,0.975))
Features<-quantile(pbmc@meta.data$nFeature_RNA,c(0.025,0.975))
percentmt<-quantile(pbmc@meta.data$percent.mt,c(0.025,0.975))
pbmc <- subset(pbmc, subset = nFeature_RNA > 287 & nFeature_RNA < 2286 & percent.mt < 14.5 & nCount_RNA < 5432 & nCount_RNA >385 & percent.mt > 1.85)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)

all.genes <- rownames(pbmc)
pbmc1 <- ScaleData(pbmc,features = all.genes)
pbmc2 <- ScaleData(pbmc,features = all.genes,model.use = 'poison')
# samle as pbmc2 pbmc3 <- ScaleData(pbmc,features = all.genes,model.use = 'negbinom')
pbmc3 <- ScaleData(pbmc,features = all.genes,vars.to.regress = "percent.mt")
pbmc4 <- ScaleData(pbmc,features = all.genes,vars.to.regress = "nFeature_RNA")
pbmc5 <- ScaleData(pbmc,features = all.genes,vars.to.regress = "nCount_RNA")
pbmc6 <- ScaleData(pbmc,features = all.genes,vars.to.regress = c("percent.mt","nCount_RNA","nFeature_RNA"))
remove(counts,Features,i,percentmt)

pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(object = pbmc1),npcs = 1000)
pbmc2 <- RunPCA(pbmc2, features = VariableFeatures(object = pbmc2),npcs = 1000)
pbmc3 <- RunPCA(pbmc3, features = VariableFeatures(object = pbmc3),npcs = 1000)  
pbmc4 <- RunPCA(pbmc4, features = VariableFeatures(object = pbmc4),npcs = 1000)
pbmc5 <- RunPCA(pbmc5, features = VariableFeatures(object = pbmc1),npcs = 1000)
pbmc6 <- RunPCA(pbmc6, features = VariableFeatures(object = pbmc6),npcs = 1000)
pbmc1_1 <- RunPCA(pbmc1, features = all.genes,npcs = 1000)
pbmc2_1 <- RunPCA(pbmc2, features = all.genes,npcs = 1000)
pbmc3_1 <- RunPCA(pbmc3, features = all.genes,npcs = 1000)  
pbmc4_1 <- RunPCA(pbmc4, features = all.genes,npcs = 1000)
pbmc5_1 <- RunPCA(pbmc5, features = all.genes,npcs = 1000)
pbmc6_1 <- RunPCA(pbmc6, features = all.genes,npcs = 1000)
S1<-data.frame(0,0,0,0)
colnames(S1)<-c("j","p","q","r")
S2<-data.frame(0,0,0,0)
colnames(S2)<-c("j","p","q","r")
S3<-data.frame(0,0,0,0)
colnames(S3)<-c("j","p","q","r")
S4<-data.frame(0,0,0,0)
colnames(S4)<-c("j","p","q","r")
S5<-data.frame(0,0,0,0)
colnames(S5)<-c("j","p","q","r")
S6<-data.frame(0,0,0,0)
colnames(S6)<-c("j","p","q","r")
S1_1<-data.frame(0,0,0,0)
colnames(S1_1)<-c("j","p","q","r")
S2_1<-data.frame(0,0,0,0)
colnames(S2_1)<-c("j","p","q","r")
S3_1<-data.frame(0,0,0,0)
colnames(S3_1)<-c("j","p","q","r")
S4_1<-data.frame(0,0,0,0)
colnames(S4_1)<-c("j","p","q","r")
S5_1<-data.frame(0,0,0,0)
colnames(S5_1)<-c("j","p","q","r")
S6_1<-data.frame(0,0,0,0)
colnames(S6_1)<-c("j","p","q","r")



for (j in 5:i) {
  for (p in c(50,100,200,300,400,500,600,700,800,900,1000)) {
    for (q in c(20,30,40,50,60,70,80,90,100)) {
      for (r in seq(0.8,10,by=0.1)) {
        
        pbmc1 <- FindNeighbors(pbmc1, dims = 1:j,n.trees = p,k.param = q)
        pbmc2 <- FindNeighbors(pbmc2, dims = 1:j,n.trees = p,k.param = q)
        pbmc3 <- FindNeighbors(pbmc3, dims = 1:j,n.trees = p,k.param = q)
        pbmc4 <- FindNeighbors(pbmc4, dims = 1:j,n.trees = p,k.param = q)
        pbmc5 <- FindNeighbors(pbmc5, dims = 1:j,n.trees = p,k.param = q)
        pbmc6 <- FindNeighbors(pbmc6, dims = 1:j,n.trees = p,k.param = q)
        pbmc1_1 <- FindNeighbors(pbmc1_1, dims = 1:j,n.trees = p,k.param = q)
        pbmc2_1 <- FindNeighbors(pbmc2_1, dims = 1:j,n.trees = p,k.param = q)
        pbmc3_1 <- FindNeighbors(pbmc3_1, dims = 1:j,n.trees = p,k.param = q)
        pbmc4_1 <- FindNeighbors(pbmc4_1, dims = 1:j,n.trees = p,k.param = q)
        pbmc5_1 <- FindNeighbors(pbmc5_1, dims = 1:j,n.trees = p,k.param = q)
        pbmc6_1 <- FindNeighbors(pbmc6_1, dims = 1:j,n.trees = p,k.param = q)
        pbmc1 <- FindClusters(pbmc1, resolution = r)
        pbmc2 <- FindClusters(pbmc2, resolution = r)
        pbmc3 <- FindClusters(pbmc3, resolution = r)
        pbmc4 <- FindClusters(pbmc4, resolution = r)
        pbmc5 <- FindClusters(pbmc5, resolution = r)
        pbmc6 <- FindClusters(pbmc6, resolution = r)
        pbmc1_1 <- FindClusters(pbmc1_1, resolution = r)
        pbmc2_1 <- FindClusters(pbmc2_1, resolution = r)
        pbmc3_1 <- FindClusters(pbmc3_1, resolution = r)
        pbmc4_1 <- FindClusters(pbmc4_1, resolution = r)
        pbmc5_1 <- FindClusters(pbmc5_1, resolution = r)
        pbmc6_1 <- FindClusters(pbmc6_1, resolution = r)
        
        pbmc1.markers <- FindAllMarkers(pbmc1, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        pbmc2.markers <- FindAllMarkers(pbmc2, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        pbmc3.markers <- FindAllMarkers(pbmc3, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        pbmc4.markers <- FindAllMarkers(pbmc4, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        pbmc5.markers <- FindAllMarkers(pbmc5, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        pbmc6.markers <- FindAllMarkers(pbmc6, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        pbmc1_1.markers <- FindAllMarkers(pbmc1_1, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        pbmc2_1.markers <- FindAllMarkers(pbmc2_1, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        pbmc3_1.markers <- FindAllMarkers(pbmc3_1, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        pbmc4_1.markers <- FindAllMarkers(pbmc4_1, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        pbmc5_1.markers <- FindAllMarkers(pbmc5_1, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        pbmc6_1.markers <- FindAllMarkers(pbmc6_1, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
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
        
        if ("Myc" %in% pbmc2.markers$gene){
          Myc2<-filter(pbmc2.markers,gene == "Myc")
          if (length(rownames(Myc2)) <=2){
            a<-unique(Myc2$cluster)
            Myc2_seurat<-subset(pbmc2,idents = a)
            Myc2_data<-FetchData(Myc2_seurat,vars = "Myc")
            b<-mean(Myc2_data$Myc)
            Myc2_data_hi<-as.data.frame(Myc2_data[Myc2_data > b,])
            c<-length(rownames(Myc2_data_hi))/length(rownames(Myc2_data))
            if (c >= 0.5){
              P2<-data.frame(j,p,q,r)
              colnames(P2)<-c("j","p","q","r")
              S2<-rbind(S2,P2)
            }
            
            
            
          }
        }
        
        if ("Myc" %in% pbmc3.markers$gene){
          Myc3<-filter(pbmc3.markers,gene == "Myc")
          if (length(rownames(Myc3)) <=2){
            a<-unique(Myc3$cluster)
            Myc3_seurat<-subset(pbmc3,idents = a)
            Myc3_data<-FetchData(Myc3_seurat,vars = "Myc")
            b<-mean(Myc3_data$Myc)
            Myc3_data_hi<-as.data.frame(Myc3_data[Myc3_data > b,])
            c<-length(rownames(Myc3_data_hi))/length(rownames(Myc3_data))
            if (c >= 0.5){
              P3<-data.frame(j,p,q,r)
              colnames(P3)<-c("j","p","q","r")
              S3<-rbind(S3,P3)
            }
            
            
            
          }
        }
        if ("Myc" %in% pbmc4.markers$gene) {
          Myc4<-filter(pbmc4.markers,gene == "Myc")
          if (length(rownames(Myc4)) <=2){
            a<-unique(Myc4$cluster)
            Myc4_seurat<-subset(pbmc4,idents = a)
            Myc4_data<-FetchData(Myc4_seurat,vars = "Myc")
            b<-mean(Myc4_data$Myc)
            Myc4_data_hi<-as.data.frame(Myc4_data[Myc4_data > b,])
            c<-length(rownames(Myc4_data_hi))/length(rownames(Myc4_data))
            if (c >= 0.5){
              P4<-data.frame(j,p,q,r)
              colnames(P4)<-c("j","p","q","r")
              S4<-rbind(S4,P4)
            }
            
            
            
          }
        }
        
        if ("Myc" %in% pbmc5.markers$gene) {
          Myc5<-filter(pbmc5.markers,gene == "Myc")
          if (length(rownames(Myc5)) <=2){
            a<-unique(Myc5$cluster)
            Myc5_seurat<-subset(pbmc5,idents = a)
            Myc5_data<-FetchData(Myc5_seurat,vars = "Myc")
            b<-mean(Myc5_data$Myc)
            Myc5_data_hi<-as.data.frame(Myc5_data[Myc5_data > b,])
            c<-length(rownames(Myc5_data_hi))/length(rownames(Myc5_data))
            if (c >= 0.5){
              P5<-data.frame(j,p,q,r)
              colnames(P5)<-c("j","p","q","r")
              S5<-rbind(S5,P5)
            }
            
            
            
          }
        }
        if ("Myc" %in% pbmc6.markers$gene) {
          Myc6<-filter(pbmc6.markers,gene == "Myc")
          if (length(rownames(Myc6)) <=2){
            a<-unique(Myc6$cluster)
            Myc6_seurat<-subset(pbmc6,idents = a)
            Myc6_data<-FetchData(Myc6_seurat,vars = "Myc")
            b<-mean(Myc6_data$Myc)
            Myc6_data_hi<-as.data.frame(Myc6_data[Myc6_data > b,])
            c<-length(rownames(Myc6_data_hi))/length(rownames(Myc6_data))
            if (c >= 0.5){
              P6<-data.frame(j,p,q,r)
              colnames(P6)<-c("j","p","q","r")
              S6<-rbind(S6,P6)
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
        if ("Myc" %in% pbmc2_1.markers$gene) {
          Myc2_1<-filter(pbmc2_1.markers,gene == "Myc")
          if (length(rownames(Myc2_1)) <=2){
            a<-unique(Myc2_1$cluster)
            Myc2_1_seurat<-subset(pbmc2_1,idents = a)
            Myc2_1_data<-FetchData(Myc2_1_seurat,vars = "Myc")
            b<-mean(Myc2_1_data$Myc)
            Myc2_1_data_hi<-as.data.frame(Myc2_1_data[Myc2_1_data > b,])
            c<-length(rownames(Myc2_1_data_hi))/length(rownames(Myc2_1_data))
            if (c >= 0.5){
              P2_1<-data.frame(j,p,q,r)
              colnames(P2_1)<-c("j","p","q","r")
              S2_1<-rbind(S2_1,P2_1)
            }
            
            
            
          }
        }
        if ("Myc" %in% pbmc3_1.markers$gene) {
          Myc3_1<-filter(pbmc3_1.markers,gene == "Myc")
          if (length(rownames(Myc3_1)) <=2){
            a<-unique(Myc3_1$cluster)
            Myc3_1_seurat<-subset(pbmc3_1,idents = a)
            Myc3_1_data<-FetchData(Myc3_1_seurat,vars = "Myc")
            b<-mean(Myc3_1_data$Myc)
            Myc3_1_data_hi<-as.data.frame(Myc3_1_data[Myc3_1_data > b,])
            c<-length(rownames(Myc3_1_data_hi))/length(rownames(Myc3_1_data))
            if (c >= 0.5){
              P3_1<-data.frame(j,p,q,r)
              colnames(P3_1)<-c("j","p","q","r")
              S3_1<-rbind(S3_1,P3_1)
            }
            
            
            
          }
        }
        if ("Myc" %in% pbmc4_1.markers$gene) {
          Myc4_1<-filter(pbmc4_1.markers,gene == "Myc")
          if (length(rownames(Myc4_1)) <=2){
            a<-unique(Myc4_1$cluster)
            Myc4_1_seurat<-subset(pbmc4_1,idents = a)
            Myc4_1_data<-FetchData(Myc4_1_seurat,vars = "Myc")
            b<-mean(Myc4_1_data$Myc)
            Myc4_1_data_hi<-as.data.frame(Myc4_1_data[Myc4_1_data > b,])
            c<-length(rownames(Myc4_1_data_hi))/length(rownames(Myc4_1_data))
            if (c >= 0.5){
              P4_1<-data.frame(j,p,q,r)
              colnames(P4_1)<-c("j","p","q","r")
              S4_1<-rbind(S4_1,P4_1)
            }
            
            
            
          }
        }
        if ("Myc" %in% pbmc5_1.markers$gene) {
          Myc5_1<-filter(pbmc5_1.markers,gene == "Myc")
          if (length(rownames(Myc5_1)) <=2){
            a<-unique(Myc5_1$cluster)
            Myc5_1_seurat<-subset(pbmc5_1,idents = a)
            Myc5_1_data<-FetchData(Myc5_1_seurat,vars = "Myc")
            b<-mean(Myc5_1_data$Myc)
            Myc5_1_data_hi<-as.data.frame(Myc5_1_data[Myc5_1_data > b,])
            c<-length(rownames(Myc5_1_data_hi))/length(rownames(Myc5_1_data))
            if (c >= 0.5){
              P5_1<-data.frame(j,p,q,r)
              colnames(P5_1)<-c("j","p","q","r")
              S5_1<-rbind(S5_1,P5_1)
            }
            
            
            
          }
        }
        if ("Myc" %in% pbmc6_1.markers$gene) {
          Myc6_1<-filter(pbmc6_1.markers,gene == "Myc")
          if (length(rownames(Myc6_1)) <=2){
            a<-unique(Myc6_1$cluster)
            Myc6_1_seurat<-subset(pbmc6_1,idents = a)
            Myc6_1_data<-FetchData(Myc6_1_seurat,vars = "Myc")
            b<-mean(Myc6_1_data$Myc)
            Myc6_1_data_hi<-as.data.frame(Myc6_1_data[Myc6_1_data > b,])
            c<-length(rownames(Myc6_1_data_hi))/length(rownames(Myc6_1_data))
            if (c >= 0.5){
              P6_1<-data.frame(j,p,q,r)
              colnames(P6_1)<-c("j","p","q","r")
              S6_1<-rbind(S6_1,P6_1)
            }
            
            
            
          }
        }
        
      }
    }
  }
}

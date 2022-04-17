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

pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(object = pbmc1),npcs = 100)



S1<-data.frame(0,0,0,0)
colnames(S1)<-c("t","u","v","w")
for (t in 50:5) {
  for (u in c(50)) {
    for (v in c(20)) {
      for (w in seq(10,0.5,by=-0.1)) {
        
        pbmc1 <- FindNeighbors(pbmc1, dims = 1:t,n.trees = u,k.param = v)
        pbmc1 <- FindClusters(pbmc1, resolution = w)
        
        Myc<-FetchData(pbmc1,vars = "Myc")
        idents<-Idents(pbmc1)[colnames(pbmc1)]
        Myc$ident<-idents
        noise<-rnorm(n=length(x=Myc[,"Myc"]))/1e+05
        Myc[,"Myc"]<-Myc[,"Myc"]+noise
        Myc$ident<-factor(Myc$ident,levels = names(x=rev(x=sort(x=tapply(X=Myc[,"Myc"],INDEX = Myc$ident,FUN = mean),decreasing = FALSE))))
        
        Irf4<-FetchData(pbmc1,vars = "Irf4")
        idents<-Idents(pbmc1)[colnames(pbmc1)]
        Irf4$ident<-idents
        noise<-rnorm(n=length(x=Irf4[,"Irf4"]))/1e+05
        Irf4[,"Irf4"]<-Irf4[,"Irf4"]+noise
        Irf4$ident<-factor(Irf4$ident,levels = names(x=rev(x=sort(x=tapply(X=Irf4[,"Irf4"],INDEX = Irf4$ident,FUN = mean),decreasing = FALSE))))
        
        Kdm6b<-FetchData(pbmc1,vars = "Kdm6b")
        idents<-Idents(pbmc1)[colnames(pbmc1)]
        Kdm6b$ident<-idents
        noise<-rnorm(n=length(x=Kdm6b[,"Kdm6b"]))/1e+05
        Kdm6b[,"Kdm6b"]<-Kdm6b[,"Kdm6b"]+noise
        Kdm6b$ident<-factor(Kdm6b$ident,levels = names(x=rev(x=sort(x=tapply(X=Kdm6b[,"Kdm6b"],INDEX = Kdm6b$ident,FUN = mean),decreasing = FALSE))))
        
        Ly75<-FetchData(pbmc1,vars = "Ly75")
        idents<-Idents(pbmc1)[colnames(pbmc1)]
        Ly75$ident<-idents
        noise<-rnorm(n=length(x=Ly75[,"Ly75"]))/1e+05
        Ly75[,"Ly75"]<-Ly75[,"Ly75"]+noise
        Ly75$ident<-factor(Ly75$ident,levels = names(x=rev(x=sort(x=tapply(X=Ly75[,"Ly75"],INDEX = Ly75$ident,FUN = mean),decreasing = FALSE))))
        
      
        if (levels(Kdm6b$ident)[1] == levels(Irf4$ident)[1]){
          if (levels(Kdm6b$ident)[1] == levels(Ly75$ident)[1]){
            if(levels(Kdm6b$ident)[1] == levels(Myc$ident)[1]){
              
              
                  a<-levels(Myc$ident)[1]
                  Myc1_seurat<-subset(pbmc1,idents = a)
                  Myc1_data<-FetchData(Myc1_seurat,vars = "Myc")
                  b<-mean(Myc1_data$Myc)
                  Myc1_data_hi<-as.data.frame(Myc1_data[Myc1_data > b,])
                  c<-length(rownames(Myc1_data_hi))/length(rownames(Myc1_data))
                  
                  d<-levels(Myc$ident)[2]
                  Myc2_seurat<-subset(pbmc1,idents = d)
                  Myc2_data<-FetchData(Myc2_seurat,vars = "Myc")
                  e<-mean(Myc2_data$Myc)
                  Myc2_data_hi<-as.data.frame(Myc2_data[Myc2_data > e,])
                  f<-length(rownames(Myc2_data_hi))/length(rownames(Myc2_data))
                  
                  g<-levels(Irf4$ident)[1]
                  Irf4_seurat<-subset(pbmc1,idents = g)
                  Irf4_data<-FetchData(Irf4_seurat,vars = "Irf4")
                  h<-mean(Irf4_data$Irf4)
                  Irf4_data_hi<-as.data.frame(Irf4_data[Irf4_data > h,])
                  r<-length(rownames(Irf4_data_hi))/length(rownames(Irf4_data))
                  
                  s<-levels(Ly75$ident)[1]
                  Ly75_seurat<-subset(pbmc1,idents = s)
                  Ly75_data<-FetchData(Ly75_seurat,vars = "Ly75")
                  k<-mean(Ly75_data$Ly75)
                  Ly75_data_hi<-as.data.frame(Ly75_data[Ly75_data > k,])
                  l<-length(rownames(Ly75_data_hi))/length(rownames(Ly75_data))
                  
                  m<-levels(Kdm6b$ident)[1]
                  Kdm6b_seurat<-subset(pbmc1,idents = m)
                  Kdm6b_data<-FetchData(Kdm6b_seurat,vars = "Kdm6b")
                  n<-mean(Kdm6b_data$Kdm6b)
                  Kdm6b_data_hi<-as.data.frame(Kdm6b_data[Kdm6b_data > n,])
                  o<-length(rownames(Kdm6b_data_hi))/length(rownames(Kdm6b_data))
                  
                  if (c >= 0.5  ){
                    P1<-data.frame(t,u,v,w)
                    colnames(P1)<-c("t","u","v","w")
                    S1<-rbind(S1,P1)
                  
                  
                }
              
            
          }
        }
        
        
      }
    }
  }
  }
}
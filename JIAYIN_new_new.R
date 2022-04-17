library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#seurat 4.1.0
options(warn=-1)
set.seed(1)
pbmc.data <- Read10X(data.dir="D:/JYNJ1911CYH01_wxm_scRNAseq/GCB_30000/outs/filtered_feature_bc_matrix/")

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
pbmc <- subset(pbmc, subset = nFeature_RNA > 770 & nFeature_RNA < 2927 & percent.mt < 16.01 & nCount_RNA < 10419 & nCount_RNA >1560 & percent.mt > 1.48)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)

all.genes <- rownames(pbmc)
pbmc1 <- ScaleData(pbmc,features = all.genes)

remove(counts,Features,i,percentmt)

pbmc1 <- RunPCA(pbmc1, features = VariableFeatures(object = pbmc1),npcs = 1000)

mtor_1<-read.table("I:/MAC/Memory_B cells/PENG_LEUCINE_DEPRIVATION_DN_new.txt",sep = "\t",header = TRUE)
mtor_1_new<-mtor_1[2:length(rownames(mtor_1)),]
mtor_1_new<-tolower(mtor_1_new)
library(Hmisc)
mtor_1_new<-capitalize(mtor_1_new)
mtor1<-FetchData(pbmc1,vars = mtor_1_new)
for (i in 1:length(rownames(mtor1))) {
  mtor1$Average[i]<-mean(as.numeric(mtor1[i,]))
}

Total<-FetchData(pbmc1,vars = rownames(pbmc1))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_mtor1<-mean(mtor1$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_mtor1)
Sortdatalo<-filter(Sortdata,S < average_mtor1)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((mtor_1_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((mtor_1_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc1,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc1@meta.data$Mtor_1<-mtor1$Average-control$Average

mtor_2<-read.table("I:/MAC/Memory_B cells/PENG_RAPAMYCIN_RESPONSE_DN_new.txt",sep = "\t",header = TRUE)
mtor_2_new<-mtor_2[2:length(rownames(mtor_2)),]
mtor_2_new<-tolower(mtor_2_new)
library(Hmisc)
mtor_2_new<-capitalize(mtor_2_new)
mtor2<-FetchData(pbmc1,vars = mtor_2_new)
for (i in 1:length(rownames(mtor2))) {
  mtor2$Average[i]<-mean(as.numeric(mtor2[i,]))
}

Total<-FetchData(pbmc1,vars = rownames(pbmc1))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_mtor2<-mean(mtor2$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_mtor2)
Sortdatalo<-filter(Sortdata,S < average_mtor2)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((mtor_2_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((mtor_2_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc1,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc1@meta.data$Mtor_2<-mtor2$Average-control$Average

mtor_3<-read.table("I:/MAC/Memory_B cells/MTORC1_hallmarks.txt",sep = "\t",header = TRUE)
mtor_3_new<-mtor_3[2:length(rownames(mtor_3)),]
mtor_3_new<-tolower(mtor_3_new)
library(Hmisc)
mtor_3_new<-capitalize(mtor_3_new)
mtor3<-FetchData(pbmc1,vars = mtor_3_new)
for (i in 1:length(rownames(mtor3))) {
  mtor3$Average[i]<-mean(as.numeric(mtor3[i,]))
}

Total<-FetchData(pbmc1,vars = rownames(pbmc1))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_mtor3<-mean(mtor3$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_mtor3)
Sortdatalo<-filter(Sortdata,S < average_mtor3)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((mtor_3_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((mtor_3_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc1,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc1@meta.data$Mtor_3<-mtor3$Average-control$Average

Glutamine_1<-read.table("I:/MAC/Memory_B cells/GO_BP_Glutamine.txt",sep = "\t",header = TRUE)
Glutamine_1_new<-Glutamine_1[2:length(rownames(Glutamine_1)),]
Glutamine_1_new<-tolower(Glutamine_1_new)
library(Hmisc)
Glutamine_1_new<-capitalize(Glutamine_1_new)
Glutamine1<-FetchData(pbmc1,vars = Glutamine_1_new)
for (i in 1:length(rownames(Glutamine1))) {
  Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
}

Total<-FetchData(pbmc1,vars = rownames(pbmc1))
Data<-data.frame()
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_Glutamine1<-mean(Glutamine1$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_Glutamine1)
Sortdatalo<-filter(Sortdata,S < average_Glutamine1)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((Glutamine_1_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((Glutamine_1_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc1,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc1@meta.data$Glutamine_1<-Glutamine1$Average-control$Average

remove(control,Data,Data_T,Glutamine_1,Glutamine1,Glutamine_1_new,mtor_1,mtor_1_new,mtor_2,mtor_2_new,mtor_3,mtor_3_new,mtor1,mtor2,mtor3)
remove(Total,Sortdata,Sortdatahi,Sortdatalo,rownamesDataT,all.genes,average_Glutamine1,average_mtor1,average_mtor2,average_mtor3,higene,logene,i,j,k,p)


S1<-data.frame(0,0,0,0)
colnames(S1)<-c("j","p","q","r")
for (j in 5:50) {
  for (p in c(50,100,200,300,400,500,600,700,800,900,1000)) {
    for (q in c(20,30,40,50,60,70,80,90,100)) {
      for (r in seq(0.5,10,by=0.1)) {
        
        pbmc1 <- FindNeighbors(pbmc1, dims = 1:j,n.trees = p,k.param = q)
        pbmc1 <- FindClusters(pbmc1, resolution = r)
        
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
        
        MTOR_1<-FetchData(pbmc1,vars = "Mtor_1")
        idents<-Idents(pbmc1)[colnames(pbmc1)]
        MTOR_1$ident<-idents
        noise<-rnorm(n=length(x=MTOR_1[,"Mtor_1"]))/1e+05
        MTOR_1[,"Mtor_1"]<-MTOR_1[,"Mtor_1"]+noise
        MTOR_1$ident<-factor(MTOR_1$ident,levels = names(x=rev(x=sort(x=tapply(X=MTOR_1[,"Mtor_1"],INDEX = MTOR_1$ident,FUN = mean),decreasing = FALSE))))
        
        GLUTAMINE_1<-FetchData(pbmc1,vars = "Glutamine_1")
        idents<-Idents(pbmc1)[colnames(pbmc1)]
        GLUTAMINE_1$ident<-idents
        noise<-rnorm(n=length(x=GLUTAMINE_1[,"Glutamine_1"]))/1e+05
        GLUTAMINE_1[,"Glutamine_1"]<-GLUTAMINE_1[,"Glutamine_1"]+noise
        GLUTAMINE_1$ident<-factor(GLUTAMINE_1$ident,levels = names(x=rev(x=sort(x=tapply(X=GLUTAMINE_1[,"Glutamine_1"],INDEX = GLUTAMINE_1$ident,FUN = mean),decreasing = FALSE))))
        
        MTOR_2<-FetchData(pbmc1,vars = "Mtor_2")
        idents<-Idents(pbmc1)[colnames(pbmc1)]
        MTOR_2$ident<-idents
        noise<-rnorm(n=length(x=MTOR_2[,"Mtor_2"]))/1e+05
        MTOR_2[,"Mtor_2"]<-MTOR_2[,"Mtor_2"]+noise
        MTOR_2$ident<-factor(MTOR_2$ident,levels = names(x=rev(x=sort(x=tapply(X=MTOR_2[,"Mtor_2"],INDEX = MTOR_2$ident,FUN = mean),decreasing = FALSE))))
        
        MTOR_3<-FetchData(pbmc1,vars = "Mtor_3")
        idents<-Idents(pbmc1)[colnames(pbmc1)]
        MTOR_3$ident<-idents
        noise<-rnorm(n=length(x=MTOR_3[,"Mtor_3"]))/1e+05
        MTOR_3[,"Mtor_3"]<-MTOR_3[,"Mtor_3"]+noise
        MTOR_3$ident<-factor(MTOR_3$ident,levels = names(x=rev(x=sort(x=tapply(X=MTOR_3[,"Mtor_3"],INDEX = MTOR_3$ident,FUN = mean),decreasing = FALSE))))
        
        if (levels(Kdm6b$ident)[1] == levels(Irf4$ident)[1]){
          if (levels(Kdm6b$ident)[1] == levels(Ly75$ident)[1]){
            if(levels(Kdm6b$ident)[1] == levels(Myc$ident)[1]){
              
              if (levels(Myc$ident)[1] == levels(MTOR_1$ident)[1]){
                if (levels(Myc$ident)[1] == levels(GLUTAMINE_1$ident)[1]){
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
                  i<-length(rownames(Irf4_data_hi))/length(rownames(Irf4_data))
                  
                  j<-levels(Ly75$ident)[1]
                  Ly75_seurat<-subset(pbmc1,idents = j)
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
                  
                  if (c >= 0.5 & f >=0.5 & i>=0.5 & l >=0.5 & o >=0.5 ){
                    P1<-data.frame(j,p,q,r)
                    colnames(P1)<-c("j","p","q","r")
                    S1<-rbind(S1,P1)
                  
                  
                }
              }
            }
          }
        }
        
        
      }
    }
  }
  }
}
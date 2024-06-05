# Analysis
## Dimentional Reduction
### PCA
Upon the completion of the preprocessing step, The Ewing V5 Seurat obj has a following dimension;
```
> Ewing.su
An object of class Seurat 
46939 features across 6110 samples within 1 assay 
Active assay: RNA (46939 features, 2500 variable features)
 82 layers present: data.1, data.2, data.3, data.4, data.5, data.6, data.7, data.8, data.9, data.10, data.11, data.12, data.13, data.14, data.15, data.16, data.17, data.18, data.19, data.20, data.21, data.22, data.23, data.24, data.25, data.26, data.27, counts.1, scale.data.1, counts.2, scale.data.2, counts.3, scale.data.3, counts.4, scale.data.4, counts.5, scale.data.5, counts.6, scale.data.6, counts.7, scale.data.7, counts.8, scale.data.8, counts.9, scale.data.9, counts.10, scale.data.10, counts.11, scale.data.11, counts.12, scale.data.12, counts.13, scale.data.13, counts.14, scale.data.14, counts.15, scale.data.15, counts.16, scale.data.16, counts.17, scale.data.17, counts.18, scale.data.18, counts.19, scale.data.19, counts.20, scale.data.20, counts.21, scale.data.21, counts.22, scale.data.22, counts.23, scale.data.23, counts.24, scale.data.24, counts.25, scale.data.25, counts.26, scale.data.26, counts.27, scale.data.27, scale.data
```
The first step of a lengthy dimentional reduction process is to estimate PCA and this is a very simple process;
```Ewing.su<-RunPCA(merged.su, reduction.name = "PCA", reduction.key = "PCA_")```
There are a few ways to check and make sure PCA was run reasonably.  One is a heatmap of PCA;
```DimHeatmap(Ewing.su, reduction = "PCA", dims = 1:5, ncol = 3)```

![PCA_Heatmap](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/bcecf5b3-9da0-4f8f-87d1-eab6c606b2a8)

By default, Seurat generate 50  PCA dimensions. and above command checks the fist 5 dimentions. As dimetions go deeper, the distinction of high and low expressions of genes becomes ambiguous, as shown in the figure above.  
Another way to check the PCA result is to examine scatter 2D(or maybe 3D) plots of PCAs;

```
expand.grid(1:5, 2:5) %>% .[!(.$Var1==.$Var2), ] ->s
rownames(s)<-1:nrow(s)
multiDimPlot<-function(i){
  x<-s[i, "Var1"]
  y<-s[i, "Var2"]
  t<-DimPlot(merged.su ,dims =c(x,y), reduction = "PCA", group.by = "orig.ident")
} 
lapply(1:16, multiDimPlot) -> t
wrap_plots(t)+plot_layout(guides = "collect")
```

![PCA_DimPlot](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/238f9456-e44f-4f2b-be2e-5b7b7d03fb7d)

Take home message from a figure above is that you do not see much of sample seprations on PC5.  One thing that stands out is that the sample, **TM348** clearly stands out from the rest of the sample.   

PCA is really a key step for the rest of analysis proedures and an one critital parameter of PCA is a number of effective PCs. This could be derived from a visual inspection of so-called **"Elbow Plot"** or less subjective by using some R packages to find effective PCs statistically. In most cases, the elbow plot works just fine. 
```ElbowPlot(Ewing.su, reduction = "PCA", ndims = 50)```

![Elbow](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/939d7297-9406-44fc-bb81-15adaa41e575)

An above plot seems to show right after PC 30, a line of points is hiting a plateau, and thus, PC 30 could be used as a maximum number of dimetions for the rest of steps.  Arguably the higher or lower PC values (than PC30) could be used but practicaly, they will not make significant differences.

![PC30](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/79a2d3a5-f5d7-4d87-af6d-394e049813fe)

![PC49andPC50](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/8ac8048d-022d-4d2c-8717-789f076d2a6c)

A heatplot PCA at PC30 shows very vague distinction betwwen high and low expression genes but epression of some of these genes are still polarized to some cells.  However, at PC49 and PC50, there is absolutely no distinguishing features in any of the cells.  

There are some interesting things about the elbow plot.  Now, it is possible to fit a curve over these points;
```
Using a mgcv package's gam function to fit the curve over those points.  
library(mgcv)
ggplot(elbow$data, aes(dims, stdev))+geom_point()+geom_smooth(method = "gam")
```
![elbow gam](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/a9456bec-6ccf-4336-a652-f636da181363)

The fitted curve looks great.  Now,lets just predict values of ```stdev``` beyond PC 50;
```
elbow.gam<-gam(stdev ~ s(dims, bs="cs"), data = elbow$data)
data.frame(dims=seq(100, 1000, by=100))->df
predict(elbow.gam, newdata = df)

> predict(elbow.gam, newdata = df)
          1           2           3           4           5           6           7           8           9          10          11 
 1.13387599 -0.06106273 -1.25600146 -2.45094019 -3.64587891 -4.84081764 -6.03575637 -7.23069509 -8.42563382 -9.62057254  4.54860938
```
As number of ```dims``` increases, ```stdev``` decreases, and this is really expected, however, right around PC 100 and beyond, ```stdev``` gets negative values, and ```stdev``` never gets the negative value.  Strage.  

Once PCA is done, the rest of the step is pretty straightforward.  A following step could be done in any order but lets' figure out number of potential clusters.  The first thing to do is to figure out ```cell neighbors```.
```
Ewing.su<-FindNeighbors(Ewing.su, reduction = "PCA", dims = 1:30)
```
Then, find clusters among cells. 
```
Ewing.su<-FindClusters(merged.su, graph.name = "RNA_snn", resolution = c(0.7, 0.8, 1.0, 1.5), algorithm = 4, method = "igraph")
```
This will create  and assign cluster IDs to each cells.  Seurat stores the cluster info by creating new metadata columns, named ```RNA_snn_res``` as follows ;
```
> str(merged.su[[]])
'data.frame':	6110 obs. of  22 variables:
 $ orig.ident         : chr  "TM338" "TM338" "TM338" "TM338" ...
 $ nCount_RNA         : num  11769 4075 5127 2343 2263 ...
 $ nFeature_RNA       : int  3123 2028 2154 1105 963 2373 1010 1341 2273 1077 ...
 $ harmony_clusters   : chr  "18" "7" "7" "28" ...
 $ seurat_clusters    : Factor w/ 24 levels "1","2","3","4",..: 13 3 3 15 2 13 16 7 7 2 ...
 $ harmony_clusters_1 : chr  "13" "5" "5" "20" ...
 $ patientID          : chr  "ES-024" "ES-024" "ES-024" "ES-024" ...
 $ tumor_site         : chr  "Femur" "Femur" "Femur" "Femur" ...
 $ percentMT          : num  7.09 11.07 9.25 23.47 6.14 ...
 $ cell.sorted        : chr  "Live" "Live" "Live" "Live" ...
 $ sample.type        : chr  "Resection" "Resection" "Resection" "Resection" ...
 $ percentRB          : num  1.22 3.04 7.41 3.16 14.18 ...
 $ percentERCC        : num  2.93 7.75 9.3 9.43 27.49 ...
 $ scDblFinder.class  : chr  "singlet" "singlet" "singlet" "singlet" ...
 $ RNA_snn_res.0.7    : Factor w/ 20 levels "1","2","3","4",..: 7 3 3 7 1 7 17 8 8 1 ...
 $ RNA_snn_res.0.8    : Factor w/ 20 levels "1","2","3","4",..: 6 3 3 6 1 6 17 8 8 1 ...
 $ RNA_snn_res.1      : Factor w/ 20 levels "1","2","3","4",..: 6 3 3 6 1 6 17 8 8 1 ...
 $ RNA_snn_res.1.5    : Factor w/ 24 levels "1","2","3","4",..: 5 3 18 5 1 5 19 7 7 1 ...
```
Now, create UMAP;
```
Ewing.su<-RunUMAP(Ewing.su, dims = 1:30, reduction = "PCA", reduction.name = "unitegrated.umap")
```
It's time to see the UMAP of Ewing.su;
```
DimPlot(Ewing.su, reduction = "unitegrated.umap", group.by = "RNA_snn_res.1.5", label = F, label.size = 7)+scale_color_manual(values = P40)+guides(color = guide_legend(title = "Sorted Sample ", size=12, ncol = 2,override.aes = list(size = 5), theme = theme(legend.title = element_text(hjust = 0.5))))+labs(title = "Ewing UMAP Unitegrated_RNA_snn_res.1.5")+theme(plot.title=element_text(hjust = 0.5, size = 20))
```
![Umap unitegrated_clusters](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/2a076a6c-1ab7-4681-9421-864e6df4f38d)

Chainging ```group_by``` argument from "seurat_clusters" to "orig.ident" recated the same UMAP with different grouping of cells;

![unitegrated Umap](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/4a0eb9d9-5827-4fe0-b48b-3f0bf2a80e37)

These UMAPs are based off the unitegrated data across different samples/batches (e.g. TM786, ...etc) and to alleviate batch and sample effects, sequencing data should be integrated, so that the exact same cell types with undistinguishable gene expression profiles from distinct samples will fall into the same position hence, the cluster in UMAP.  Despite the fact that this UMAP is not based on the integrated data, except some samples, a majoriy of samples are "mixed" well in clusters.  There are better more obvious exmaples of UMAP based on the unitegrated data in the Seurat tutorial and elsewhwere. 

Integrating differnet layers (distinct batches and samples (in this case)) is pretty simple under a Seurat V5's new data integration pipe line.  
```
Ewing.su<-IntegrateLayers(Ewing.su, method = HarmonyIntegration, orig.reduction = "PCA", new.reduction="integrated.harmony")
```
It is used be a bit more involved in the previous versions but now the integration can be run in the one liner.  Also, the same line can be used to apply different integration methods (e.g. CCA, Harmony, RPCA, ...).  
You can run Harmony manually without using the Seurat's integration line.  All needed is the PCA.  
```
RunHarmony(Ewing.su[["PCA"]]@cell.embeddings, merged.su[[]], "orig.ident", early_stop=F,  lambda=NULL, plot_convergence = F, nclust=50)->tm1
```
An advantage  of this step is to run a quick diagnostic on the peformance of Harmony by setting ```early_stop=T``` and ```plot_convergence = TRUE```. This termiates ```RunHarmony`` at 10the cycle. A resulting plot looks like below;
![harmony](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/1fd64702-7439-4d20-9a4c-c9078b2dea1c)
As seen in the plot, after 6 cycle, integration seems to be converged to a single entity, meaning ```nclust=50``` running 50 cycles of integration steps will not necessary give a superior results. A manually generated Hamrnomy integration set is embedded into the Seurat obj as a reduction obj as follows;
```
Ewing.su[["harmony.pca"]]<-CreateDimReducObject(embeddings = "tm1", key = "harmonyPCA_", assay = DefaultAssay(Ewing.su))
```
After this, it is necessary to run ```FindNeighbors```, ```FindClusters```, and ```RunUMAP``` again, and plot the new UMAP based on the Harmony integration by Dimplot specifying a name of the DimRed created by ```RunUMAP```.  A following mod DimPlot will generate the UMAP below.  
```
DimPlot(merged.su, reduction = "integrated.harmony.umap", group.by = "orig.ident")+scale_color_manual(values = p40)+guides(col = guide_legend(
title = "Orig.Ident",
size=12,
ncol = 1,
override.aes = list(size = 5),theme = theme(legend.title = element_text(hjust = 0.5))))+
labs(title = "Ewing UMAP Integrated Orig.Ident")+
theme(plot.title=element_text(hjust = 0.5, size = 20))
```
![HamonryIntegratedUmap](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/45ec22c8-8fa3-4c83-b16a-5986d9dd5f1c)

How about a UMAP by the cluster ? 
![Harmony_Umap_with clusters](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/6af2a928-2be3-4d01-b9ec-ce1b44163736)
As noticed, a figure above is different from others created by Suerat's ```DimPlot```.  Actually, this figure is manually created by using ```ggplot2```and ```ggrepel```.  The reason why this was done is to simply to accentuate location of each clusters on the UMAP.  In many occasions, it is not always straightforward to see boundaries of clusters and **3D ball** shape helps to see the tangible pile-up of the point.  At any rate, the number indicates the cluster levels, and there are 23 clusters.  In the previous **unitegrated** UMAP, there are 24 clusters.  Integration definitely influences the neighboring and clustering procedures.  A script for creating this figure is as follows;
```
Ewing.su@reductions$integrated.harmony.umap@cell.embeddings %>% as.data.table(, keep.rownames = T) -> umap.dt #extract the umap embedding and create data.table
umap.dt[, .(median_x=median(integratedharmonyumap_1), median_y=median(integratedharmonyumap_2)), by=harmony.pca_snn_res.1.5] ->umap.median #median values will be used as UMAP coordinates for ggrepel text annotation
ggplot(umap.dt, aes(integratedharmonyumap_1,integratedharmonyumap_2))+
  geom_point(
    shape = 21,
    color = "black",
    aes(fill=harmony.pca_snn_res.1.5))+
  geom_label_repel(data=umap.median, 
                   aes(median_x , median_y, 
                       label = harmony.pca_snn_res.1.5),
                   fill="white", 
                   color="red",
                   max.overlaps = Inf,
                   arrow = arrow(length = unit(0.015, "npc"), type = "closed", ends = "last"),
                   size=3, 
                   min.segment.length = 0,
                   nudge_x =0, 
                   nudge_y = 2,
                   box.padding = 0.5)+
  scale_fill_manual(values = p40)+
  theme_classic()+
  labs(title = "Harmony Integration with Clusters")+ #change the title accordingly
  theme(plot.title = element_text(family = "Arial",hjust = 0.5, size = 14))+
  theme(legend.position = "none") #comment this out if a legend will be included
```
This pretty much concludes this sectioin.  The next section will be about ```differential gene expression```. 

## Differential Gene Expression
Differnetial gene expression (DGE) in a SueratV5 obj could be easily done but the ```layers``` of the obj must be joined before DGE is calculated.  Joining layers of the V5 obj is very simple;

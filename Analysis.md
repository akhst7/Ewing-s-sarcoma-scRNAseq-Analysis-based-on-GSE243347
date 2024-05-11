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

Right around PC 30, the plot seems to hit a plateau, which could be used as a maximum number of dimetions for the rest of steps.  The higher PC values could be used but practicaly, they will not make significant differences.  


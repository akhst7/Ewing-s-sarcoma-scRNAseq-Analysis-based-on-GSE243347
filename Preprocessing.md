# Preprocessing 
## Introduction
Ewing sarcoma is a form of a pediatric sarcoma that develops in the bone and surrounding soft tissue. In a majority of cases, it is initiated by chromosomal translocation of between chromosomes 11 and 22, t(11,22), resulting in the fusion of the Ewing Sarcoma Breakpoint Region 1 (EWSR1) and he Friend Leukemia Virus Integration 1 (FLI1) gene, resulting chimeric protein that drives tumor promotion (https://en.wikipedia.org/wiki/Ewing_sarcoma).   

There are quite a few NGS data on GEO site, a majority of which belong to the human vitro studies.  Apparently, there is a significant hurdle to develop a preclinical mouse model, simply because a EWS-FLI1 fusion protein is often lethal  (Methods Mol Biol. 2021:2226:183-189.  doi: 10.1007/978-1-0716-1020-6_14), and mouse do not naturally develop bone sarcoma that resemble Ewing. Also for human clinical studies, biopsy samples from tumor lesion may not be enough to do a relaible NGS study, which typically is the case of any clinical studies.  Nonetheless, as of 2024, there are more than dozens publications that have performed NGS from biopsies and I came across one of the recent study with scRNAseq data posted in GEO,  Visser LL, Bleijs M, Margaritis T, van de Wetering M et al. Ewing Sarcoma Single-cell Transcriptome Analysis Reveals Functionally Impaired Antigen-presenting Cells. Cancer Res Commun 2023 Oct 24;3(10):2158-2169. PMID: 37823774.  GEO accession # is GSE243347

As its title suggests, this study looks at a cross sectional state of Ewing sarcoma single cells transcriptomes based on 18 Ewing sarcoma biopsies from 11 patients.  I am sure the biopsy samples were minuscule and was pain staking processes to generate clean viable single cell suspension.  The authors employed CEL-Seq2 rather than 10x V3.  At the end of cell extraction and purification step, they did not have enough cells to run the pricy V3 10X library prep methods, which require many more cells.  
Anyway, a following is the abstract; 
>Novel therapeutic strategies are urgently needed for patients with high-risk Ewing sarcoma and for the reduction of severe side effects for all patients. Immunotherapy may fill this need, but its successful application has been hampered by a lack of knowledge on the composition and function of the Ewing sarcoma immune microenvironment. Here, we explore the immune microenvironment of Ewing sarcoma, by single-cell RNA sequencing of 18 Ewing sarcoma primary tissue samples. Ewing sarcoma is infiltrated by natural killer, T, and B cells, dendritic cells, and immunosuppressive macrophages. Ewing sarcoma–associated T cells show various degrees of dysfunction. The antigen-presenting cells found in Ewing sarcoma lack costimulatory gene expression, implying functional impairment. Interaction analysis reveals a clear role for Ewing sarcoma tumor cells in turning the Ewing sarcoma immune microenvironment into an immunosuppressive niche. These results provide novel insights into the functional state of immune cells in the Ewing sarcoma tumor microenvironment and suggest mechanisms by which Ewing sarcoma tumor cells interact with, and shape, the immune microenvironment.

Now, GSE243347 has 27 samples as SRA (fastq files) and two supplemental files;

<img width="334" alt="image" src="https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/6d89ed13-b89e-4925-9d8d-e46dcae9ef97">

Two supplemental files have the processed data; ```GSE243347_RAW.tar``` is a tar archive of count matrix and ```GSE243347_data_annotation_txt.gz``` contains meta data.  
There are 27 gziped count matrices in ```GSE243347_RAW.tar```;

![image](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/fcf04175-50a0-4000-8293-b1e87f1ef2b8)

## Data loading
These count matrices are imported to R and used to create a Seurat V5 object.  This process needs a bit of imagination but it is not as complicated, owing to the advanced feature implemented in the Seurat V5.  A Seurat v5 object can store a series of individual matrix as a layer and exist as an one giant file with multiple  merged layers; which makes handling and processing of a large pile of data. 
```
library(data.table)
library(fs)
tmp<-lapply(dir_ls("/GEO/GSE243347_RAW/"), function(x) fread(x, header = T))
names(tmp)<-list.files("/GEO/GSE243347_RAW/")

library(Seurat)
#making sure a Seurat default assay is V5
options("Seurat.object.assay.version")
$Seurat.object.assay.version
[1] "v5"

suobjfromcount<-function(x){
y<-x[, -1] #using x as a counter
z<-y[, 1:ncol(y):=lapply(.SD, as.integer), .SDcols = 1:ncol(y)]
tmp.mx<-as.matrix(z, rownames = x$GENEID) #text files must be coverted to matrix prior to creating a Seurat obj
w<-CreateSeuratObject(counts = tmp.mx)
return(w)
}
su.list<-lapply(tmp, suobjfromcount) #creating individual Seurat objs
su.merged<-merge(su.list$TM338, su.list[-1], add.cell.ids=names(su.list), project="Ewing") 
su.merged
An object of class Seurat 
48257 features across 10395 samples within 1 assay 
Active assay: RNA (48257 features, 0 variable features)
 27 layers present: counts.1, counts.2, counts.3, counts.4, counts.5, counts.6, counts.7, counts.8, counts.9, counts.10, counts.11, counts.12, counts.13, counts.14, counts.15, counts.16, counts.17, counts.18, counts.19, counts.20, counts.21, counts.22, counts.23, counts.24, counts.25, counts.26, counts.27
Layers(su.merged) #you should see all the matrices imported.  

Layers(su.merged)
 [1] "counts.1"  "counts.2"  "counts.3"  "counts.4"  "counts.5"  "counts.6"  "counts.7"  "counts.8"  "counts.9"  "counts.10" "counts.11" "counts.12" "counts.13" "counts.14" "counts.15"
[16] "counts.16" "counts.17" "counts.18" "counts.19" "counts.20" "counts.21" "counts.22" "counts.23" "counts.24" "counts.25" "counts.26" "counts.27"

Layers(su.merged) %>% length()
[1] 27

dim(su.merged)
[1] 48257 10395
```
A merged seurat obj has 27 samples with total of combined 48275 genes and 10395 cells;
```
rownames(su.merged) %>% head()
[1] "ENSG00000000003--TSPAN6"   "ENSG00000000419--DPM1"     "ENSG00000000457--SCYL3"    "ENSG00000000460--C1orf112" "ENSG00000000938--FGR"      "ENSG00000000971--CFH" 

colnames(su.merged) %>% head()
[1] "TM338_UNK" "TM338_A1"  "TM338_A2"  "TM338_A3"  "TM338_A4"  "TM338_A5" 
```
Interestingly, authors pasted Ensembl gene ids to gene names, which is sort of unorthodox approach but I think a reason why they did this to avoid having duplicated gene names/symbol,  which Suerat will most definitely complain. Column name indicates the cell ids.   Descriptions of all the ids are well documented in a companion meta data file but briefly, cell ids are composed of sample id_well# corresponding the single cell library, and there are 385 cells sorted in to 385 wells individually  per sample as shown below; 

<img width="600" alt="image" src="https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/70541c21-3af6-4749-a549-437da4e2b926">

The next  is a quality control step where low quality cells and presumptive doublets are removed by mainly using Bioconductor packages, ```Scuttle and ScDblFinder```.   For this, su.merged must be first converted to a ```SingleCellExperiment``` obj ,and then take advantage of ```scuttle``` ’s functionality for using **median absolute deviation (MAD)** of RNA counts for easy and reliable filtering.  Typically,  MAD value is set to 3 which serves as a lower (as well as an upper) cutoff value, followed by a  ```scDblFinder``` package to identify potential doublets with high high accuracy.  Additional quality control step is performed possibly by the spike-in control denoted as ERCC.  Also, the spike-in control could potentially be used to normalize counts among different samples (wells).  For bulk RNAseq,  it is somewhat controversial to use the spike-in normalization.  

****Creating SingleCellExperimemnt obj****

```Layers(Ewing.su)[31:54] -> layer
lapply(layer, FUN = function(x){LayerData(Ewing.su, layer = x)})-> list
tmp<-do.call(cbind, args = list)
sce<-SingleCellExperiment(assays=list(counts=tmp))
ERCC<-sce[grep("ERCC", rownames(sce)), ] #subsetting an ERCC SingleCellExperiment obj from sce
altExp(sce, “ERCC")<-ERCC #incorporating the ERCC SingleCellExperiment obj into the altExp slot in sce
sce
class: SingleCellExperiment 
dim: 48257 10395 
metadata(0):
assays(2): counts logcounts
rownames(48257): ENSG00000000003--TSPAN6 ENSG00000000419--DPM1 ... ENSG00000284452--RP11-625N16.1 ENSG00000284522--RP11-444J21.6
rowData names(0):
colnames(10395): TM338_UNK TM338_A1 ... TM770_P23 TM770_P24
colData names(16): cell.ID orig.ident ... sizeFactor sizeFactor.ERCC
reducedDimNames(0):
mainExpName: NULL
altExpNames(1): ERCC
```
Just how many spike-in controls are there ? 

```
sce[grep("ERCC", rownames(sce)), ] %>% rownames()
  [1] "ENSG00000012061--ERCC1"         "ENSG00000049167--ERCC8"         "ENSG00000104884--ERCC2"         "ENSG00000134899--ERCC5-cluster" "ENSG00000163161--ERCC3"        
  [6] "ENSG00000175595--ERCC4"         "ENSG00000182150--ERCC6L2"       "ENSG00000225830--ERCC6"         "ERCC-00002"                     "ERCC-00003"                    
 [11] "ERCC-00004"                     "ERCC-00009"                     "ERCC-00014"                     "ERCC-00019"                     "ERCC-00022"                    
 [16] "ERCC-00025"                     "ERCC-00028"                     "ERCC-00033"                     "ERCC-00034"                     "ERCC-00035"                    
 [21] "ERCC-00039"                     "ERCC-00040"                     "ERCC-00042"                     "ERCC-00043"                     "ERCC-00044"                    
 [26] "ERCC-00046"                     "ERCC-00051"                     "ERCC-00053"                     "ERCC-00054"                     "ERCC-00058"                    
 [31] "ERCC-00059"                     "ERCC-00060"                     "ERCC-00062"                     "ERCC-00067"                     "ERCC-00069"                    
 [36] "ERCC-00071"                     "ERCC-00073"                     "ERCC-00074"                     "ERCC-00076"                     "ERCC-00077"                    
 [41] "ERCC-00078"                     "ERCC-00079"                     "ERCC-00084"                     "ERCC-00085"                     "ERCC-00092"                    
 [46] "ERCC-00095"                     "ERCC-00096"                     "ERCC-00097"                     "ERCC-00099"                     "ERCC-00104"                    
 [51] "ERCC-00108"                     "ERCC-00111"                     "ERCC-00112"                     "ERCC-00113"                     "ERCC-00116"                    
 [56] "ERCC-00120"                     "ERCC-00123"                     "ERCC-00126"                     "ERCC-00130"                     "ERCC-00131"                    
 [61] "ERCC-00134"                     "ERCC-00136"                     "ERCC-00138"                     "ERCC-00143"                     "ERCC-00144"                    
 [66] "ERCC-00145"                     "ERCC-00148"                     "ERCC-00150"                     "ERCC-00154"                     "ERCC-00157"                    
 [71] "ERCC-00158"                     "ERCC-00160"                     "ERCC-00162"                     "ERCC-00163"                     "ERCC-00164"                    
 [76] "ERCC-00165"                     "ERCC-00170"                     "ERCC-00171"                     "ENSG00000134899--ERCC5"         "ENSG00000186871--ERCC6L"       
 [81] "ERCC-00013"                     "ERCC-00017"                     "ERCC-00031"                     "ERCC-00086"                     "ERCC-00137"                    
 [86] "ERCC-00147"                     "ERCC-00109"                     "ERCC-00142"                     "ERCC-00156"                     "ERCC-00168"                    
 [91] "ERCC-00057"                     "ERCC-00117"                     "ERCC-00041"                     "ERCC-00024"                     "ERCC-00016"                    
 [96] "ERCC-00075"                     "ERCC-00012"                     "ERCC-00098"                     "ERCC-00081"                     "ERCC-00083"                    
[101] "ERCC-00061" 
```
So there are 101 spike in controls they threw in and the size factor of ERCC in each sample are used as a cutoff to remove unwanted cells.  In addition, % mitochondria (MT) transcripts (per cell) needs to be determined prior to removing low quality cells but how many MT transcripts are there in the data set ? 
```
rownames(sce)[grep(glob2rx("*MT-*"), rownames(sce))]
 [1] "ENSG00000198695--MT-ND6"          "ENSG00000198712--MT-CO2"          "ENSG00000198727--MT-CYB"          "ENSG00000198763--MT-ND2"          "ENSG00000198786--MT-ND5"         
 [6] "ENSG00000198804--MT-CO1"          "ENSG00000198840--MT-ND3"          "ENSG00000198886--MT-ND4"          "ENSG00000198886--MT-ND4-cluster"  "ENSG00000198888--MT-ND1"         
[11] "ENSG00000198899--MT-ATP6"         "ENSG00000198899--MT-ATP6-cluster" "ENSG00000198938--MT-CO3"          "ENSG00000209082--MT-TL1"          "ENSG00000210049--MT-TF"          
[16] "ENSG00000210082--MT-RNR2"         "ENSG00000210100--MT-TI"           "ENSG00000210107--MT-TQ"           "ENSG00000210112--MT-TM"           "ENSG00000210127--MT-TA"          
[21] "ENSG00000210140--MT-TC"           "ENSG00000210144--MT-TY"           "ENSG00000210154--MT-TD"           "ENSG00000210194--MT-TE"           "ENSG00000210195--MT-TT"          
[26] "ENSG00000210196--MT-TP"           "ENSG00000211459--MT-RNR1"         "ENSG00000212907--MT-ND4L"         "ENSG00000228253--MT-ATP8"         "ENSG00000210151--MT-TS1"         
[31] "ENSG00000210117--MT-TW"           "ENSG00000210191--MT-TL2"          "ENSG00000210176--MT-TH"           "ENSG00000210164--MT-TG"           "ENSG00000210184--MT-TS2"         
[36] "ENSG00000210077--MT-TV"           "ENSG00000210156--MT-TK"           “ENSG00000210135--MT-TN
```
As expected, a majority of these genes encode mitochondrial tRNA
```
library(AnnotationHub)
ah<-AnnotationHub()
ah
query(ah, c("Ensembl", "Homo sapiens"))
Ensembl111Hs<-ah[["AH116291"]]
grep(glob2rx("*MT-*"), rownames(sce), value = T) %>% str_extract(., pattern="\\w+") ->mt
genes(Hs, filter = ~ gene_name=="^MT-" ) ->grl
df<-as.data.frame(grl)
df[, c("gene_name", "description")]
```

<img width="546" alt="image" src="https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/e100bad1-077f-4484-98f9-fea594a8a0a8">

Next, % mitochondrial transcript in each cells is determined  The same could be done for ribosome transcripts ;
```
mitogenes<-grep(glob2rx("*MT-*"), rownames(sce), value = T) #this could also be done by Ewing.su
rownames(sce)[grep("ENSG\\d+--RP[SLF]", rownames(sce))] -> rbgenes #this could also be done by Ewing.su
cbind(colData(sce), perCellQCMetrics(sce, subsets=list(Mito=mitogenes, Ribo=rbgene))) -> colData(sce)
```
You can do the same with Seurat;
```
pattern <- "ENSG\\d+--MT" # perl=TRUE
Erwing.su<-PercentageFeatureSet(Erwing.su, pattern = pattern, col.name = "percentMT")
PercentageFeatureSet(Ewing.su, pattern = "ENSG\\d+--RP[SLF]", col.name = "percentRB") ->Ewing.su
```
Now it is time to examine these parameters. . Typically, scatter plots of parameters such as transcript counts and types are employed.
```
In Seurat;
FeatureScatter(Ewing.su, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")+scale_color_manual(values = p30)
FeatureScatter(Ewing.su, feature1 = "nCount_RNA", feature2 = "percentMT", group.by = "orig.ident")+scale_color_manual(values = p30)
```
![image](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/6684a83e-543c-4188-85ed-fd6648006da3)
![image](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/4e0a7807-c3ed-47d5-b2d1-e402e9ba600f)

There is a  large number of “low quality” (low RNA count and RNA features) cells and cells with the high MT content (anything above 20-40% but this have to be decided case by case).  Visualizing the MT content of cells from each samples could be done by a Violin plot. 
```
VlnPlot(Ewing.su, features = "percentMT", group.by = "orig.ident")+NoLegend()
```
![image](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/effa1744-c520-491d-9b67-ace3e5866af5)
As seen in the above violin plot, a majority of cells in all the samples do not exhibit extreme %MT values which are under 40%, suggesting a few dead cells to be removed.  based on this criteria. .  However, this based solely on the MT content for low quality cell filtration may be a bit too naive.  By looking at a violin plot of RNA feature, it is clear that cells in TM768 have quite a low RNA feature counts  and this suggests that quality of these cells may be low enough to be removed, though again this must be confirmed by the further QC analysis.
```
VlnPlot(Ewing.su, features = "nFeature_RNA", group.by = "orig.ident")+NoLegend()+scale_y_continuous(name = "nFeature_RNA_Log10", transform = "log10")
```
![image](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/6e63b315-4e48-4bee-8648-a0db7b67e589)

As mentioned above, the spike-in-control (ERCC) can provide more clear threshold of low quality cells.  It is generally accepted that the high % ERCC is associated less desirable low quality cells.  Since ERCC parameters are stored in the sce obj, the ```scater``` and ```scuttle``` packages are used to visualize %ERCC  and generate MAD values.
```
plotColData(sce, y = "percentERCC", x = "orig.ident")
```
![image](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/951d4b1a-072a-481c-8530-b94603f43fee)

It appears that a majority of cells in TM768 are extremely low quality. If a majority of cells are to be  removed, contributions of TM768 as a group certainly become insginificant.  It may be easier to remove all the cells in TM768 but at this point, they will be included.  
If 3 MAD is used as a cutoff value, almost a half of the cell  will be removed, and only 4768 cell are retained. 
```
discard.ercc <- isOutlier(colData(sce.filtered.1), type = 'higher', nmads = 3)
sce$discard.ercc.selected.nmad3<-discard.ercc
plotColData(sce, y = "percentERCC", x = "orig.ident", color_by = "discard.ercc.selected.nmad3")+theme(axis.text.x =element_text(angle = 45, vjust = 1, hjust = 1 ))+ggtitle("MAD with Selected Batch Correction")+theme( plot.title = element_text(hjust = 0.5, size = 20))
```
![Rplot01](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/263dcadb-40b6-4c15-8c38-87554b07a9e4)

This simple approach basically is to determine outliers that are beyound 3 MAD.  Noticeably, the same cutoff %ERCC is applied to all the sample groups, although MAD values should be quite different among some groups by looking at the distribution of cells in each sample group.  Evidently, ```isOutlier``` without a particular argument does not distinguish batches/samples, thus all the cells will be merged into the one batch/sample.  To obtain more accuate batch/sample specific MAD values, ```batch=``` argument must be supplied.   
```
discard.ercc <- isOutlier(sce$altexps_ERCC_percent,type="higher", batch=sce$orig.ident, nmads = 3)
```
![Rplot 2](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/36015148-a397-4de9-81b7-c1d31d78c694)

The figure above shows potential outliers only in certain samples, particularly the samples that have relatively small proportion of cells with high %ERCC, indicating that sequecing in  a majority of cells in these samples are more than adequate. In contrast, other samples show distribution of cells with less than desirable %ERCC, yet these cells are not cosidered **outliers**. This happens because ```isOutlier``` assumes all the samples contain high quality cells.  In reality, this is not always the case.   To fix this uneven identification of outliers, ```isOutlier``` requires yet another argument, ```subset=``` which selects samples with visibly **cleaner** distributions of cells lower %ERCC, in other words, samples with high proportion of high quality cells, in this case, samples, ```"TM416", "TM417", "TM505", "TM548","TM549", "TM574", and "TM736."```  

```
discard.ercc <- isOutlier(sce$percentERCC,type="higher", batch=sce$orig.ident, nmads = 3, subset = sce$orig.ident %in%  c("TM416", "TM417", "TM505", "TM548","TM549", "TM574", "TM736"))
```
![Rplot02](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/355f37e6-26c6-4890-835d-944f65208cf7)


This reuslts in 4768 cells. Removing almost half of the total cells may be too excessibe, even though it sounds statiscally reasonable.  According to the first author, 9,537 wells(cell equivalnet) were recovered after the ERCC QC step, which is almost double the number of cells done here.  The discrepancy can be explained by the fact that the first author used the raw ERCC counts rather than %ERCC(normalized) and employed ```The Sharq pipeline``` https://www.biorxiv.org/content/10.1101/250811v2.full.  I have absolutely no issues with them using ```The Sharq pipeline``` which is designed particularly for flow soreted well based scRNAseq and reviewers did not have any issues, yet cells with high values of %ERCC or even raw counts suggest that substancial numbers of cells needed tobe removed.  The ERCC qc filter that they used seems a lot softer and MAD 3 filterraton demostrated here seems a bit too extreme in removeing cells.  

I came up with unscientific yet visually acceptable cutoff valuee based on %ERCC using MAD 3.  By looking at a plot of RNA features and counts, a majority of low quality cells (TRUE) cells appear to be piling up at very low end of scales, agreeing with a fact that low quality cells typically have relatively lower number of RNA counts and types with higher %ERCC or ERCC counts.  

![Rplot03](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/38c95b32-3e13-4457-90d1-fda9a7c7edb7)

Now, a violin plot of RNA features or counts  vs MAD 3 reject criterium shows a little hump at the bottom  of "FALSE" cells and a large disc at the bottom of the TRUE cells with a clear boundary.  If line is drwawn, a corresponding RNA feature of an intercept of the line is 500 features.  

![Rplot 3](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/fbbe5436-2c6e-4bb9-a418-74853ac26b10)

Using RNA feature less than 500 as a filter criterium removes roughly 4000 cells and now the total number of the cell is 6434, which is acceptable for now. 

```
sce.500<-sce[, !sce$nFeature_RNA <=500]
Ewing.su[, !Ewing.su$nFeature_RNA <= 500] ->Ewing.500.su
```

A next step is to remove doublets.  As mentioned above, doublet removal is done by the biocondactor pacakge, ```scDBlFinder```.  

```
sce.500<-scDblFinder(sce.500, samples = sce.500$orig.ident, BPPARAM = bp)
```
```scDblFinder``` comes with a ton of customerization but there are a few very important ones.  Those are listed in https://www.bioconductor.org/packages/release/bioc/html/scDblFinder.html 
You can play with these parameters if you know what you are doing, thhough as suggested, default parameters of ```scDBlFinder`` works fine.  Of note, a ```sample``` parameter instructs ```scDBlFinder``` to treat the data file as a merged multiple sample file, while multiple samples increase the accuracy and precision of detecting doublets.  A ```bp``` paremter informs ```scDBlFinder``` to utilize a multi-thread or multi-core capbility of the computer. 

```
library(BiocParallel)
MulticoreParam(18, RNGseed = 1234, progressbar = T)-> bp
```
There are 5 metadata column created by ```scDBlFinder```, of which the most relevant is ```scDblFinder.class```.  A column under ```scDblFinder.class``` denotes singlet vs doublet classification, and used to filter out doublets.  Values in other columns are for optimization of doublet detection. 

```
sce.500.nodbl<-sce.500[, sce.500$scDblFinder.class=="singlet"] or sce.500.nodbl<-sce.500[, !sce.500$scDblFinder.class=="doublet"]

For Seurat:
colData(sce.500)[,18]-> tmp
names(tmp)<-colnames(sce.500)
Ewing.500.su<-AddMetaData(object = Ewing.500.su, metadata = tmp, col.name = "scDblFinder.class")
Ewing.500.nodbl.su<-Ewing.500.su[, Ewing.500.su$scDblFinder.class == "singlet"]
```
A scatter plot of RNA features vs counts show some odd cells. These cells have excessively large RNA features and/or counts.  For some reasons, these cells were not considered to be doublets, however, these will be removed for the rest of the analysis for obvious reasons.  

![Rplot 5](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/f9099606-1ab8-49cd-a25e-46c85384378a)

```
`%notin%` <- Negate(`%in%`)

For SingleCellExperiment:
sce.500.nodbl[, colnames(sce.500.nodbl) %notin% c("TM506_A18", "TM338_UNK", "TM506_UNK")] ->sce.500.nodbl.nooddcells
For Seurat
Ewing.500.nodbl[, colnames(Ewing.500.nodbl) %notin% c("TM506_A18", "TM338_UNK", "TM506_UNK")] ->Ewing.500.nodbl.nooddcells
```
Beacause of the curse of **S3 obj**, the simplest approach, ```sce.500.nodbl.su[, !c("TM506_A18", "TM338_UNK", "TM506_UNK")] or Ewing.500.nodbl.su[, !c("TM506_A18", "TM338_UNK", "TM506_UNK")] does not work. 

The preprocessing does not end here. It tuned out that the data contained a significant number of duplicates, which happened to be the "lncRNA, and rRNA"  in the category of **Ensembl biotype** known as **misc_RNA**.  But, the first thing is to take a look at gene names.  
```
genename<-rownames(Ewing.500.nodbl.nooddcells)
length(genename)
[1] 48156
```
Now, names of genes are composed of the Ensembl ID and gene symbol, for example  "ENSG00000000003--TSPAN6".  This is rather unorthodox gene names since in most cases, there are separate slots for the Ensembl ID and the symbol.  I guess authors did this to avoid gene symbol duplication but this type of designating gene symbol is troubling for downstream analysis particularly cell annotation and DEGs. In fact, in trial DEGs, there is gene named, **"Y-RNA"** showed up multiple times in the list go top20 DEGs per clusters. It turns out that "Y-RNA" is,  **ncRNA** and in Ensembl, is in the category of **Misc_RNA**.  These genes are non-coding genes including **lncRNA** that play roles in expression regulation of protein coding genes.  As in the example below, it has an Ensemble gene designation, and there are many of them but exactly how many are there in the data set ?
![Screenshot 2024-05-07 at 1 37 36 PM](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/096429a6-05c8-4ac8-9a5b-b0d4fba60b83)

Within the total genes registered in Ewing.500.nodbl.nooddcells, 48156, **"Y-RNA"** count is 583, meaning each of these has **"Y-RNA"** as gene symbols but has distinct Ensembl IDs.  Besides **"Y-RNA"**, there are other 133 genes similar to **"Y-RNA"**, belonging to Ensembl's Misc_RNA and most of them are **lncRNA**;  
```
which(str_extract(genename, pattern) %>% duplicated()) %>% genename[.] %>% str_extract(., pattern) %>% unique() -> genename.unique
```
A ```genename.unique``` is a collection of unique gene symbols of all the duplicates; 
```
> genename.unique
  [1] "Y-RNA"           "SNORA2"          "SNORA72"         "SNORA26"         "Metazoa-SRP"     "RGS5"            "U3"             
  [8] "snoU13"          "ACA64"           "SNORA31"         "SNORA67"         "SNORA73"         "PDXDC2P"         "LINC01481"      
 [15] "U1"              "BMS1P4"          "CYB561D2"        "COG8"            "RNU11"           "U6"              "SNORD45"        
 [22] "SNORA27"         "SNORA25"         "SNORA4"          "SNORA51"         "SNORA40"         "DNAJC9-AS1"      "SPATA13"        
 [29] "LINC01598"       "LINC01422"       "SNORD63"         "GOLGA8M"         "EMG1"            "POLR2J4"         "uc-338"         
 [36] "SNORA63"         "SNORA7"          "SNORA24"         "SNORD38"         "SNORA48"         "SNORD113"        "SNORA79"        
 [43] "LINC00484"       "SNORD59"         "SCARNA20"        "SOGA3"           "SFTA3"           "LINC01238"       "LINC01115"      
 [50] "7SK"             "U7"              "DLEU2-2"         "DGCR5"           "RABGEF1"         "C2orf81"         "SNORA9"         
 [57] "U8"              "SNORD46"         "SNORD42"         "SNORA62"         "SNORA19"         "LINC01347"       "SNORA11"        
 [64] "snoU109"         "SNORA58"         "SCARNA21"        "SNORD56"         "SNORD74"         "SNORD67"         "NBPF13P"        
 [71] "SCARNA15"        "SNORA81"         "FAS-AS1"         "MATR3"           "SNORA75"         "SNORD112"        "SNORA64"        
 [78] "SNORD70"         "SNORD11"         "RPS23P5"         "SNORA33"         "SNORA20"         "SNORA74"         "U4"             
 [85] "SNORD116"        "SNORA22"         "SNORA8"          "SNORA42"         "5S-rRNA"         "SNORA3"          "SNORD77"        
 [92] "MAL2"            "C2orf61"         "SNORA70"         "SNORD81"         "SNORA69"         "SNORA12"         "SNORA77"        
 [99] "SNORD75"         "SNORA57"         "SCARNA18"        "SCARNA6"         "SCARNA16"        "OR7E47P"         "DLEU2-6"        
[106] "DLEU2-1"         "DLEU2-5"         "RAET1E-AS1"      "SNORD30"         "U2"              "SNORA43"         "SNORA30"        
[113] "SNORA1"          "SNORA68"         "SNORD5"          "SNORA35"         "SCARNA17"        "SNORD39"         "snoMe28S-Am2634"
[120] "PROX1-AS1"       "SNORA38"         "SNORD66"         "snoU2-19"        "SNORD65"         "SCARNA11"        "HYMAI"          
[127] "SNORA76"         "SNORA46"         "5-8S-rRNA"       "ALG1L9P"         "CTSLP2"          "SNORA18"         "Vault"
```
A followin plot shows the extent of duplication of these;
```
sapply(genename.unique, function(x){grep(x, rownames(Ewing.join.noERCC.su)) %>% length()})->genename.unique.counts
data.table(counts=genename.unique.counts, name=names(genename.unique.counts))->tmp
ggplot(tmp)+geom_bar(aes(x=reorder(name, counts, decreasing=T), y=counts), stat="identity")+theme(axis.text.x = element_text(face="bold", angle = 90, size = 5, hjust = 0.95, vjust = 0.2))+scale_y_continuous(expand = c(0, 0))
```
![Rplot01 4](https://github.com/akhst7/Ewing-s-sarcoma-scRNAseq-Analysis-based-on-GSE243347/assets/3075799/20586a90-8f95-4440-8ae8-0571f1398b11)
Again, as mentioned earlier, these are not quite multiplicates of the same gene  but rather they are various isforms of the same gene; for instance, there are approximately 800 isoforms of U6RNA.  For the downstream analysis, particularly the DGE analysis, they will predominate over low abundant protein conding RNAs of significance.  They are not  ***noises*** in a biological sense and they may well have some biological importance.  However, it may be best to analyze them separately from the protein conding RNAs.

Before removing these genes from the latest Ewing Seurat obj, ```Ewing.500.nodbl.nooddcells``` there are a few protein coding RNA, mRNAs in ```gene.unique```, and this gets a bit complicated as follows;
```
library(AnnotationHub)
ah<-AnnotationHub()
> ah
AnnotationHub with 71308 records
# snapshotDate(): 2024-04-29
# $dataprovider: Ensembl, BroadInstitute, UCSC, ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/, Haemcode, FANTOM5,DLRP,IUPHAR,HPRD,STRING,SWISSPROT,...
# $species: Homo sapiens, Mus musculus, Drosophila melanogaster, Rattus norvegicus, Bos taurus, Pan troglodytes, Danio rerio, Gallus gallus...
# $rdataclass: GRanges, TwoBitFile, BigWigFile, EnsDb, Rle, OrgDb, SQLiteFile, ChainFile, TxDb, Inparanoid8Db
# additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
#   rdatapath, sourceurl, sourcetype 
# retrieve records with, e.g., 'object[["AH5012"]]' 

             title                                 
  AH5012   | Chromosome Band                       
  AH5013   | STS Markers                           
  AH5014   | FISH Clones                           
  AH5015   | Recomb Rate                           
  AH5016   | ENCODE Pilot                          
  ...        ...                                   
  AH116725 | TENET_consensus_open_chromatin_regions
  AH116726 | TENET_consensus_promoter_regions      
  AH116727 | ENCODE_dELS_regions                   
  AH116728 | ENCODE_pELS_regions                   
  AH116729 | ENCODE_PLS_regions
```
I am looking at an Ensemble database;
```
> query(ah, c("Ensembl", "Homo sapiens"))
AnnotationHub with 287 records
# snapshotDate(): 2024-04-29
# $dataprovider: Ensembl, BioMart, FANTOM5,DLRP,IUPHAR,HPRD,STRING,SWISSPROT,TREMBL,ENSEMBL,CELLPHONEDB,BADERLAB,SINGLECELLSIGNALR,HOMOLOGE...
# $species: Homo sapiens, homo sapiens
# $rdataclass: GRanges, TwoBitFile, EnsDb, data.frame, SQLiteFile, list, OrgDb
# additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
#   rdatapath, sourceurl, sourcetype 
# retrieve records with, e.g., 'object[["AH5046"]]' 

             title                                  
  AH5046   | Ensembl Genes                          
  AH5160   | Ensembl Genes                          
  AH5311   | Ensembl Genes                          
  AH5434   | Ensembl Genes                          
  AH5435   | Ensembl EST Genes                      
  ...        ...                                    
  AH113665 | Ensembl 110 EnsDb for Homo sapiens     
  AH113914 | LRBaseDb for Homo sapiens (Human, v006)
  AH116291 | Ensembl 111 EnsDb for Homo sapiens     
  AH116542 | LRBaseDb for Homo sapiens (Human, v007)
  AH116710 | org.Hs.eg.db.sqlite
```
Ensembles 111 is the latest database for human genes, and getting EnsID and gene symbol can be down as follows;
```
select(ensHu111, keys = GenebiotypeFilter("protein_coding"), columns = "SYMBOL") -> protein.name
head(protein.name)
   SYMBOL    GENEBIOTYPE
1         protein_coding
2    A1BG protein_coding
3    A1CF protein_coding
4     A2M protein_coding
5   A2ML1 protein_coding
6 A3GALT2 protein_coding
> nrow(protein.name)
[1] 19458
```
There are almost 20,000 mRNA entries in the latest EnsDb, and to find those shared by ```protein.name``` and ```genename.unique``` is doen by a following script;
```
intersect(genename.unique, protein.name$SYMBOL) ->s
> s
 [1] "RGS5"     "CYB561D2" "COG8"     "SPATA13"  "GOLGA8M"  "EMG1"     "SOGA3"    "SFTA3"    "RABGEF1"  "C2orf81"  "MATR3"    "MAL2"
```
12 coding genes are present in ```genename.unique``` and they must be removed.  
```
genename.unique[which(genename.unique %notin% s)] ->u
genename[which(str_extract(genename, pattern) %notin% u)] ->v
Ewing.join.noERCC.su[rownames(Ewing.join.noERCC.su) %in% v, ] ->Ewing.join.noERCC.protincoding.su 
```

The preprocessing step is now finished.  One final comment in this section is that authors did not go through these steps described here.  The preprocessing step described here is somewhat more rigrous than what they did in their publications, however, this is by no means that their preprocessig step is inadequate.
Of note, ribosome RNAs are deleted from their downstream analysis.  This was intentionally done because like mitochondrial transcripts, quite a number of ribosome transcripts are significantly more abundant than coding transcipts which otherwise interfere with and have a significant impact on a following DEG analysis.  Although this may be a common practice for some groups, I would suggest otherwise.  Arguably, not all ribosome transcripts are housekeeping and are by no means ubiquitously expressed; they do exibit the element of cell type specific expression. Thus,  importance of ribosome transcripts as DEGs should not be disregarded. 

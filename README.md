# omcplr
An R package for identifing genomic regions harboring complex genomic rearrangements and deconstructing their hierarchical organization based on the optical genome mapping technology.

## Getting Started

We recommend an R version 3.6.0 or more recent in order to use the package.

## Introduction
Genomic instability and mutagenesis are hallmarks of cancer, and most tumors have rearranged genomes. Most genomic rearrangements in tumors are structurally simple, but complex genomic rearrangements (CGR) i.e. genomic reorganization events involving multiple intertwined structural rearrangements are common in all major cancer types. We utilized optical mapping to identify genomic rearrangements in cancer genomes, including those that are complex. we developed omcplR to resolve the CGR architectures and predict their genesis using the genome-graph framework. 

## Availability and Installation
The development version of omcplr package is available at https://github.com/sjdlabgroup/omcplr and can be installed as
```
install.packages("devtools")
devtools::install_github("sjdlabgroup/omcplr",build_vignettes = FALSE )
```
## omcplr workflow

Step 1. Aligned each mapID contig with the reference genome and binned the mapping mapID contig with a resolution of 10kb 

Step 2. Retrieve all mapID contigs that mapping to the same seed bin and clustered mapID contigs and defined a complex haplotype as having 
at least four mapID contigs to the same reference chromosome position.

Step 3. Utilized the BNLearn algorithm to trace the origin and trajectory of the segment junctions.

```

## Step 1. Aligned each mapID contig with the reference genome and binned the mapping mapID contig with a
resolution of 10kb 
    
Details

Aligned each mapID contig with the reference genome and binned the mapping mapID contig with a resolution of 10kb,
utilizing the same version of the reference genome as used in the bionano analysis. The reference genome was binned
with a 10kb resolution. The output includes all the bins that are covered by each mapID, along with the
corresponding mapIDs that have overlapping regions in the reference genome. The output consists of two columns:
the first column displays the included reference bins, and the second column lists the mapIDs associated with
those reference bins.
    
Examples
 
library(GenomicRanges)
library(data.table)
genome_bin <- read.table("data/genome.10kb.bed",header=T)
xmap <- read.table("data/xmap",header=T)
combine_tmp <- overlap_bin_mapID(genome_bin,xmap,window_size=10000)

head(combine_tmp)
     V1        type_merge
1  C1_1 26572,26622,17811
4  C1_2 26572,26622,17811
7  C1_3 26572,26622,17811
10 C1_4 26572,26622,17811
13 C1_5 26572,26622,17811
16 C1_6 26572,26622,17811
 
## Step 2. Retrieve all mapID contigs that mapping to the same seed bin and clustered mapID contigs and defined
a complex haplotype as having at least four mapID contigs to the same reference chromosome position.
    
Details
 
Retrieve all mapID contigs that mapping to the same seed bin and clustered mapID contigs and defined a complex
haplotype as having at least four mapID contigs to the same reference chromosome position.

Examples

library(stringr)
complex_event<-cplx_events(combine_tmp,4)

head(complex_event$blocks[[1]])
[1] "C1_1634" "C1_1635" "C1_1636" "C1_1637" "C1_1638" "C1_1639"

head(complex_event$contigs[[1]])
[1] "152"   "172"   "22361" "23532" "23541" "26061"

## Step 3. Utilized the BNLearn algorithm to trace the origin and trajectory of the segment junctions.

Details

Complex genomic regions are frequent occurrences that drive the development of cancer. However, identifying
somatic complex genomic regions has proven to be a challenging task. By employing Optical mapping, structural
variations (SVs) can be successfully detected on the same DNA fiber, indicating their tandem occurrence on a
single/mutiple chromosomes. This enabled to pinpoint the segment junctions between these tandem SVs. During
the investigation using BNG analysis, we observed the presence of multiple heterogeneous DNA contigs,
referred to as complex regions (inferred by step 2), all aligning to the same location on the reference genome.
We were able to extract all the segment junctions from these regions, and these junctions appeared to be closely
associated with the initiation of large-scale rearrangement events across multiple chromosomes.

i) extract segment junctions between tandem mapID

Examples

library(GenomicRanges)
library(dplyr)
library(stringr)
library(gtools)
library(tidyr)

xmap <- read.table("data/xmap")
cov <- read.table("data/xmap_coverage")
colnames(cov)<-c("mapID","cov")

event=26  
single_complex_event<-complex_event$contigs[[event]]
single_complex_df<-xmap[which(xmap$V1 %in% single_complex_event),] 
single_complex_junction<-extract_junctions(single_complex_df)

head(single_complex_junction)
  mapID              note
1 22671 C7_14420|C7_14429
3 22671 C7_14430|C7_14438
4 22672 C7_14423|C7_14429
5 22672 C7_14420|C7_14429
6 22672 C7_14428|C7_14438
7 22672 C7_14432|C7_14438

ii) add coverage data for each mapID and generate trajectory plot

#extract junction for each complex event,The vector 'single_complex_junction' must not be empty

single_complex_junction<-separate(single_complex_junction, col=note, into=c('junc1', 'junc2'), sep='\\|')
colnames(single_complex_junction)<-c("V1","V2","V3")

#extract coverage for each complex event

single_complex_cov<-unique(cov[which(cov$mapID %in% unique(single_complex_junction$V1)),c("mapID","cov")])

#extract bins and mapIDs for each complex event

single_complex_combine_tmp<-overlap_bin_mapID(genome_bin,single_complex_df, window_size=10000)
one_complex_event<-cplx_events(single_complex_combine_tmp,4)
res<-order.of.events(cplx_list=one_complex_event,bpdf=single_complex_junction,covdf=single_complex_cov,
cplx.event=1, optimization="pc.stable")

#make trajectory plot
traj_plot <- ooe.plot(res)
```

## Other functions 

### junction_circos_plot

#### Description
Created circos plot highlighting the most prevalent junctions
```
library(reshape2)
library(pals)
library(circlize)
library(RColorBrewer)
junction_circos_plot(one_complex_ta_cov,link_cut_off = 2)
```
#### Value
one_complex_ta_cov is a dataframe with "mapID junction coverage"
link_cut_off define the link color, if shared junctions are less than 4, the color define as lightgrey

#### Example
```
one_complex_ta_cov <- read.table("data/one_complex_ta_cov.txt",sep='\t',header=T)
junction_circos_plot(one_complex_ta_cov,link_cut_off = 2)
```
### coor_dataframe

#### Description
convert junction bin to a dataframe including the first and second junction positions
```
coor_dataframe(coordinate,window=10000)
```
#### Value
coordinates is the junction bins
window is the bin size

#### Examples
```
junction <- "C7.14419kb.C7.14432kb"

coor_dataframe(junction,window=10000)
[[1]]
   chr    start      end
1 chr7 14419000 14429000

[[2]]
   chr    start      end
1 chr7 14432000 14442000
```
 
### complex_block_annotation

#### Description
After finding the interesting complex regions, one can proceed to map them to known genes and annotate these bins 
with complex events
```
library(tidyr)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
complex_block_annotation(block,window_size=10000,USCS_gene="TxDb.Hsapiens.UCSC.hg38.knownGene")
```
#### Value
block is a dataframe with the complex region bin column
window is the bin size
USCS_gene is the reference gene that is utilized 

#### Examples
```
block<-read.table("data/complex_block.txt")
block_gene<-complex_block_annotation(block,window_size=10000,USCS_gene="TxDb.Hsapiens.UCSC.hg38.knownGene")
block_gene
[1] "TPK1"      "ARHGEF35"  "ARHGEF34P"
```

### combine_rows

#### Description
Merge rows that have any shared elements
```
combine_rows(data, sep=';')
```
#### Value
data is a dataframe with shared element columns
sep is the seperator for each element and can be modified as needed

#### Examples
```
df <- read.table("data/df_combine_row.txt")      
combine_rows(df,sep=';')
```
### merge_range

#### Description
merged SVs that satisfy two criteria: (i) if both breakends of the SVs fall within a 5kb range, and (ii) the SVs 
share the same type, they are combined into a single SV
```
merge_range(data)
```
#### Value
data is a dataframe with columns "chr1 start1 end1 chr2  start2 end2 BNGTYPE type"

#### Example
```
merge_data<-read.table("data/merge_range.txt")
after_merge_range<-do.call(rbind,merge_range(merge_data))
```







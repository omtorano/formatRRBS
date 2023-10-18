# formatRRBS
This is a package to format RRBS methylation data from bismark files for multi-omic integration or other downstream processing. The package contains two functions import_RRBS(), see R/import_RRBS.R, for importing bismark files into R, and format_RRBS(), see R/format_RRBS.R, for formatting individual files into a data matrix to be used for multi-omics integration or other application. Regular usage of format_RRBS() assumes RRBS feature annotation to the EPA fathead minnow genome GCF_016745375.1_EPA_FHM_2.0 https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016745375.1/ For notes on other usage see Formatting for other species. 

# Import
import_RRBS() has three parameter inputs

- ctrl: Character(s) denoting control treatment in file names
- trts: List of character(s) denoting treatment groups in file names
- dir: Directory with bismark.cov.gz files

Example usage
```
import_RRBS("FL-0", c("ML-0","ML-10"), dir = "RRBS/methylation_coverage")
```
Where RRBS/methylation_coverage contains the following files:
FL-0-1_bismark_bt2.bismark.cov.gz
FL-0-2_bismark_bt2.bismark.cov.gz
ML-0-1_bismark_bt2.bismark.cov.gz
ML-0-2_bismark_bt2.bismark.cov.gz
ML-10-1_bismark_bt2.bismark.cov.gz
ML-10-2_bismark_bt2.bismark.cov.gz

The directory parameter dir will default to the current working directory. The characters denoting control and treatment groups must be unique up to the character(s) denoting replicate. In this example n=2 for control and each treatment groups. The characters "FL-F"are sufficient to classify files 1 & 2 as control, "ML-0" classify files 2 & 3 as treatment 1, and "ML-10" classify files 4 & 5 as treatment 2. Note in this example including the proceeding "-" e.g. "FL-F-" would also correctly classify the samples, but it is not necessary and omission makes for cleaner file names in the final output. The entries for ctrl and trts parameters are saved in the global environment upon running import_RRBS() and are required for format_RRBS(). 
The imported bismark files will be saved to the global environment as a list called "x". Users should note that this will overwrite any existing variable of the same name (x) if already present in the global environment. This list will have n data frame elements equal to the number of files imported. Each data frame contains the contents of the individual bismark files, 6 columns by n rows equaling the number of methylation loci. 
The columns in each data frame are as follows: 1) chromosome/scaffold, 2) start of methylation location, 3) end of methylation location, 4) percent methylation, 4) methylated C's, and 5) unmethylated C's. For information about the generation of bismark.cov.gz files see the bismark users guide https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf. 
image
# Format
format_RRBS() has five parameter inputs
- x: List of data frames containing bismark.cov.gz file information, generated during import_RRBS()
- meta: A metadata file with rownames matching names of x and a TRT column containing treatment group information
- total_count: Minimum counts (methylated + unmethylated) at each feature, default = 10
- coverage: Feature coverage within and across treatment levels required for retainment, default = 0.8 (i.e. features not present in at least 80% of samples are removed)
- window: Window of interest in bases upstream and downstream of genes, requires two values, default =  1000, 1000
Examples of usage
```
format_RRBS(x, meta = meta, coverage = 0.8, window = c(1000, 1000)) #default parameters
format_RRBS(x, meta = proj_x_meta, coverage = 1, window = c(2000, 1000)) #will return values for features present in all samples falling between 2000 bp upstream and 1000 bp downstream of genes
```
## Notes on metadata 
Only the rownames and TRT column of metadata will be considered. Metadata can contain other columns but none will be used. The rownames of meta must match the names of x, the list of data frames containing RRBS info. This can be confirmed using the following code, meta should be replaced with the name of your metadata.
```
names(x) == rownames(meta)
```
The TRT column should contain values representing the treatment group of the associated samples. This usage of "TRT" was chosen to match the MIB differential gene expression pipeline.
## Notes on total counts
Filtering out methylation sites with a low number of total counts, methylated + unmentylated, is recommended. In edgeR the conservative rule of thumb is for each site to have a minimum of 8 counts, see https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf. Our default is slightly more conservative at 10.
## Notes on coverage
The features present in each sample varies greatly. Since we are interested in concatenating samples into a single data frame for downstream analysis, it is necessary to limit features to those shared between samples. The default parameters of the format_RRBS() function will limit features to those shared between 80% of samples. This is carried out so that the distribution of coverage is equal between treatment groups. This is done to ensure no features are kept that are absent from a given treatment group. 
## Notes on window
The window of interest is defined as base pairs upstream and downstream of a gene window. Setting the window of interest to 0,0 would therefore limit features to those only falling within genes. There is currently no functionality that allows the exclusion of the gene window. Gene annotations are stored internally (see data-raw/DATASET.R) and are parsed from GCF_016745375.1_EPA_FHM_2.0_genomic.gtf.gz
## Function output
Running format_RRBS() will result in four files saved to the current working directory of the R session. The concatenated RRBS data will be stored in a .RDS file, the contents of this file are a data frame with rows corresponding to features (chromosome-loci) and columns corresponding to sample. The values of the data frame will be either beta or M-values as noted in the file name. Beta values represent "percent methylation" at a given methylation loci and are calculated with the following equation: (100*(M/(M+U)). These values follow a beta distribution. M-values: (log2((M+1)/(U+1)) are unbounded and follow a Gaussian distribution. Using these values for downstream applications including differential methylation and multi-omics integration is recommended. 
The third output is a stats.csv file containing information from each filtering step, specifically the number of features in each raw bismark.cov.gz file, features remaining after total_count filtering, features remaining after % coverage filtering, features remaining after window filtering, and features remaining after filtering non-variable methylation loci (features with the same methylation state across all samples). If the number of features is zero after percent_coverage filteration choosing a less stringent percent coverage is recommended.
The fourth output is a Unique_Features_to_Genes.csv file which provides gene information for each unique methylation feature included in analysis. This file has five columns, the first two correspond to the chromosome and location of each methylation feature, column four is the gene associated with the methylation feature, and columns four and five are the gene start and end locations on the chromosome. For a methylation feature to be associated with a gene it must fall within the specified basepair window, i.e. (chrom_end + window[2)] >= loc >=(chrom_start - window[1]). Note therefore that these associations will vary depending on the specified window. In addition, if the specified window is greater than the spacing of genes along the chromosome this location will be associated with the downstream gene. This will be modified in v2 of this package.
# Formatting for other species
import_RRBS() will work for any bismark files, regardless of study species. format_RRBS() assumes methylation features will be annotated to the EPA fathead minnow genome GCF_016745375.1_EPA_FHM_2.0.
If you are working with human data please see https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf for importing and analyzing data. The edgeR RRBS import and formatting methods differ from formatRRBS in several ways. First, edgeR readBismark2DGE() concatenates all features across samples. This necessitates the user filtering for percent coverage after concatenation. This is not ideal for fathead minnow (and potentially other non-human species) because of the variability in features between samples. Concatenating all features from all samples results in data frames with millions of features with large amounts of missing data, requiring large amounts of storage. Second, edgeR readBismark2DGE() assumes numeric chromosome annotations, which is not true for fathead minnow scaffolds. Character chromosomes/scaffolds will not throw errors but will be misinterpreted at line 38 of readBismark2DGE().
If you are working with non-human data there are two recommended methods to make format_RRBS() functional. The first is to make a copy of the formatRRBS function, edit data-raw/DATASET.R line 2 by replacing the FHM gtf.gz with your species of interest, and run the use_data() command on line 8. As long as the replacement gtf follows standard format, this replacement should allow running format_RRBS() as described above. devtools::load_all('path/to/formatRRBS/') will load formatRRBS as a function, see https://r-pkgs.org/introduction.html for more details on making & editing R packages.
The second option is to copy lines 12-173 of format_RRBS() and run to save as a function in your global environment. Then copy lines 2-3 of data-raw/DATASET.R replacing the FHM gtf.gz with your species of interest, run to save in your global environment. Then format_RRBS() can be run as a function outside the formatRRBS package.

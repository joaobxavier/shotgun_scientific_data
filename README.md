# shotgun_scientific_data
 Code from the Scientific Data manuscript to publish the shotgun metagenomics of allo-HCT patients. We also demonstrate functional study of those data such as antibiotic resistance genes. 

## - Figures
#### Figure1/scFigure1.m
 
<ul>
 <li>Create the t-SNE plot calculated from Bray-curtris dissimilarity matrix using 16S and plot the samples with metagenome shotgun sequencing.  </li> 
  <ul> * CAUTION: running this function requires >30min! A saved result is available in 'savedMat' so it can be loaded directly. A parameter named 'calculateTSNE' is set to 1 when the script will recalculate the t-SNE score; else the script will load the saved data.</li>
 </ul>
 <li>Compares read counts in different stool consistency. </li>
</ul> 

#### Figure2/scFigure2.m

<ul>
   <li>Plot and compare the stool composition between 16S and metagenome of the samples from a single patient. </li> 
   <ul> *CAUTION: This script load the output from kraken2 of patric database and create a .mat file for each sample. This step is time consuming and takes >1 hour. When this is done, the variable 'rewriteShotgunAbundances' can be set to a value != 1, the script will read from a saved csv file. 
 </ul>
    </li> 
    <ul>*NOTE: This script does not include the (U)nclassified reads from PATRIC output.
   </ul>
   <li> Compare the relative abundance in different stool consistency of each taxa. </li>
</ul>

#### Figure3/scFigure3.m
<ul>
   <li>Calculate and plot the correlation between 16S and metagenomic sequencing; 
   <li>calculate the alpha diversity using Shannon Index and compare among stool consistency.  </li>
</ul>

#### Figure4/scFigure4.m
<ul>
   <li>Plot the vanA PCR result in the t-SNE plot (the same as in Figure 1) using the saved .mat file. 
   <li>Compare the relative abundance of vanA gene in the PCR(+) and PCR(-) groups. Examine the correlation between vanA  and vanB gene in shotgun metagenomes. </li>
</ul> 

#### Figure5/scFigure5.m
<ul>
   <li>Plot the phylogenetic tree built from Enterococcus isolated from stool of a HCT patient and metagenome assembled genomes from the same samples.  </li>
</ul> 

#### makeCARDtbl/scCreateCardTable2.m
<ul>
   <li>Make the output cardTbl.csv file.  </li>
</ul> 

#### makeVFDBtable/scProcessVFDB.m
<ul>
   <li>Make the output vfdbTbl_2021.csv file.  </li>
</ul> 


### - Metagenome data
#### metagenome_data/
<ul>
 Containing .csv files used for metagenome analysis
</ul>

#### savedMat/
<ul>
Containing saved data during figure generation to avoid re-calculation. This direcotry contains the results of t-SNE score calculated from Bray-curtris dissimilarity matrix using 16S, the correlation between 16S and shotgun, and the Shannon Index. 
</ul>

#### deidentified_data_tables/
<ul>
 Files from previous 16S data paper used in this study for data comparison. The tblASVsampes from <a href="https://github.com/liaochen1988/MSKCC_Microbiome_SD2021_Scripts">our other project</a> with an additional columns, AccessionShotgun, which is the SRA accession for the shotgun fastq files (only available for 395 samples).</ul>

#### PATRIC_output/
<ul>
Directories containing the kraken2, CARD and VFDB output of all shotgun samples by PATRIC. Only the files used for data analysis was included due to the limit of size limitation. </ul>

#### PATRIC_output_10samples/
<ul> 
 Directories with full output of the kraken2 for 10 samples, including kraken2, CARD and VFDB, and a .txt file containing the name of the 10 samples. This folder provides a complete view of the patric output and allows to try data analysis with a small sized data output. </ul>

## - utils
<ul>
 To run our scripts smoothly, the following functions are in use: 
 </ul>
 <ul> - violinplot: Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project
https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.4559847
 </ul>

 <ul> - beeswarm: https://github.com/ihstevenson/beeswarm
 </ul>
  <ul> - hex2rgb: https://www.mathworks.com/matlabcentral/fileexchange/46289-rgb2hex-and-hex2rgb
 </ul>
 <ul> - distinguishable_colors: https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
 </ul>
 <ul> - brewermap: https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
 </ul>
 <ul> - bh-tsne: https://lvdmaaten.github.io/tsne/

 </ul>
## - color-legends_16S.pdf
<ul>
  Color legends used for the major taxa in the 16S gene and shotgun metagenome sequencing.
 </ul>

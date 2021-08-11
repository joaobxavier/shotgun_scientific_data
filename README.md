# shotgun_scientific_data
 Code from the Scientific Data manuscript to publish the shotgun metagenomics of allo-HCT patients. We also demonstrate how those data can be used to study antibiotic resistance genes using PATRIC. Figure1-5 are directories containing the matlab code scripts to generate the figures in our metagenome shotgun manuscript (href=).
 
<ul>
  <li><b>Figure1:</b> scFigure1_tsne.m: Create the t-SNE plot calculated from Bray-curtris dissimilarity matrix using 16S and plot the samples with metagenome shotgun sequencing. The script also compares read counts in different stool consistency. (CAUTION: running this function requires >30min! And a saved result is available in 'savedMat' so it can be loaded directly. A parameter named 'calculateTSNE' is provided in the matlab code. When calculateTSNE==1, the script will recalculate the t-SNE score; when calculateTSNE != 1, the script will load the saved data.). </li>
 
   <li><b>Figure2:</b> scFigure2.m: Plot and compare the stool composition between 16S and metagenome of the samples from a single patient. This script load the output from kraken2 of patric database and create a .mat file for each sample. This step is time consuming and takes >30 min. When this is done, the variable 'rewriteShotgunAbundances' can be set to a value != 1, the script will read from a saved excel file. The relative abundance is also compared in each taxa. </li>
 
   <li><b>Figure3:</b> scFigure3.m: Calculate and plot the correlation between 16S and metagenomic sequencing; calculate the alpha diversity using Shannon Index and compare among stool consistency.  </li>
 
   <li><b>Figure4:</b> scFigure4.m: Plot the vanA PCR result in the t-SNE plot the same as in Figure 1 using the saved .mat file. Compare the relative abundance of vanA gene in the PCR(+) and PCR(-) groups. Examine the correlation between vanA  and vanB gene in shotgun metagenomes.  </li>
 
    <li><b>Figure5:</b> scFigure5.m: Plot the phylogenetic tree built from Enterococcus isolated from stool of a HCT patient and metagenome assembled genomes from the same samples. The tree was output by PATRIC using 1000 genes.  </li>
 
  <li><b>savedMat:</b> Containing some saved data during figure generation to avoid re-calculation. This also containing the results of t-SNE score calculated from Bray-curtris dissimilarity matrix using 16S, which requires >30 min to run in matlab.  </li>
 
   <li><b>PATRIC_output:</b> Directories with copies of the kraken2, CARD and VFDB analyses ran by PATRIC. Only the files used for data analysis was included due to the limit of size allowed to upload to github. </li>
 
  <li><b>PATRIC_output_10samples:</b> Directories with output of full report of the kraken2, CARD and VFDB analyses ran by PATRIC for 10 metagenome samples. A .txt file containing the name of samples was provided. This folder provides a complete view of the patric output and allows to try data analysis with a small sized data output. </li>
 
  <li><b>compareShotgunAnd16S:</b> Has Matlab script to import Kraken2 of each sample, create a .mat file in each subdirectory (to avoid having to import again), and compare the Kraken2 result with the 16S result. Running this code requires the files from <a href="https://github.com/liaochen1988/MSKCC_Microbiome_SD2021_Scripts">our other project</a>.</li>
  <li><b>tempFiles:</b> Directory usefull to keep files needed temporarily while running scripts.</li>
   importKraken2Output.m
  <li><b>importKraken2Output.m:</b> function to parse Kraken2 output into a Matlab table.</li>
  <li><b>scExampleImportKraken2Output.m:</b> An example script to import Kraken2 into Matlab.></li>
  <li><b>tblASVsamplesUpdatedWithShotgun.csv:</b> The tblASVsampes from <a href="https://github.com/liaochen1988/MSKCC_Microbiome_SD2021_Scripts">our other project</a> with an additional columns, AccessionShotgun, which is the SRA accession for the shotgun fastq files (only available for ### samples).</li>
   <li><b>updateFigshareSampleTableWithShotgunAccessionNumbers.m:</b> Script used to create the files tempFiles/tblASVsamplesWithShotgun04072021.csv and tempFiles/tblIsolateSamples04072021.csv.</li>
</ul>

# shotgun_scientific_data
 Code from the Scientific Data manuscript to publish the shotgun metagenomics of allo-HCT patients. We also demonstrate how those data can be used to study antibiotic resistance genes using PATRIC.
<ul>
  <li><b>PATRIC_output:</b> Directories with copies of the kraken2, CARD and VFDB analyses ran by PATRIC. Only the files used for data analysis was included due to the limit of size allowed to upload to github. </li>
  <li><b>PATRIC_output_10samples:</b> Directories with output of full report of the kraken2, CARD and VFDB analyses ran by PATRIC for 10 metagenome samples. A .txt file containing the name of samples was provided. This folder provides a complete view of the patric output and allows to try data analysis with a small sized data output. </li>
  <li><b>Figure1-5:</b> Directories containing the matlab code scripts to generate the figures in our metagenome shotgun manuscript (href=). </li>
  <li><b>savedMat:</b> Containing some saved data during figure generation to avoid re-calculation. This also containing the results of t-SNE score calculated from Bray-curtris dissimilarity matrix using 16S, which requires >30 min to run in matlab.  </li>
  <li><b>compareShotgunAnd16S:</b> Has Matlab script to import Kraken2 of each sample, create a .mat file in each subdirectory (to avoid having to import again), and compare the Kraken2 result with the 16S result. Running this code requires the files from <a href="https://github.com/liaochen1988/MSKCC_Microbiome_SD2021_Scripts">our other project</a>.</li>
  <li><b>tempFiles:</b> Directory usefull to keep files needed temporarily while running scripts.</li>
   importKraken2Output.m
  <li><b>importKraken2Output.m:</b> function to parse Kraken2 output into a Matlab table.</li>
  <li><b>scExampleImportKraken2Output.m:</b> An example script to import Kraken2 into Matlab.></li>
  <li><b>tblASVsamplesUpdatedWithShotgun.csv:</b> The tblASVsampes from <a href="https://github.com/liaochen1988/MSKCC_Microbiome_SD2021_Scripts">our other project</a> with an additional columns, AccessionShotgun, which is the SRA accession for the shotgun fastq files (only available for ### samples).</li>
   <li><b>updateFigshareSampleTableWithShotgunAccessionNumbers.m:</b> Script used to create the files tempFiles/tblASVsamplesWithShotgun04072021.csv and tempFiles/tblIsolateSamples04072021.csv.</li>
</ul>

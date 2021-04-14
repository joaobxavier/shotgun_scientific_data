# shotgun_scientific_data
 Code from the Scientific Data manuscript to bpuclish the shotgun metagenomics of allo-HCT patients. We will also demonstrate how those data can be used to study antibiotic resistance genes using PATRIC.
<ul>
  <li><b>PATRIC_output:</b> Directories with copies of the kraken2, CARD and VFDB analyses ran by PATRIC.</li>
  <li><b>compareShotgunAnd16S:</b> Has Matlab script to import Kraken2 of each sample, create a .mat file in each subdirectory (to avoid having to import again), and compare the Kraken2 result with the 16S result. Running this code requires the files from <a href="https://github.com/liaochen1988/MSKCC_Microbiome_SD2021_Scripts">our other project</a>.</li>
  <li><b>tempFiles:</b> Directory usefull to keep files needed temporarily while running scripts.</li>
   importKraken2Output.m
  <li><b>importKraken2Output.m:</b> function to parse Kraken2 output into a Matlab table.</li>
  <li><b>scExampleImportKraken2Output.m:</b> An example script to import Kraken2 into Matlab.></li>
  <li><b>tblASVsamplesUpdatedWithShotgun.csv:</b> The tblASVsampes from <a href="https://github.com/liaochen1988/MSKCC_Microbiome_SD2021_Scripts">our other project</a> with an additional columns, AccessionShotgun, which is the SRA accession for the shotgun fastq files (only available for ### samples).</li>
   <li><b>updateFigshareSampleTableWithShotgunAccessionNumbers.m:</b> Script used to create the files tempFiles/tblASVsamplesWithShotgun04072021.csv and tempFiles/tblIsolateSamples04072021.csv.</li>
</ul>

# Ndt80_LLR

Files and scripts associated with preprint 2023-03-12:

Meiotic prophase length modulates Tel1-dependent DNA double-strand break interference 

Luz M. López Ruiz1, Rachal M. Allison, William H. Gittens, Dominic Johnson, George Brown, and Matthew J. Neale1

Genome Damage and Stability Centre, School of Life Sciences, University of Sussex, UK

For correspondence: m.neale@sussex.ac.uk; L.Lpez-Ruiz@sussex.ac.uk


ABSTRACT 

During meiosis, genetic recombination is initiated by the formation of many DNA double-strand breaks (DSBs) catalysed by the evolutionary conserved topoisomerase-like enzyme, Spo11, in preferred genomic sites known as hotspots. DSB formation activates the Tel1/ATM DNA damage responsive (DDR) kinase, causing local inhibition of Spo11 activity in the chromosomal vicinity of the activating DSB—a process referred to as DSB interference. Intriguingly, in S. cerevisiae, over short genomic distances (<15 kb), Spo11 activity displays characteristics of concerted activity or clustering, wherein the frequency of DSBs that arise in adjacent hotspots is greater than expected by chance. We have proposed that clustering may arise due to a limited number of sub-chromosomal domains becoming primed for DSB formation.

Here, we demonstrate that whilst the process of broad-range DSB interference is retained, the characteristics of DSB clustering that arise between adjacent hotspots at short range is abolished when meiotic prophase timing is extended via deletion of the Ndt80 transcription factor. We propose that the extension of meiotic prophase in this manner, therefore, enables most cells, and therefore most chromosomal domains within them, to reach an equilibrium state of similar Spo11-DSB potential, reducing the impact that priming of a subpopulation of subchromosomal domains has on estimates of coincident DSB formation. Consistent with this view, genome-wide maps of Spo11-DSB formation generated in the absence of TEL1 are skewed towards regions that load pro-DSB factors early— revealing regions of preferential priming—but this effect is abolished when NDT80 is deleted.

Collectively, our work reinforces the view that the stochastic nature of Spo11-DSB formation in individual cells within the limited temporal window of meiotic prophase is the cause of localised DSB clustering—a phenomenon that is exacerbated in tel1∆ cells due to the dual role that Tel1 plays in the propagation of both DSB interference and meiotic prophase checkpoint control.

HOTSPOT TABLES:
Averages of the biological replicates, plus hotspot template used to measure hotspot-specific signals in each library:
Hotspot_Template_0.193V2_UPDATE_3473HS.txt
0.125_Hotspot.Table.spo11_yf_NorDNA.txt
V2.193_Hotspot.Table.sae2D_12A357BC.txt
V2.193_Hotspot.Table.sae2Dndt80D_D1D2TC10TC17.txt
V2.193_Hotspot.Table.sae2Dndt80Dtel1D_D1D2TC5TC10TC17.txt
V2.193_Hotspot.Table.sae2Dtel1D_12A357BC.txt

R SCRIPTS:
Averaging_FullMap_tables_V1: This script averages individual FullMap replicates into a combined FullMap where the sum of HpM equals 1 million.

Calculating background reads_V1: This script measures the percentage of signal registered within the 50 largest genes—regions of presumed Spo11 inactivity—on the S. cerevisiae genome as an estimate of the background noise.

Hotspot_analysis_V1: This script performs pairwise comparisons between datasets to study the degree of overlap, specificity and density of the identified hotspots (venn diagram and histograms).

Hotspot_identification_V1: This script identifies position and length of hotspots on single or multiple Spo11-DSB libraries. The total HpM signal is smoothed with a 201 Hann window. A cutoff of 0.1 HpM is then applied to remove the background noise. Hotspots are defined setting a minimum length of 25 bp and a minimum number of reads of 25. Hotspots separated by < 200 bp are merged and considered as a single hotspot. Hotspots are defined in each library separately and then combined to produce a single hotspot template that defines the position of every hotspot identified on the libraries.

Hotspot_Smooth_ratios_V1:This script calculates and represents the hotspots' fold changes between two Libraries (NormHpM ratio). 
 
Hotspot_table_V1: This script calculates the HpM signal included within each hotspot. Detailed description of the term reference list included in X.

NormHpM_V1: This script generates DSB maps representing the position and frequency of hotspots (NormHpM or NormHpChr). 

Pearson_correlation_V1:This script analyses the correlation between the hotspots strength of different datasets (NormHpM and NormHpChr Pearson correlation).

Ratio_heatmaps_V1: This script calculates and represents the hotspots' fold changes between two libraries (NormHpM) at 50 kb bin intervals on a per chromosome base ranked by chromosome size and centred at the centromere. 

Spo11 mapping Totals_V1: This script represents the position and frequency of the Spo11-DSBs signal (Total HpM) along the chromosome.

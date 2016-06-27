# Peak-assignment
Assign metabolite/protein peaks from mass spectrometer data

0) "giant.csv" file contains all the chemical compound information of arabidopsis. 

1) Files include "pos" in the files' name is for positive mode, "neg" is for negative mode.

2) "assigned_peak_pos.csv" includes the peaks assigned in the postive mode.
Every row has "x1_pos" indicate the number of the assigned peak, and the peak value.
The rows underneath that peak indicate the row and column number in "giant.csv" file.

3) "output_pos.csv" stores all the matching metabolites in positive mode.  This is raw data of analysis with redundancy.

4) "peak_assign_pos.csv" stores all the matching metabolites without redundancy.

5) "peak_assign_pos1.csv" stores all the matching metabolites with different possible ions.  The first column indicates the row number in "output_pos.csv".

6) The final report of peak assigment for arabidopsis leaves from mass spectrometer is in "report_peak_assignment.pdf".

7) This is program used Sweave to create dynamic reports usig LaTex.  knitr package adds new capability to Sweave and is fully supported by RStudio.

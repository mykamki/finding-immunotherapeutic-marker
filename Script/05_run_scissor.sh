
mkdir $outdir
Rscript run_scissor_rm_165_2.R ~/ANALYSIS/singlecell/Rdata/ ~/ANALYSIS/singlecell/Rdata/integ_3000_mt/dim30/rm_epi/preprocessed_data1.RData null ~/ANALYSIS/singlecell/RESULT/GSE176307/BACI165_2/
[1] "|**************************************************|"
[1] "Performing quality-check for the correlations"
[1] "The five-number summary of correlations:"
        0%        25%        50%        75%       100%
0.05230147 0.22688970 0.31980264 0.35015502 0.58440731
[1] "|**************************************************|"
[1] "Perform cox regression on the given clinical outcomes:"
[1] "alpha = 0.005"
[1] "Scissor identified 263 Scissor+ cells and 730 Scissor- cells."
[1] "The percentage of selected cell is: 65.458%"

[1] "alpha = 0.01"
[1] "Scissor identified 273 Scissor+ cells and 556 Scissor- cells."
[1] "The percentage of selected cell is: 54.647%"

[1] "alpha = 0.05"
[1] "Scissor identified 233 Scissor+ cells and 259 Scissor- cells."
[1] "The percentage of selected cell is: 32.432%"

[1] "alpha = 0.1"
[1] "Scissor identified 174 Scissor+ cells and 163 Scissor- cells."
[1] "The percentage of selected cell is: 22.215%"

[1] "alpha = 0.2"
[1] "Scissor identified 57 Scissor+ cells and 77 Scissor- cells."
[1] "The percentage of selected cell is: 8.833%"
[1] "|**************************************************|"

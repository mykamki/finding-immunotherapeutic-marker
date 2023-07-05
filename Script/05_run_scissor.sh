
mkdir ../Output/SCISSOR/
Rscript run_scissor.R ../Output/BULK/ ../Output/INTEGRATION/IMMUNE/preprocessed_data1.RData ../Output/SCISSOR/

#[1] "|**************************************************|"
#[1] "Performing quality-check for the correlations"
#[1] "The five-number summary of correlations:"
#       0%       25%       50%       75%      100% 
#0.1239261 0.3059452 0.3404038 0.3746918 0.6112234 
#[1] "|**************************************************|"
#[1] "Perform cox regression on the given clinical outcomes:"
#[1] "alpha = 0.005"
#[1] "Scissor identified 361 Scissor+ cells and 665 Scissor- cells."
#[1] "The percentage of selected cell is: 67.633%"

#[1] "alpha = 0.01"
#[1] "Scissor identified 312 Scissor+ cells and 456 Scissor- cells."
#[1] "The percentage of selected cell is: 50.626%"

#[1] "alpha = 0.05"
#[1] "Scissor identified 252 Scissor+ cells and 177 Scissor- cells."
#[1] "The percentage of selected cell is: 28.279%"

#[1] "alpha = 0.1"
#[1] "Scissor identified 173 Scissor+ cells and 115 Scissor- cells."
#[1] "The percentage of selected cell is: 18.985%"
#[1] "|**************************************************|"


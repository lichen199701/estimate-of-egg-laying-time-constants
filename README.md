This script is to create histgrams of C. elegans egg-laying log intervals, along with the estimate time constants based on the observed intervals. 
The model was described in Zhou, Schafer, and Schafer, IEEE% Trans. on Signal Processing, vol. 46, No. 10, October 1998.
This script was revised by Dr. William Schafer from the original script written by Dr. G. Tong Zhou, Georgia Tech.
 
 To use the script:
 1. download the whole folder "estimate-of-egg-laying-time-constants" with all the files
 2. add this folder into the path in MATLAB
 3. open the file "ml_analysis.m."
 4. import your dataset with the observed intervals. 
 5. run [p,lambda1,lambda2]=ml_analysis(x,[name]) -- replace the 'x' with the the name of your imported dataset, replace [name] with the column name of the observed intervals for analysis, in the example dataset, it is 'intervals'. 

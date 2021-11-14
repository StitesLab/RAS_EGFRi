% Running this script will generate the data
% evaluated in Figure 4 and S4

% The matrices of Random factors that are needed for this script
% are each > 100MB and exceed the GitHub file size limits.
% Please contact the corresponding author if you would like files for the exact set of
% random numbers used.
% Alternatively, one could create their own set of random numbers where all
% are within 1, 2, or 3 orders of magnitude, and where the dependent
% parameters that are functions of the independent parameters are also
% within the same order of magnitude.

% The output of the script is ultimately five text files.
% Copies of these five output text files are included.
% They are: Sens_Rand_1e1.txt, Sens_Rand_1e2.txt, Sens_Rand_1e3.txt,
% Sens_Rand2_1e2.txt, and Sens_Rand2_1e3.txt.

load RANDfactor_1e1
RANDfactors=RANDfactor_1e1;
save RANDfactors RANDfactors;
clear all;
MillMuts_CRCMut;
clear all;
Find_Muts_In_Range_UNSTIM;
clear all;
MillMuts_CRCMut_STIM;
clear all;
Find_Muts_In_Range_STIM;
clear all;
DRs_RAND_MUTS;
clear all;
load DRs_RAND_MUTS_output
save DRs_RAND_MUTS_output_1e1;
clear all;
Analyze_RAND_DRs_1e1;
clear all;

load RANDfactor_1e2
RANDfactors=RANDfactor_1e2;
save RANDfactors RANDfactors;
clear all;
MillMuts_CRCMut;
clear all;
Find_Muts_In_Range_UNSTIM;
clear all;
MillMuts_CRCMut_STIM;
clear all;
Find_Muts_In_Range_STIM;
clear all;
DRs_RAND_MUTS;
clear all;
load DRs_RAND_MUTS_output
save DRs_RAND_MUTS_output_1e2;
clear all;
Analyze_RAND_DRs_1e2;
clear all;

load RANDfactor_1e3
RANDfactors=RANDfactor_1e3;
save RANDfactors RANDfactors;
clear all;
MillMuts_CRCMut;
clear all;
Find_Muts_In_Range_UNSTIM;
clear all;
MillMuts_CRCMut_STIM;
clear all;
Find_Muts_In_Range_STIM;
clear all;
DRs_RAND_MUTS;
clear all;
load DRs_RAND_MUTS_output
save DRs_RAND_MUTS_output_1e3;
clear all;
Analyze_RAND_DRs_1e3;
clear all;
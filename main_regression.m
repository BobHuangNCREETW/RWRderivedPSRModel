%-----
%      This is an example regression for the S1 global model using all earthquakes listed in Table S1.
%      The regression's functional form is constructed using the Beta cumulative
%        distribution function.
%      Variables used in this regression analysis include the Rup_Wid_Ratio (RWR) and
%        surface rupture status.
%      Inputs for regression:
%        Rup_Wid_Ratio:   Calculated as the rupture width divided by the fault width.
%        Sur_Rup_Status:  0 indicates no surface rupture, and 1 indicates a rupture
%                           reaching the surface.
%      Outputs for S1-PSR modeling term:
%        alpha_fit_S1     =alpha;
%        var_alpha        =variance of alpha
%        ci_alpha         =confidence interval of alpha
%                                          by Bob J.Y. Huang in September 2026
%-----
load TableS1.mat;
index_surface_rup=find(flag_SurRupYN1Yes0No==1);
index_buried=find(flag_SurRupYN1Yes0No==0);
Sur_Rup_Status(index_surface_rup,1)=1;Sur_Rup_Status(index_buried,1)=0;
% Regression analysis for the S1-PSR model
[alpha_fit_S1,var_alpha,ci_alpha]=sub_fitBetaMLE_fixbeta(Rup_Wid_Ratio,Sur_Rup_Status);

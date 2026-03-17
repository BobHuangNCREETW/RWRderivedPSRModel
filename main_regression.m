%-----
%      This is an example regression for the backbone model (S1) and the proposed 
%        stage 2 (S2) model using all earthquakes listed in Table 1.
%      The regression's functional form is constructed using the Beta cumulative 
%        distribution function.
%      Variables used in this regression analysis include the Rup_Wid_Ratio (RWR) and 
%        surface rupture status. 
%      Inputs for regression: 
%        Rup_Wid_Ratio:   Calculated as the rupture width divided by the fault width.
%        Sur_Rup_Status:  0 indicates no surface rupture, and 1 indicates a rupture 
%                           reaching the surface. 
%        NZbar_Hypo:      The regional or fault-specific normalized mean hypocenter
%                           depth.
%      Outputs for S1-PSR modeling term: 
%        alpha_fit_S1     =alpha;
%        var_alpha        =variance of alpha
%        ci_alpha         =confidence interval of alpha
%      Outputs for S2-PSR modeling term:
%        params_fit_S2(1) =c0; 
%        params_fit_S2(2) =c1
%        cov_mat_S2       =The covariance matrix derived from maximum likelihood 
%                           estimation.
%        ci_for_params    =confidence interval of c0 and c1
%                                          by Bob J.Y. Huang in March 2026
%-----
load Table1.mat;
index_surface_rup=find(flag_SurRupYN1Yes0No==1);
index_buried=find(flag_SurRupYN1Yes0No==0);
Sur_Rup_Status(index_surface_rup,1)=1;Sur_Rup_Status(index_buried,1)=0;
% Regression analysis for the S1-PSR model
[alpha_fit_S1,var_alpha,ci_alpha]=sub_fitBetaMLE_fixbeta(Rup_Wid_Ratio,Sur_Rup_Status);
% Regression analysis for the S2-PSR model
[params_fit_S2,cov_mat_S2,ci_for_params]=sub_fitBetaMLE_fixbeta_addedNMZHypo(Rup_Wid_Ratio,Sur_Rup_Status,N_Zbar_Hypo) 



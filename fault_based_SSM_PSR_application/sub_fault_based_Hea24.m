function [pre_S4_Asp,pre_S3_Area,MBP_RWReq1,pre_Wrup_mean,Sigma_al_Wrup_upper,Sigma_al_Wrup_lower,M_RWReq076]=sub_fault_based_Hea24(Mw,dip,NZbor,Evtlon,Evtlat,FZS_WFlt)
%-----
%      This is the translated fault-based source scaling models (SSM) derived from the original 
%        regional-scale SSM for Huang et al. (SRL, 2024) (denoted as Hea24). The detailed process
%        of the translation is developed and presented in the Appendix of Huang and Abrahamson (in preparation, 2026).
%      Inputs for the SSM:
%      Mw:            Moment magnitude.
%      dip:           Fault dip angle (degrees).
%      NZbor:         The normalized depth of the fault rupture bottom, which uses the seismogenic zone 
%                       bottom depths of 19.3 km, 27.4 km, 20.4 km, and 35.0 km for the CA, TW, JP, and 
%                       other regions, respectively.
%      Evtlon:        Source longitude; Evtlat: Source latitude. (degrees)
%      FZS_WFlt:      The fault-zone-specific fault width.
%      Outputs:
%      pre_S4_Asp:    Mean prediction of the aspect ratio derived from the fault-based Hea24 SSM.
%      pre_S3_Area:   Mean prediction of the area derived from the fault-based Hea24 SSM.    
%      MBP_RWReq1:    The break magnitude at the condition where the rupture-width ratio (RWR) equals 1.
%                       (Eq. A4 in Huang and Abrahamson, 2026)
%      pre_Wrup_mean: Mean prediction of rupture width (W_{Rup}) obtained from the average number of the 
%                       saturated Monte Carlo process.
%      Sigma_al_Wrup: The 1.96 times aleatory variability (\sigma_{al}) from the saturated Monte Carlo 
%                       process, with its upper and lower limits.
%      sigma_Asp:     The \sigma_{al} of the aspect ratio model in Hea24.
%      sigma_A:       The \sigma_{al} of the area model in Hea24.
%                                     by Bob J.Y. Huang in October 2026
%%%%%%
% c0: Constant for M<MBP; c1: slope of M>MBP; 
% MBP: break magnitude for fault rupture reaches the bottom of the seismogenic zone, which is dip-dependent; This is replaced with MBP_RWReq1 
c0=0.1139; c1=0.532;
% c2; c3 : intercept and slope for dip-dependent MBP  
c2=7.17; c3=-0.0105;
% Coefficients c4 and c5 are for Stage 3, representing depth-dependent relations
c4=0.126;c5=-0.196;
% Area scaling
d0=-0.215;d1=-0.068;d2=1.371;
d3_AreaS3_wo_FS=-0.2435;d4_AreaS3_wo_FS=0.3127; % Combined two groups for application used because the user does not know if the fault rupture will be a surface rupture or a buried fault;
% identify region 
minLON_CA=-125;maxLON_CA=-112;minLAT_CA=31;maxLAT_CA=42;
minLON_TW=120;maxLON_TW=123;minLAT_TW=21.7;maxLAT_TW=25.5;
minLON_JP=128.88;maxLON_JP=143;minLAT_JP=29.79;maxLAT_JP=41.85;
if(Evtlon>=minLON_CA&Evtlon<=maxLON_CA&Evtlat>=minLAT_CA&Evtlat<=maxLAT_CA)
  region='CA';region_flag=1; % 1=CA; 2=TW; 3=JP; 4=other
elseif(Evtlon>=minLON_TW&Evtlon<=maxLON_TW&Evtlat>=minLAT_TW&Evtlat<=maxLAT_TW);
  region='TW';region_flag=2;
elseif(Evtlon>=minLON_JP&Evtlon<=maxLON_JP&Evtlat>=minLAT_JP&Evtlat<=maxLAT_JP);
  region='JP';region_flag=3;
else
  region='other';region_flag=4;
end
% Stage 2 model in Hea24: A bilinear model with varying break magnitudes (dip-dependent MBP).
%------ Calculate Asp and A when touching the ZBSZ meaning (RWR=1) that without suffering width-limit effect (Asp^(b_WL) and A^(b_WL)).
c0_FZS_WFlt=4-((d0*(-1)*d2)+d1)-d3_AreaS3_wo_FS+c0+c4;
MBP_RWReq1=2*(log10(FZS_WFlt))-((d4_AreaS3_wo_FS-c5)*NZbor)+c0_FZS_WFlt; 
%------
%------ Replace the original MBP(dip) [MBP_ASPS] with MBP_RWReq1.
if(Mw<=MBP_RWReq1) % Mw smaller than or equal to MBP
  pre_S2_Asp=10.^(c0);
elseif(Mw>MBP_RWReq1) % Mw larger than MBP
  pre_S2_Asp=10.^(c0+c1.*(Mw-MBP_RWReq1));
end
% predicted aspect ratio in S3 model
pre_S3_Asp=10.^(log10(pre_S2_Asp)+c4+(c5.*NZbor));
%-----
pre_S4_Asp=pre_S3_Asp; % turn off RT because we consider Z_BSZ^FZS directly
%%%%%
% Area scaling in Hea24
base_model_logA=Mw-4; % use typical constant 4 as a base model in S1_A(M)
MBP_AreaS2=MBP_RWReq1;
MBPL_AreaS2=MBP_RWReq1+d2;
pre_S2_Area=base_model_logA+max([d0.*(max([Mw,MBP_AreaS2])-MBPL_AreaS2)+d1,d1]);
pre_S2_Area=10.^pre_S2_Area;
pre_S3_Area=10.^(log10(pre_S2_Area)+(d3_AreaS3_wo_FS+d4_AreaS3_wo_FS.*NZbor));
%--- all the above belong to median model, below deal with the sigma model
%================================
%--- Set M dependent aleatory variability for Wrup
sigma_Wrup_global=0.12; % derived from a Monte Carlo sampling of A and Asp using a multivariate log-normal distribution: (sigma_Asp_global=0.1534; sigma_A_global=0.2411; rho=0.3379), which are from residuals between the earthquake data in Table 1 of Hea24 (observation) and the Hea24 model (S3_Asp and S3_A)
small_num_for_sigma=0.02;
M_RWReq076=2*(log10(FZS_WFlt*0.76))-((d4_AreaS3_wo_FS-c5)*NZbor)+c0_FZS_WFlt;
M_RWReq076P1=M_RWReq076+1;
%--- find W at M_{BP}^{(RWR=0,76)+1}
pre_S2_Asp_M_RWReq076P1=10.^(c0+c1.*(M_RWReq076P1-MBP_RWReq1));
pre_S4_Asp_M_RWReq076P1=10.^(log10(pre_S2_Asp_M_RWReq076P1)+c4+(c5.*NZbor));
base_model_logA_M_RWReq076P1=M_RWReq076P1-4; % use typical constant 4 as a base model in S1_A(M)
pre_S2_Area_M_RWReq076P1=base_model_logA_M_RWReq076P1+max([d0.*(max([M_RWReq076P1,MBP_AreaS2])-MBPL_AreaS2)+d1,d1]);
pre_S2_Area_M_RWReq076P1=10.^pre_S2_Area_M_RWReq076P1;
pre_S3_Area_M_RWReq076P1=10.^(log10(pre_S2_Area_M_RWReq076P1)+(d3_AreaS3_wo_FS+d4_AreaS3_wo_FS.*NZbor));
%---
Sigma_taper=min(max(Mw-M_RWReq076,0),1);
Mdep_sigma_Wrup=sigma_Wrup_global-((sigma_Wrup_global-small_num_for_sigma)*Sigma_taper);
%--- Monte Carlo
N=10000;
mu_Wrup_log10=log10(sqrt(pre_S3_Area/pre_S4_Asp));
sigma_Wrup_log10=Mdep_sigma_Wrup;
% convert from log10 to natural log
mu_Wrup_ln=mu_Wrup_log10*log(10);
sigma_Wrup_ln=sigma_Wrup_log10*log(10);
% generate samples
samples=lognrnd(mu_Wrup_ln,sigma_Wrup_ln,N,1);
Wrupsim=samples;
% set width saturation at FZS_WFlt
Wrup_saturated=FZS_WFlt;
Wrupsim_unsaturated=Wrupsim;
index_saturated_Wrup=find(Wrupsim>Wrup_saturated);
Wrupsim(index_saturated_Wrup)=Wrup_saturated;
%---
% Final statistics
pre_Wrup_mean=10^(mean(log10(Wrupsim)));
Sigma_al_Wrup_upper=prctile(Wrupsim,95); % two times sigma of the aleatory variability
Sigma_al_Wrup_lower=prctile(Wrupsim,5);

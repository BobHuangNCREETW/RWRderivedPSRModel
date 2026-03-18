function [pre_S2_PSR,pre_S2_PSR_W_add_Sigma_al,pre_S2_PSR_W_minus_Sigma_al,ASP_MBPRWReq1_NZborMS_woFS,Asp_mean_add_Sigma,Asp_mean_minus_Sigma,W_MBPRWReq1_NZborMS_woFS,W_add_Sigma_al,W_minus_Sigma_al,A_MBPRWReq1_NZborMS_woFS,A_mean_add_Sigma,A_mean_minus_Sigma,RWR_used,NZbor_regional_mean_SSM_used,MBP_RWReq1]=sub_makeComparecurve(Lon,Lat,Mw,dip,ZBSZ);
%                                     by Bob J.Y. Huang in March 2026
%%%% check regional from Hypolon & Hypolat
minLON_CA=-125;maxLON_CA=-112;minLAT_CA=31;maxLAT_CA=42; % California
minLON_TW=119;maxLON_TW=123;minLAT_TW=21.5;maxLAT_TW=25.5; % Taiwan
minLON_JP=128.88;maxLON_JP=145;minLAT_JP=29.79;maxLAT_JP=46; % Japan
minLON_NZ=162;maxLON_NZ=180;minLAT_NZ=-48;maxLAT_NZ=-33; % New Zealand
minLON_TR=26;maxLON_TR=47;minLAT_TR=34;maxLAT_TR=46; % Turkiye
minLON_IR=47;maxLON_IR=64.5;minLAT_IR=24.36;maxLAT_IR=42.83; % Iran
minLON_NWCN=80;maxLON_NWCN=112;minLAT_NWCN=32;maxLAT_NWCN=55; % NorthWest China
minLON_AF=64.5;maxLON_AF=80;minLAT_AF=29.2;maxLAT_AF=46; % Afghanistan
minLON_SEU=10;maxLON_SEU=26;minLAT_SEU=33;maxLAT_SEU=46; % South Europe SEU (combined Greece and Italy)
minLON_SHI=71.26;maxLON_SHI=92;minLAT_SHI=24.36;maxLAT_SHI=37.5; %South Himalayas region
minLON_SA=92;maxLON_SA=109.29;minLAT_SA=17;maxLAT_SA=32; % South Asia
region_nm_forsort=[{'SHI'},{'CA'},{'TW'},{'JP'},{'NZ'},{'TR'},{'IR'},{'NWCN'},{'AF'},{'SEU'},{'SA'}];
gotRegionYN=0;
for j=1:length(region_nm_forsort)
  eval(['minLON_forcheck=minLON_',region_nm_forsort{j},';']);
  eval(['maxLON_forcheck=maxLON_',region_nm_forsort{j},';']);
  eval(['minLAT_forcheck=minLAT_',region_nm_forsort{j},';']);
  eval(['maxLAT_forcheck=maxLAT_',region_nm_forsort{j},';']);
  if(Lon>=minLON_forcheck&Lon<=maxLON_forcheck&Lat>=minLAT_forcheck&Lat<=maxLAT_forcheck)
    eval(['Region_assign=''',region_nm_forsort{j},''';']);
    gotRegionYN=1;
    break
  end
end
if(gotRegionYN==0)
  Region_assign='others' % means other regions
end
if(strcmp(Region_assign,'CA')==1)
  NMZHypo_used=0.5122;
  NZbor_regional_mean_SSM_used=0.7433; % mean regional NZbor from Table1 in Huang et al. (2024); CA region
elseif(strcmp(Region_assign,'SHI')==1)
  NMZHypo_used=0.4055;
  NZbor_regional_mean_SSM_used=0.5734; % mean regional NZbor from Table1 in Huang et al. (2024); other region 
elseif(strcmp(Region_assign,'TW')==1)
  NMZHypo_used=0.5198;
  NZbor_regional_mean_SSM_used=0.7138; % mean regional NZbor from Table1 in Huang et al. (2024); TW region
elseif(strcmp(Region_assign,'JP')==1)
  NMZHypo_used=0.5036;
  NZbor_regional_mean_SSM_used=0.7038; % mean regional NZbor from Table1 in Huang et al. (2024); JP region
elseif(strcmp(Region_assign,'NZ')==1)
  NMZHypo_used=0.4695;
  NZbor_regional_mean_SSM_used=0.5734; % mean regional NZbor from Table1 in Huang et al. (2024); other region
elseif(strcmp(Region_assign,'TR')==1)
  NMZHypo_used=0.4842;
  NZbor_regional_mean_SSM_used=0.5734; % mean regional NZbor from Table1 in Huang et al. (2024); other region
elseif(strcmp(Region_assign,'IR')==1)
  NMZHypo_used=0.5483;
  NZbor_regional_mean_SSM_used=0.5734; % mean regional NZbor from Table1 in Huang et al. (2024); other region
elseif(strcmp(Region_assign,'NWCN')==1)
  NMZHypo_used=0.3504;
  NZbor_regional_mean_SSM_used=0.5734; % mean regional NZbor from Table1 in Huang et al. (2024); other region
elseif(strcmp(Region_assign,'AF')==1)
  NMZHypo_used=0.3760;
  NZbor_regional_mean_SSM_used=0.5734; % mean regional NZbor from Table1 in Huang et al. (2024); other region
elseif(strcmp(Region_assign,'SEU')==1)
  NMZHypo_used=0.6116;
  NZbor_regional_mean_SSM_used=0.5734; % mean regional NZbor from Table1 in Huang et al. (2024); other region
elseif(strcmp(Region_assign,'SA')==1)
  NMZHypo_used=0.3657;
  NZbor_regional_mean_SSM_used=0.5734; % mean regional NZbor from Table1 in Huang et al. (2024); other region
end
%%%%
%---- Use regional mean scaling (MS) NZbor directly, w/o free surface effect
%     Applied M_{BP}^{RWR=1} instead of M_{BP}, which correlated with both fault zone specific 
%       W_{Flt} (FZS_WFlt) and dip angle.
FZS_WFlt=ZBSZ/sind(dip); 
[ASP_MBPRWReq1_NZborMS_woFS,A_MBPRWReq1_NZborMS_woFS,MBP_RWReq1,W_MBPRWReq1_NZborMS_woFS,sigma_al_Wrup_upper,sigma_al_Wrup_lower,sigma_Asp,sigma_A]=sub_fault_based_Hea24(Mw,dip,NZbor_regional_mean_SSM_used,Lon,Lat,FZS_WFlt);
Asp_mean_add_Sigma=10.^(log10(ASP_MBPRWReq1_NZborMS_woFS)+(2.*sigma_Asp));
Asp_mean_minus_Sigma=10.^(log10(ASP_MBPRWReq1_NZborMS_woFS)-(2.*sigma_Asp));
A_mean_add_Sigma=10.^(log10(A_MBPRWReq1_NZborMS_woFS)+(2.*sigma_A));
A_mean_minus_Sigma=10.^(log10(A_MBPRWReq1_NZborMS_woFS)-(2.*sigma_A));
RWR_used=W_MBPRWReq1_NZborMS_woFS./(FZS_WFlt.*1.1);
if(RWR_used>1) RWR_used=1; end
%----
[pre_S2_PSR,y_lo_PSR,y_high_PSR]=sub_RWR_derived_PSRmodel(RWR_used,NMZHypo_used);
% considering 1.96 times of aleatory variability (5% to 95% confident interval) of Wrup
W_add_Sigma_al=sigma_al_Wrup_upper;
RWR_used_add_Sigma_al=W_add_Sigma_al./(FZS_WFlt.*1.1); % added 10% for conservative consideration of the Wrup can existed a bit of WFlt
if(RWR_used_add_Sigma_al>1) RWR_used_add_Sigma_al=1; end
[pre_S2_PSR_W_add_Sigma_al,y_lo_PSR_W_add_Sigma_al,y_high_PSR_W_add_Sigma_al]=sub_RWR_derived_PSRmodel(RWR_used_add_Sigma_al,NMZHypo_used);
W_minus_Sigma_al=sigma_al_Wrup_lower;
RWR_used_minus_Sigma_al=W_minus_Sigma_al./(FZS_WFlt.*1.1);
if(RWR_used_minus_Sigma_al>1) RWR_used_minus_Sigma_al=1; end
[pre_S2_PSR_W_minus_Sigma_al,y_lo_PSR_W_minus_Sigma_al,y_high_PSR_W_minus_Sigma_al]=sub_RWR_derived_PSRmodel(RWR_used_minus_Sigma_al,NMZHypo_used);





%-----
%      This program demonstrates an example of translating the regional-scale-derived source scaling models (SSM)
%        of Huang et al. (SRL, 2024) (denoted as Hea24) into a fault-based SSM, which is combined with the Beta 
%        cumulative distribution function (Beta-CDF) derived stage 2 probability model of surface rupture (S2-PSR)
%        model. The translation procedure and the S2-PSR are developed in Huang and Abrahamson (in preparation, 2026). 
%      The example scenarios are demonstrated based on various bottom depths of the seismogenic zone (Z_BSZ), or 
%        equivalently, on various fault-zone-specific fault widths (W_{Flt}^{FZS}) in a fixed fault dip angle (\nu) 
%        condition.  
%                                                               by Bob J.Y. Huang in March 2026 
%-----
Lon_CA=-122.235;Lat_CA=37.858; % for identification only, no use because the regional term (RT(M)) is turning off in Hea24
M_CA_used=4.8:0.1:8;
dip_CA=80; % an example fixed fault dip angle
ZBSZ_CA_Comp1=4.924;ZBSZ_CA_Comp2=9.848;ZBSZ_CA_Comp3=19.6962; % equivelent to W_{Flt}^{FZS} equal to 5, 10, and 20 km
%--- for legend labels and variables in plot
W_flt_Comp1=ZBSZ_CA_Comp1/sind(dip_CA);W_flt_Comp2=ZBSZ_CA_Comp2/sind(dip_CA);
W_flt_Comp3=ZBSZ_CA_Comp3/sind(dip_CA);
legend_in1=['Z_B_S_Z=',num2str(round(ZBSZ_CA_Comp1*10)/10),'km; Hea24 (CA, w/o FS,',sprintf('\n'),'        W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp1*10)/10),'km, '];
legend_in2=['Z_B_S_Z=',num2str(round(ZBSZ_CA_Comp2*10)/10),'km; Hea24 (CA, w/o FS,',sprintf('\n'),'        W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp2*10)/10),'km, '];
legend_in3=['Z_B_S_Z=',num2str(round(ZBSZ_CA_Comp3*10)/10),'km; Hea24 (CA, w/o FS,',sprintf('\n'),'        W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp3*10)/10),'km, '];
legendshiftx=0.57;
%--- Calculate and plot the translated progress for Hea24 combined with the S2-PSR model (Fig. 13).
title_in=['\nu=',num2str(dip_CA),'^\circ'];
outnm='Progress_CA_ZBSZ_scaling';
[pre_S2_PSR1_CA_FZS_WFlt5,pre_S2_PSR1_CA_FZS_WFlt5_Psigmaal,pre_S2_PSR1_CA_FZS_WFlt5_Msigmaal,pre_S2_PSR1_CA_FZS_WFlt10,pre_S2_PSR1_CA_FZS_WFlt10_Psigmaal,pre_S2_PSR1_CA_FZS_WFlt10_Msigmaal,pre_S2_PSR1_CA_FZS_WFlt20,pre_S2_PSR1_CA_FZS_WFlt20_Psigmaal,pre_S2_PSR1_CA_FZS_WFlt20_Msigmaal,RWR_used1_CA_FZS_WFlt5,RWR_used2_CA_FZS_WFlt10,RWR_used3_CA_FZS_WFlt20]=sub_pt_progress_comp_threecurves(Lon_CA,Lat_CA,M_CA_used,dip_CA,ZBSZ_CA_Comp1,dip_CA,ZBSZ_CA_Comp2,dip_CA,ZBSZ_CA_Comp3,legend_in1,legend_in2,legend_in3,title_in,legendshiftx,outnm);
%--- Calculate CPSR: Mea24 from Mammarella et al. EQS (2024), doi: 10.1177/87552930241293570, for comparison.
if(ZBSZ_CA_Comp1>=10)
  ZBSZ_CA_Comp1_star=round(((ZBSZ_CA_Comp1-9.5193)/0.64431)*10)/10; % from Fig. S5 (in the electronic supplement) because the ZBSZ used in Mea24 is from the distribution of ZHypo, so the transfer equation is required (Eq. 20)
else
  ZBSZ_CA_Comp1_star=round((ZBSZ_CA_Comp1*(2/3))*10)/10; % considering the 2/3 relation to avoid negative ZHypo-derived ZBSZ^* (Eq. 21)
end
if(ZBSZ_CA_Comp2>=10)
  ZBSZ_CA_Comp2_star=round(((ZBSZ_CA_Comp2-9.5193)/0.64431)*10)/10;
else
  ZBSZ_CA_Comp2_star=round((ZBSZ_CA_Comp2*(2/3))*10)/10;
end
if(ZBSZ_CA_Comp3>=10)
  ZBSZ_CA_Comp3_star=round(((ZBSZ_CA_Comp3-9.5193)/0.64431)*10)/10;
else
  ZBSZ_CA_Comp3_star=round((ZBSZ_CA_Comp3*(2/3))*10)/10;
end
%---- Legend text for CPSR: Mea24
legend_inCPSR1=['CPSR: Mea24 (T17, CA, SS, \nu=',num2str(dip_CA),'^\circ, Z_B_S_Z*=',num2str(ZBSZ_CA_Comp1_star),'km)'];
legend_inCPSR2=['CPSR: Mea24 (T17, CA, SS, \nu=',num2str(dip_CA),'^\circ, Z_B_S_Z*=',num2str(ZBSZ_CA_Comp2_star),'km)'];
legend_inCPSR3=['CPSR: Mea24 (T17, CA, SS, \nu=',num2str(dip_CA),'^\circ, Z_B_S_Z*=',num2str(ZBSZ_CA_Comp3_star),'km)'];
%----
outfignm='Mea24_case1';
[Mw_Mea24_case1,P_Mea24_case1]=sub_CPSR_Mea24(2,5,'CA_S',80,2,2.5,ZBSZ_CA_Comp1_star,2,1,outfignm);
outfignm='Mea24_case2';
[Mw_Mea24_case2,P_Mea24_case2]=sub_CPSR_Mea24(2,5,'CA_S',80,2,2.5,ZBSZ_CA_Comp2_star,2,1,outfignm);
outfignm='Mea24_case3';
[Mw_Mea24_case3,P_Mea24_case3]=sub_CPSR_Mea24(2,5,'CA_S',80,2,2.5,ZBSZ_CA_Comp3_star,2,1,outfignm);
% plot for comparison (Fig. 14 (a))
legend_PSR_in1=['S2-PSR; Z_B_S_Z=',num2str(round(ZBSZ_CA_Comp1*10)/10),'km; W_R_u_p: Hea24 (CA, \nu=',num2str(dip_CA),'^\circ, W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp1*10)/10),'km)'];
legend_PSR_in2=['S2-PSR; Z_B_S_Z=',num2str(round(ZBSZ_CA_Comp2*10)/10),'km; W_R_u_p: Hea24 (CA, \nu=',num2str(dip_CA),'^\circ, W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp2*10)/10),'km)'];
legend_PSR_in3=['S2-PSR; Z_B_S_Z=',num2str(round(ZBSZ_CA_Comp3*10)/10),'km; W_R_u_p: Hea24 (CA, \nu=',num2str(dip_CA),'^\circ, W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp3*10)/10),'km)'];
outnm='Comp_CA_ZBSZ_scaling';Comp_SoF='SS';
[nu]=sub_pt_RWRPSR_w_other_models(pre_S2_PSR1_CA_FZS_WFlt5,pre_S2_PSR1_CA_FZS_WFlt5_Psigmaal,pre_S2_PSR1_CA_FZS_WFlt5_Msigmaal,pre_S2_PSR1_CA_FZS_WFlt10,pre_S2_PSR1_CA_FZS_WFlt10_Psigmaal,pre_S2_PSR1_CA_FZS_WFlt10_Msigmaal,pre_S2_PSR1_CA_FZS_WFlt20,pre_S2_PSR1_CA_FZS_WFlt20_Psigmaal,pre_S2_PSR1_CA_FZS_WFlt20_Msigmaal,M_CA_used,M_CA_used,M_CA_used,P_Mea24_case1,P_Mea24_case2,P_Mea24_case3,Mw_Mea24_case1,Mw_Mea24_case2,Mw_Mea24_case2,legend_PSR_in1,legend_PSR_in2,legend_PSR_in3,legend_inCPSR1,legend_inCPSR2,legend_inCPSR3,outnm,Comp_SoF);









%-----
%      This program demonstrates an example of translating the regional-scale-derived source scaling models (SSM)
%        of Huang et al. (SRL, 2024) (denoted as Hea24) into a fault-based SSM, which is
%        combined with the Beta cumulative distribution function (Beta-CDF) derived SoF 
%        probability model of surface rupture (using a global SS type as an example).
%      The translation procedure and the global SoF PSR models are developed in 
%        Huang and Abrahamson (in preparation, 2026).
%      The example scenarios are demonstrated based on various bottom depths of the 
%        seismogenic zone (Z_BSZ), or equivalently, on various fault-zone-specific 
%        fault widths (W_{Flt}^{FZS}) in a fixed fault dip angle (\nu) condition.
%                                                    by Bob J.Y. Huang in September 2026
%-----
Lon_global=87.4;Lat_global=52.8; % A random place for specify the other region used in Hea24
M_global_used=4.8:0.1:8;
dip_global=80; 
ZBSZ_global_Comp1=9.8481;ZBSZ_global_Comp2=14.7721;ZBSZ_global_Comp3=19.6962; % equivelent to W_{Flt}^{FZS} equal to 10, 15, and 20 km
Rake_global_used=0; % set SS
%--- for legend labels and variables in plot
W_flt_Comp1=ZBSZ_global_Comp1/sind(dip_global);W_flt_Comp2=ZBSZ_global_Comp2/sind(dip_global);
W_flt_Comp3=ZBSZ_global_Comp3/sind(dip_global);
legend_in1=['Hea24 (other, W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp1*10)/10),'km, ',sprintf('\n'),'  w/o FS, '];
legend_in2=['Hea24 (other, W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp2*10)/10),'km, ',sprintf('\n'),'  w/o FS, '];
legend_in3=['Hea24 (other, W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp3*10)/10),'km, ',sprintf('\n'),'  w/o FS, '];
%--- Calculate and plot the translated progress for Hea24 combined with the global SoF PSR model (Fig. A2 (b) in Appendix).
title_in=['\nu=',num2str(dip_global),'^\circ'];
outnm='Progress_global_ZBSZ_scaling_SS';legendshiftx=0.57;
[pre_SoF_PSR1_global_FZS_WFlt5,pre_SoF_PSR1_global_FZS_WFlt5_Psigmaal,pre_SoF_PSR1_global_FZS_WFlt5_Msigmaal,pre_SoF_PSR1_global_FZS_WFlt10,pre_SoF_PSR1_global_FZS_WFlt10_Psigmaal,pre_SoF_PSR1_global_FZS_WFlt10_Msigmaal,pre_SoF_PSR1_global_FZS_WFlt20,pre_SoF_PSR1_global_FZS_WFlt20_Psigmaal,pre_SoF_PSR1_global_FZS_WFlt20_Msigmaal,RWR_used1_global_FZS_WFlt5,RWR_used2_global_FZS_WFlt10,RWR_used3_global_FZS_WFlt20]=sub_pt_progress_comp_threecurves(Lon_global,Lat_global,Rake_global_used,M_global_used,dip_global,ZBSZ_global_Comp1,dip_global,ZBSZ_global_Comp2,dip_global,ZBSZ_global_Comp3,legend_in1,legend_in2,legend_in3,title_in,legendshiftx,outnm);
%--- Calculate Mea25 from Mammarella et al. EQS (2025), doi: 10.1177/87552930241293570, for comparison.
%---- Legend text for Mea25
legend_inMea25_1=['Mea25, AGG, SS, \nu=',num2str(dip_global),'^\circ, W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp1*10)/10),'km'];
legend_inMea25_2=['Mea25, AGG, SS, \nu=',num2str(dip_global),'^\circ, W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp2*10)/10),'km'];
legend_inMea25_3=['Mea25, AGG, SS, \nu=',num2str(dip_global),'^\circ, W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp3*10)/10),'km'];
%----
outfignm='Mea25_case1';
[Mw_Mea25_case1,P_Mea25_case1]=sub_Mea25(2,5,'AGG_S',80,2,2.5,ZBSZ_global_Comp1,2,1,outfignm);
outfignm='Mea25_case2';
[Mw_Mea25_case2,P_Mea25_case2]=sub_Mea25(2,5,'AGG_S',80,2,2.5,ZBSZ_global_Comp2,2,1,outfignm);
outfignm='Mea25_case3';
[Mw_Mea25_case3,P_Mea25_case3]=sub_Mea25(2,5,'AGG_S',80,2,2.5,ZBSZ_global_Comp3,2,1,outfignm);
% plot for comparison (Fig. 12 (a))
legend_PSR_in1=['This study, SS; W_R_u_p: Hea24 (W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp1*10)/10),'km)'];
legend_PSR_in2=['This study, SS; W_R_u_p: Hea24 (W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp2*10)/10),'km)'];
legend_PSR_in3=['This study, SS; W_R_u_p: Hea24 (W_F_l_t^F^Z^S=',num2str(round(W_flt_Comp3*10)/10),'km)'];
outnm='Comp_global_ZBSZ_scaling_SS';Comp_SoF='SS';
[nu]=sub_pt_RWRPSR_w_other_models(pre_SoF_PSR1_global_FZS_WFlt5,pre_SoF_PSR1_global_FZS_WFlt5_Psigmaal,pre_SoF_PSR1_global_FZS_WFlt5_Msigmaal,pre_SoF_PSR1_global_FZS_WFlt10,pre_SoF_PSR1_global_FZS_WFlt10_Psigmaal,pre_SoF_PSR1_global_FZS_WFlt10_Msigmaal,pre_SoF_PSR1_global_FZS_WFlt20,pre_SoF_PSR1_global_FZS_WFlt20_Psigmaal,pre_SoF_PSR1_global_FZS_WFlt20_Msigmaal,M_global_used,M_global_used,M_global_used,P_Mea25_case1,P_Mea25_case2,P_Mea25_case3,Mw_Mea25_case1,Mw_Mea25_case2,Mw_Mea25_case2,legend_PSR_in1,legend_PSR_in2,legend_PSR_in3,legend_inMea25_1,legend_inMea25_2,legend_inMea25_3,outnm,Comp_SoF);





function [pre_S4_Asp,pre_S3_Area,MBP_RWReq1,pre_Wrup_mean,Sigma_al_Wrup_upper,Sigma_al_Wrup_lower,sigma_Asp,sigma_A]=sub_fault_based_Hea24(Mw,dip,NZbor,Evtlon,Evtlat,FZS_WFlt)
%-----
%      This is the translated fault-based source scaling models (SSM) derived from the original 
%        regional-scale SSM for Huang et al. (SRL, 2024) (denoted as Hea24). The detailed process
%        of the translation is developed and presented in Huang and Abrahamson (in preparation, 2026).
%      Inputs for the SSM:
%      Mw:            Moment magnitude.
%      dip:           Fault dip angle (degree).
%      NZbor:         The normalized depth of the fault rupture bottom, which uses the seismogenic zone 
%                       bottom depths of 19.3 km, 27.4 km, 20.4 km, and 35.0 km for the CA, TW, JP, and 
%                       other regions, respectively.
%      Evtlon:        Source longitude; Evtlat: Source latitude. (degree)
%      FZS_WFlt:      The fault-zone-specific fault width.
%      Outputs:
%      pre_S4_Asp:    Mean prediction of the aspect ratio derived from the fault-based Hea24 SSM.
%      pre_S3_Area:   Mean prediction of the area derived from the fault-based Hea24 SSM.    
%      MBP_RWReq1:    The break magnitude at the condition where the rupture-width ratio (RWR) equals 1.
%                       (Eq. 15 in Huang and Abrahamson, 2026)
%      pre_Wrup_mean: Mean prediction of rupture width (W_{Rup}) obtained from the average number of the 
%                       saturated Monte Carlo process.
%      Sigma_al_Wrup: The 1.96 times aleatory variability (\sigma_{al}) from the saturated Monte Carlo 
%                       process, with its upper and lower limits.
%      sigma_Asp:     The \sigma_{al} of the aspect ratio model in Hea24.
%      sigma_A:       The \sigma_{al} of the area model in Hea24.
%                                     by Bob J.Y. Huang in March 2026
%%%%%%
% c0: Constant for M<MBP; c1: slope of M>MBP; 
% MBP: break magnitude for fault rupture reaches the bottom of the seismogenic zone, which is dip dependent; This is replaced with MBP_RWReq1 
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
MBP_RWReq1=2*(log10(FZS_WFlt*0.9))-((d4_AreaS3_wo_FS-c5)*NZbor)+c0_FZS_WFlt; % Applied 90% of FZS_WFlt because the rupture energy will continue rupturing down-dip slightly below the width limit when touching BSZ, as seen in empirical data. We set a conservative value for 90% of FZS_BSZ as an initial point that contacts the BSZ.
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
%---- Calculate the aleatory variability from the sigma of residual Asp, sigma of residual A, and residual correlation coefficient (rho) using the delta method. 
sigma_Asp=0.1534;sigma_A=0.2411;rho=0.3379; % from residuals between earthquake data in Table 1 of Hea24 (observation) and the Hea24 model (S3_Asp and S3_A)
%----Multivariate monte Carlo method for estimating aleatory variability
% Residual covariance matrix
Res_Cov_Mat=[sigma_A^2, rho*sigma_A*sigma_Asp;
             rho*sigma_Asp*sigma_A, sigma_Asp^2];
% Monte Carlo
N=10000;
% set saturation for W_Rup (upper bound for RWR=1)
Wrup_saturated=FZS_WFlt*1.1; % Rupture width saturates at the fault-zone-specific fault width, with the rupture-width ratio equal to 1.
mu_Asp=pre_S4_Asp;mu_A=pre_S3_Area;
mu=[mu_A mu_Asp];
% Monte Carlo sampling for normal distribution in log10 units 
samples=mvnrnd(log10(mu),Res_Cov_Mat,N);
Asim=10.^(samples(:,1));
Aspsim=10.^(samples(:,2));
% Calculate the rupture width of the simulated samples.
Wrupsim=sqrt(Asim./Aspsim);
Wrupsim_unsaturated=Wrupsim;
index_saturated_Wrup=find(Wrupsim>Wrup_saturated); 
Wrupsim(index_saturated_Wrup)=Wrup_saturated; % set W_Rup saturation
% Final statistics
pre_Wrup_mean=10^(mean(log10(Wrupsim)));  
Sigma_al_Wrup_upper=prctile(Wrupsim,95); % two times aleatory variability
Sigma_al_Wrup_lower=prctile(Wrupsim,5);
%----
index_Wrup_gt_110per_FZS_WFlt=find(pre_Wrup_mean>FZS_WFlt*1.1); % Determine the specific fault width of the fault zone, then set a limit where the fault rupture can still go slightly deeper by 10%, indicating that a larger rupture energy can still break through a bit at a higher geothermal or surrounding pressure environment (within the brittle-ductile transition zone).
pre_Wrup_mean(index_Wrup_gt_110per_FZS_WFlt)=FZS_WFlt*1.1;
pre_S3_Area(index_Wrup_gt_110per_FZS_WFlt)=pre_S4_Asp(index_Wrup_gt_110per_FZS_WFlt).*(pre_Wrup_mean(index_Wrup_gt_110per_FZS_WFlt).^2);
%----
%plotsliceYN='Y';
plotsliceYN='N';
if(plotsliceYN=='Y')
figure(100);
H100=subplot(2,3,1);
  semilogy(Mw,mu_Asp,'sr ','linewidth',2,'markersize',8);hold on;
  Mw_array(1:length(Aspsim))=Mw;
  semilogy(Mw,Aspsim,'ok ','linewidth',2,'markersize',6);
  semilogy(Mw,mu_Asp,'sr ','linewidth',2,'markersize',8);
  xlabel('M','fontsize',13,'fontweight','bold');
  ylabel('Asp','fontsize',13,'fontweight','bold');
  ylim([10^-1,10^2]);
  grid on;
  if(Mw<=7)
    xlim([Mw-1 Mw+1]);
  else
    xlim([Mw-1 Mw+0.5]);
  end
H101=subplot(2,3,4);
  [Asphist_y,Asphist_x]=hist(log10(Aspsim));
  Lbar_Asp=bar(Asphist_x,Asphist_y,'facecolor','k','barwidth',1);hold on;
  xlabel('Asp (log_1_0 unit)','fontsize',13,'fontweight','bold');
  ylabel('count','fontsize',13,'fontweight','bold');
  grid on;
  if(Mw==7&FZS_WFlt==20)
    xlim([0 1]);
  end
H102=subplot(2,3,2);
  semilogy(Mw,mu_A,'sr ','linewidth',2,'markersize',8);hold on;  
  semilogy(Mw,Asim,'ok ','linewidth',2,'markersize',6);
  semilogy(Mw,mu_A,'sr ','linewidth',2,'markersize',8);
  xlabel('M','fontsize',13,'fontweight','bold');
  ylabel('A (km^2)','fontsize',13,'fontweight','bold');
  if(Mw<=7)
    xlim([Mw-1 Mw+1]);
  else
    xlim([Mw-1 Mw+0.5]);
  end
  grid on;
H103=subplot(2,3,5);
  [Ahist_y,Ahist_x]=hist(log10(Asim));
  Lbar_A=bar(Ahist_x,Ahist_y,'facecolor','k','barwidth',1);hold on;
  xlabel('A (km^2, log_1_0 unit)','fontsize',13,'fontweight','bold');
  ylabel('count','fontsize',13,'fontweight','bold');
  grid on;
  if(Mw==7&FZS_WFlt==20)
    xlim([2 4]);
  end
H104=subplot(2,3,3);
  semilogy(Mw,pre_Wrup_mean,'sr ','linewidth',2,'markersize',8);hold on;
  semilogy(Mw,Wrupsim(1),'ok ','linewidth',2,'markersize',6);
  semilogy(Mw+0.15,mean(Wrupsim_unsaturated),'s','color',[0.3,0.3,0.3],'linewidth',2,'markersize',8);
  semilogy(Mw+0.15,Wrupsim_unsaturated(1),'o','color',[0.7,0.7,0.7],'linewidth',2,'markersize',6); 
  lH=legend('\mu (saturated)','samples (saturated)','\mu (unsaturated)','samples (unsaturated)',2);
  if(Mw>7)
    po=get(lH,'position');set(lH,'position',[po(1)+0.017,po(2)+0.12,po(3),po(4)],'fontsize',10,'fontweight','bold');
  else
    po=get(lH,'position');set(lH,'position',[po(1)+0.01,po(2)+0.12,po(3),po(4)],'fontsize',10,'fontweight','bold');
  end
  semilogy(Mw,pre_Wrup_mean,'sr ','linewidth',2,'markersize',8);hold on;
  semilogy(Mw,Wrupsim,'ok ','linewidth',2,'markersize',6);
  semilogy(Mw+0.15,mean(Wrupsim_unsaturated),'s','color',[0.3,0.3,0.3],'linewidth',2,'markersize',8);
  semilogy(Mw+0.15,Wrupsim_unsaturated,'o','color',[0.7,0.7,0.7],'linewidth',2,'markersize',6);
  semilogy(Mw+0.15,mean(Wrupsim_unsaturated),'s','color',[0.3,0.3,0.3],'linewidth',2,'markersize',8);
  semilogy(Mw,pre_Wrup_mean,'sr ','linewidth',2,'markersize',8);
  if(Mw<7)
    ylim([10^0 10^2]);
  else
    ylim([0.5 500]);
  end
  if(Mw<=7)
    xlim([Mw-1 Mw+1]);
  else
    xlim([Mw-1 Mw+0.5]);
  end
  grid on;
  xlabel('M','fontsize',13,'fontweight','bold');
  ylabel('W_R_u_p (km)','fontsize',13,'fontweight','bold');
  set(H104,'xtick',[4:0.5:8]);
H105=subplot(2,3,6);
  hold on;
  [Nhist_y1,Nhist_x1]=hist((Wrupsim));
  [Nhist_y2,Nhist_x2]=hist((Wrupsim_unsaturated));
  Lbar2=bar(Nhist_x2,Nhist_y2,'facecolor',[0.7,0.7,0.7],'barwidth',1);
  Lbar1=bar(Nhist_x1,Nhist_y1,'k','barwidth',1);
  xlabel('W_R_u_p (km)','fontsize',13,'fontweight','bold');
  ylabel('count','fontsize',13,'fontweight','bold');
  Mw_print=num2str(Mw);index_dot=find(Mw_print=='.');Mw_print(index_dot)='_';
  max_Wrupsim=max(max([Wrupsim;Wrupsim_unsaturated]));
  max_x_H105=ceil(max_Wrupsim/5)*5;
  if(max_x_H105>50)
    max_x_H105=50;
  end
  xlim([0 max_x_H105]);
  grid on;
fontsize_used=14;
set(H100,'linewidth',2,'fontsize',fontsize_used,'fontweight','bold','xtick',[4:0.5:8.5]);
set(H101,'linewidth',2,'fontsize',fontsize_used,'fontweight','bold');
set(H102,'linewidth',2,'fontsize',fontsize_used,'fontweight','bold','xtick',[4:0.5:8.5]);
set(H103,'linewidth',2,'fontsize',fontsize_used,'fontweight','bold');
set(H104,'linewidth',2,'fontsize',fontsize_used,'fontweight','bold','xtick',[4:0.5:8.5],'ytick',[10^-1,10^0,10^1,10^2]);
if(max_x_H105<=30)
  set(H105,'linewidth',2,'fontsize',fontsize_used,'fontweight','bold','xtick',[0:5:100]);
else
  set(H105,'linewidth',2,'fontsize',fontsize_used,'fontweight','bold','xtick',[0:10:100]);
end
  outnm=['Distribution_slice_Mw',Mw_print,'_FZS_WFlt_',num2str(FZS_WFlt),'_saturated'];
  print(gcf,outnm,'-deps2c');unix(['convert -density 200 ',outnm,'.eps ',outnm,'.png']);
close(100);
end








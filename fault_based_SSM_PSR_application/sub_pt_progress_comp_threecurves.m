function [pre_S2_PSR1,pre_S2_PSR_add2sigmaal1,pre_S2_PSR_minus2sigmaal1,pre_S2_PSR2,pre_S2_PSR_add2sigmaal2,pre_S2_PSR_minus2sigmaal2,pre_S2_PSR3,pre_S2_PSR_add2sigmaal3,pre_S2_PSR_minus2sigmaal3,RWR_used1,RWR_used2,RWR_used3]=sub_pt_progress_comp_threecurves(Lon,Lat,Mw,dip1,ZBSZ1,dip2,ZBSZ2,dip3,ZBSZ3,legend_in1,legend_in2,legend_in3,title_in,legendshiftx,outnm);
%                                     by Bob J.Y. Huang in March 2026
%%% evaluate regions
index_dash=find(outnm=='_');
region=outnm(index_dash(1)+1:index_dash(2)-1);
Comptype=outnm(index_dash(2)+1:index_dash(3)-1);
%%%
%--- Additional plot with small and large M for all scaling.
if(strcmp(region,'CA')==1&strcmp(Comptype,'dip')==1)
  tar_LargeM=6.6;tar_SmallM=6;
elseif(strcmp(region,'CA')==1&strcmp(Comptype,'ZBSZ')==1)
  tar_SmallM=5.6;
  tar_LargeM=0; % for run, did not plot for this situation
elseif(strcmp(region,'TW')==1)%&strcmp(Comptype,'ZBSZ')==1)
  tar_LargeM=6.65;tar_SmallM=6;
elseif(strcmp(region,'SHI')==1)
  tar_LargeM=6.7;tar_SmallM=6;
else
  tar_LargeM=6.8;tar_SmallM=6.2;
end
%---
fid_out=fopen('Para_used.csv','w');
fprintf(fid_out,'%s\n','Mw,Asp_case1,WRup_case1,RWR_case1,MBP_FZS_WFlt_case1,Asp_case2,WRup_case2,RWR_case2,MBP_FZS_WFlt_case2,Asp_case3,WRup_case3,RWR_case3,MBP_FZS_WFlt_case3');
W_flt_FZS1=ZBSZ1/sind(dip1);W_flt_FZS2=ZBSZ2/sind(dip2);W_flt_FZS3=ZBSZ3/sind(dip3);
for i=1:length(Mw)
  [pre_S2_PSR1(i),pre_S2_PSR_add2sigmaal1(i),pre_S2_PSR_minus2sigmaal1(i),ASP_used1(i),ASP_addsigma1(i),ASP_minussigma1(i),W_used1(i),W_add2sigmaal1(i),W_minus2sigmaal1(i),A_used1(i),A_addsigma1(i),A_minussigma1(i),RWR_used1(i),NZbor_regional_mean_SSM_used1,MBP_RWReq1_1]=sub_makeComparecurve(Lon,Lat,Mw(i),dip1,ZBSZ1);
  [pre_S2_PSR2(i),pre_S2_PSR_add2sigmaal2(i),pre_S2_PSR_minus2sigmaal2(i),ASP_used2(i),ASP_addsigma2(i),ASP_minussigma2(i),W_used2(i),W_add2sigmaal2(i),W_minus2sigmaal2(i),A_used2(i),A_addsigma2(i),A_minussigma2(i),RWR_used2(i),NZbor_regional_mean_SSM_used2,MBP_RWReq1_2]=sub_makeComparecurve(Lon,Lat,Mw(i),dip2,ZBSZ2);
  [pre_S2_PSR3(i),pre_S2_PSR_add2sigmaal3(i),pre_S2_PSR_minus2sigmaal3(i),ASP_used3(i),ASP_addsigma3(i),ASP_minussigma3(i),W_used3(i),W_add2sigmaal3(i),W_minus2sigmaal3(i),A_used3(i),A_addsigma3(i),A_minussigma3(i),RWR_used3(i),NZbor_regional_mean_SSM_used3,MBP_RWReq1_3]=sub_makeComparecurve(Lon,Lat,Mw(i),dip3,ZBSZ3);
  fprintf(fid_out,'%4.2f,%6.2f,%6.2f,%5.3f,%4.2f,%6.2f,%6.2f,%5.3f,%4.2f,%6.2f,%6.2f,%5.3f,%4.2f\n',Mw(i),ASP_used1(i),W_used1(i),RWR_used1(i),MBP_RWReq1_1,ASP_used2(i),W_used2(i),RWR_used2(i),MBP_RWReq1_2,ASP_used3(i),W_used3(i),RWR_used3(i),MBP_RWReq1_3);
end
fclose(fid_out); 
index_M55=find(abs(Mw-5.5)==min(abs(Mw-5.5)));
pre_S2_PSR1_dashed=pre_S2_PSR1;ASP_used1_dashed=ASP_used1;W_used1_dashed=W_used1;A_used1_dashed=A_used1;RWR_used1_dashed=RWR_used1;
pre_S2_PSR2_dashed=pre_S2_PSR2;ASP_used2_dashed=ASP_used2;W_used2_dashed=W_used2;A_used2_dashed=A_used2;RWR_used2_dashed=RWR_used2;
pre_S2_PSR3_dashed=pre_S2_PSR3;ASP_used3_dashed=ASP_used3;W_used3_dashed=W_used3;A_used3_dashed=A_used3;RWR_used3_dashed=RWR_used3;
Mw_dashed=Mw;
index_tar_LargeM=find(abs(Mw-tar_LargeM)==min(abs(Mw-tar_LargeM)));
index_tar_SmallM=find(abs(Mw-tar_SmallM)==min(abs(Mw-tar_SmallM)));
markersize_tarRWR=7;
figure(1);
H1=subplot(3,3,1);
  plot(Mw(index_M55:end),log10(ASP_used1(index_M55:end)),'k- ','linewidth',2);hold on;
  plot(Mw(index_M55:end),log10(ASP_used2(index_M55:end)),'linestyle','-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(Mw(index_M55:end),log10(ASP_used3(index_M55:end)),'m- ','linewidth',1.5);
  title(title_in,'fontsize',14,'fontweight','bold');
  legend_in1=[legend_in1,'NZ_b_o_r^M^S=',num2str(NZbor_regional_mean_SSM_used1),')'];
  legend_in2=[legend_in2,'NZ_b_o_r^M^S=',num2str(NZbor_regional_mean_SSM_used2),')'];
  legend_in3=[legend_in3,'NZ_b_o_r^M^S=',num2str(NZbor_regional_mean_SSM_used3),')'];
  HL3=legend(legend_in1,legend_in2,legend_in3,3);
  po=get(HL3,'position');
  if(strcmp(Comptype,'ZBSZ')==1)
    if(strcmp(region,'SHI')==1)
      set(HL3,'position',[po(1)+0.51,po(2)-0.027,po(3)*0.7,po(4)*0.7],'FontSize',10,'fontweight','bold');
    else
      set(HL3,'position',[po(1)+0.516,po(2)-0.027,po(3)*0.7,po(4)*0.7],'FontSize',11,'fontweight','bold');
    end
  elseif(strcmp(Comptype,'dip')==1)
    if(strcmp(region,'SHI')==1)
      set(HL3,'position',[po(1)+0.512,po(2)-0.027,po(3)*0.7,po(4)*0.7],'FontSize',10,'fontweight','bold');
    else
      set(HL3,'position',[po(1)+0.52,po(2)-0.027,po(3)*0.7,po(4)*0.7],'FontSize',11,'fontweight','bold');
    end
  end
  plot(Mw_dashed,log10(ASP_used1_dashed),'k: ','linewidth',2);hold on;
  plot(Mw_dashed,log10(ASP_used2_dashed),'linestyle',':','color',[0/255 140/255 255/255],'linewidth',1);
  plot(Mw_dashed,log10(ASP_used3_dashed),'m: ','linewidth',1.5);
  plot(Mw(index_M55:end),log10(ASP_used1(index_M55:end)),'k- ','linewidth',2);hold on;
  plot(Mw(index_M55:end),log10(ASP_used2(index_M55:end)),'linestyle','-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(Mw(index_M55:end),log10(ASP_used3(index_M55:end)),'m- ','linewidth',1.5);  
%---- write M_{BP}^{RWR=1} on fig
  text1=['M_B_P^R^W^R^=^1:'];
  text2=[num2str(round(MBP_RWReq1_1*100)/100),sprintf('\n'),'{\color[rgb]{0 0.549 1}',num2str(round(MBP_RWReq1_2*100)/100),'}',sprintf('\n'),'{\color{magenta}',num2str(round(MBP_RWReq1_3*100)/100),'}'];
  text(0.06,0.8,text1,'units','normalized','fontweight','bold','fontsize',10,'backgroundcolor','white','edgecolor','none','margin',1);
  text(0.06,0.53,text2,'units','normalized','fontweight','bold','fontsize',10,'backgroundcolor','white','edgecolor','none','margin',1);
%---- plot aleatory variability of the SSM model
  colorcode(1,:)=[0/255 0/255 0/255];colorcode(2,:)=[0/255 140/255 255/255];colorcode(3,:)=[245/255 50/255 190/255];
  if(strcmp(region,'CA')==1&strcmp(Comptype,'ZBSZ')==1)
    target_M=[7.8,7.2,6.6,5.9,5.1];
  elseif(strcmp(region,'TW')==1&strcmp(Comptype,'dip')==1)
    target_M=[7.9,7.0,6.3,5.6,5.1];
  elseif(strcmp(region,'SHI')==1&strcmp(Comptype,'dip')==1)
    target_M=[7.9,7.1,6.3,5.7,5.1];
  elseif(strcmp(region,'SHI')==1&strcmp(Comptype,'ZBSZ')==1)
    target_M=[7.9,7.1,6.3,5.7,5.1];
  else
    target_M=[7.9,7.3,6.3,5.7,5.1];
  end
  M_interval=0.12;bar_long=0.05;
  for kk=1:length(target_M)
      target_M_plot(1,kk)=target_M(kk)-M_interval;target_M_plot(2,kk)=target_M(kk);
      target_M_plot(3,kk)=target_M(kk)+M_interval;
      if(target_M_plot(3,kk)>8)
        target_M_plot(3,kk)=8;
      end
  end
  for ik=1:length(target_M)
    for jk=1:3 
      index_Mtar_loc=find(abs(Mw-target_M_plot(jk,ik))==min(abs(Mw-target_M_plot(jk,ik))));
      eval(['target_AV_range(1)=log10(ASP_addsigma',num2str(jk),'(index_Mtar_loc));']);
      eval(['target_AV_range(2)=log10(ASP_minussigma',num2str(jk),'(index_Mtar_loc));']);
      eval(['[nu]=sub_pt_AV_range(target_M_plot(jk,ik),target_AV_range,colorcode(',num2str(jk),',:),bar_long);']);
    end
  end
%----
  xlabel('M','fontsize',13,'fontweight','bold');
  ylabel('Asp (log_1_0 unit)','fontsize',13,'fontweight','bold');
  xlim([4.9 8.2]);
  ylim([log10(0.5) log10(100)]);
  grid on;
H2=subplot(3,3,2);
  plot(Mw(index_M55:end),W_used1(index_M55:end),'k- ','linewidth',2);hold on;
  plot(Mw(index_M55:end),W_used2(index_M55:end),'linestyle','-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(Mw(index_M55:end),W_used3(index_M55:end),'m- ','linewidth',1.5);
%-- replot for evaluating extrapolating part
  plot(Mw_dashed,W_used1_dashed,'k: ','linewidth',2);hold on;
  plot(Mw_dashed,W_used2_dashed,'linestyle',':','color',[0/255 140/255 255/255],'linewidth',1);
  plot(Mw_dashed,W_used3_dashed,'m: ','linewidth',1.5);
  plot(Mw(index_M55:end),W_used1(index_M55:end),'k- ','linewidth',2);hold on;
  plot(Mw(index_M55:end),W_used2(index_M55:end),'linestyle','-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(Mw(index_M55:end),W_used3(index_M55:end),'m- ','linewidth',1.5);
%--
  if(strcmp(region,'CA')==1&strcmp(Comptype,'ZBSZ')==1)
  elseif(strcmp(region,'CA')==1&strcmp(Comptype,'dip')==1)
    plot(Mw(index_tar_LargeM),W_used1(index_tar_LargeM),'ko ','markersize',markersize_tarRWR,'linewidth',1);
    plot(Mw(index_tar_LargeM),W_used2(index_tar_LargeM),'o','color',[0/255 140/255 255/255],'markersize',markersize_tarRWR,'linewidth',1);
    plot(Mw(index_tar_LargeM),W_used3(index_tar_LargeM),'mo ','markersize',markersize_tarRWR,'linewidth',1);
  elseif(strcmp(region,'SHI')==1&strcmp(Comptype,'dip')==1)
    plot(Mw(index_tar_LargeM),W_used1(index_tar_LargeM),'ko ','markersize',markersize_tarRWR,'linewidth',1);
    plot(Mw(index_tar_LargeM),W_used2(index_tar_LargeM),'o','color',[0/255 140/255 255/255],'markersize',markersize_tarRWR,'linewidth',1);
    plot(Mw(index_tar_LargeM),W_used3(index_tar_LargeM),'mo ','markersize',markersize_tarRWR,'linewidth',1);
  elseif(strcmp(region,'SHI')==1&strcmp(Comptype,'ZBSZ')==1)
    plot(Mw(index_tar_LargeM),W_used1(index_tar_LargeM),'ko ','markersize',markersize_tarRWR,'linewidth',1);
    plot(Mw(index_tar_LargeM),W_used2(index_tar_LargeM),'o','color',[0/255 140/255 255/255],'markersize',markersize_tarRWR+1,'linewidth',1);
    plot(Mw(index_tar_LargeM),W_used3(index_tar_LargeM),'mo ','markersize',markersize_tarRWR-2,'linewidth',1);
  elseif(strcmp(region,'TW')==1&strcmp(Comptype,'ZBSZ')==1)
    plot(Mw(index_tar_LargeM),W_used1(index_tar_LargeM),'ko ','markersize',markersize_tarRWR,'linewidth',1);
    plot(Mw(index_tar_LargeM),W_used2(index_tar_LargeM),'o','color',[0/255 140/255 255/255],'markersize',markersize_tarRWR+1,'linewidth',1);
    plot(Mw(index_tar_LargeM),W_used3(index_tar_LargeM),'mo ','markersize',markersize_tarRWR-2,'linewidth',1);
  else
    plot(Mw(index_tar_LargeM),W_used1(index_tar_LargeM),'ko ','markersize',markersize_tarRWR+4,'linewidth',1);
    plot(Mw(index_tar_LargeM),W_used2(index_tar_LargeM),'o','color',[0/255 140/255 255/255],'markersize',markersize_tarRWR+1,'linewidth',1);
    plot(Mw(index_tar_LargeM),W_used3(index_tar_LargeM),'mo ','markersize',markersize_tarRWR-2,'linewidth',1);
  end
  if(strcmp(region,'CA')==1&strcmp(Comptype,'ZBSZ')==1)
    plot(Mw(index_tar_SmallM),W_used1(index_tar_SmallM),'k^ ','markersize',markersize_tarRWR,'linewidth',1);
    plot(Mw(index_tar_SmallM),W_used2(index_tar_SmallM),'^','color',[0/255 140/255 255/255],'markersize',markersize_tarRWR+1,'linewidth',1);
    plot(Mw(index_tar_SmallM),W_used3(index_tar_SmallM),'m^ ','markersize',markersize_tarRWR-2,'linewidth',1);
  else
    plot(Mw(index_tar_SmallM),W_used1(index_tar_SmallM),'k^ ','markersize',markersize_tarRWR+4,'linewidth',1);
    plot(Mw(index_tar_SmallM),W_used2(index_tar_SmallM),'^','color',[0/255 140/255 255/255],'markersize',markersize_tarRWR+1,'linewidth',1);
    plot(Mw(index_tar_SmallM),W_used3(index_tar_SmallM),'m^ ','markersize',markersize_tarRWR-2,'linewidth',1);
  end
  xlabel('M','fontsize',13,'fontweight','bold');
  ylabel('W_R_u_p (km)','fontsize',13,'fontweight','bold');
%---- plot the aleatory variability of the SSM model
  bar_long=0.05;
  for ik=1:length(target_M)
    for jk=1:3
      index_Mtar_loc=find(abs(Mw-target_M_plot(jk,ik))==min(abs(Mw-target_M_plot(jk,ik))));
      eval(['target_AV_range(1)=W_add2sigmaal',num2str(jk),'(index_Mtar_loc);']);
      eval(['target_AV_range(2)=W_minus2sigmaal',num2str(jk),'(index_Mtar_loc);']);
      eval(['[nu]=sub_pt_AV_range(target_M_plot(jk,ik),target_AV_range,colorcode(',num2str(jk),',:),bar_long);']);
    end
  end
%----
  xlim([4.9 8.2]);%ylim([0 30]);
  ymax=max(max([W_add2sigmaal1;W_add2sigmaal2;W_add2sigmaal3]));
  ymax_round=round(ymax/5);
  ylim([0 (ymax_round*5)+5]);
  grid on
H4=axes('position',[0.15,0.44,0.195,0.195]);
  plot(Mw(index_M55:end),log10(A_used1(index_M55:end)),'k- ','linewidth',2);hold on;
  plot(Mw(index_M55:end),log10(A_used3(index_M55:end)),'m- ','linewidth',1.5);
  plot(Mw(index_M55:end),log10(A_used2(index_M55:end)),'linestyle','-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(Mw_dashed,log10(A_used1_dashed),'k: ','linewidth',2);
  plot(Mw_dashed,log10(A_used3_dashed),'m: ','linewidth',1.5);
  plot(Mw_dashed,log10(A_used2_dashed),'linestyle',':','color',[0/255 140/255 255/255],'linewidth',1);
  plot(Mw(index_M55:end),log10(A_used1(index_M55:end)),'k- ','linewidth',2);hold on;
  plot(Mw(index_M55:end),log10(A_used3(index_M55:end)),'m- ','linewidth',1.5);
  plot(Mw(index_M55:end),log10(A_used2(index_M55:end)),'linestyle','-','color',[0/255 140/255 255/255],'linewidth',1);
  xlabel('M','fontsize',13,'fontweight','bold');
  ylabel('A (km^2, log_1_0 unit)','fontsize',13,'fontweight','bold');
  xlim([4.9 8.2]);
  ylim([0 5]); %log10 unit
%---- plot the aleatory variability of the SSM model
  bar_long=0.05;
  for ik=1:length(target_M)
    for jk=1:3
      index_Mtar_loc=find(abs(Mw-target_M_plot(jk,ik))==min(abs(Mw-target_M_plot(jk,ik))));
      eval(['target_AV_range(1)=log10(A_addsigma',num2str(jk),'(index_Mtar_loc));']);
      eval(['target_AV_range(2)=log10(A_minussigma',num2str(jk),'(index_Mtar_loc));']);
      eval(['[nu]=sub_pt_AV_range(target_M_plot(jk,ik),target_AV_range,colorcode(',num2str(jk),',:),bar_long);']);
    end
  end
%----  
  grid on;
H5=axes('position',[0.44,0.15,0.46,0.46]);
  plot(Mw(index_M55:end),log10(pre_S2_PSR1(index_M55:end)),'k- ','linewidth',2);hold on;
  plot(Mw(index_M55:end),log10(pre_S2_PSR2(index_M55:end)),'linestyle','-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(Mw(index_M55:end),log10(pre_S2_PSR3(index_M55:end)),'m- ','linewidth',1.5);
  plot(Mw_dashed,log10(pre_S2_PSR1_dashed),'k: ','linewidth',2);hold on;
  plot(Mw_dashed,log10(pre_S2_PSR2_dashed),'linestyle',':','color',[0/255 140/255 255/255],'linewidth',1);
  plot(Mw_dashed,log10(pre_S2_PSR3_dashed),'m: ','linewidth',1.5);
  plot(Mw(index_M55:end),log10(pre_S2_PSR1(index_M55:end)),'k- ','linewidth',2);hold on;
  plot(Mw(index_M55:end),log10(pre_S2_PSR2(index_M55:end)),'linestyle','-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(Mw(index_M55:end),log10(pre_S2_PSR3(index_M55:end)),'m- ','linewidth',1.5);
  if(strcmp(region,'TW')==1|strcmp(region,'SHI')==1)
    HL52=legend(['W_R_u_p: Hea24,',sprintf('\n'),'          regional NZ_H_y_p_o,',sprintf('\n'),'          W_F_l_t^F^Z^S=',num2str(round(W_flt_FZS1*10)/10),'km'],['...     , W_F_l_t^F^Z^S=',num2str(round(W_flt_FZS2*10)/10),'km'],['...     , W_F_l_t^F^Z^S=',num2str(round(W_flt_FZS3*10)/10),'km'],4);
    po=get(HL52,'position');set(HL52,'position',[po(1)+0.2163,po(2)+0.025,po(3)*0.6,po(4)*0.6],'FontSize',11,'fontweight','bold');
    annotation('line',[0.9513 0.962],[0.359 0.359],'color','k','linewidth',1.5);
  else
    HL52=legend(['W_R_u_p: Hea24, regional NZ_H_y_p_o,',sprintf('\n'),'          W_F_l_t^F^Z^S=',num2str(round(W_flt_FZS1*10)/10),'km'],['...     , W_F_l_t^F^Z^S=',num2str(round(W_flt_FZS2*10)/10),'km'],['...     , W_F_l_t^F^Z^S=',num2str(round(W_flt_FZS3*10)/10),'km'],4);
    po=get(HL52,'position');set(HL52,'position',[po(1)+0.23,po(2)+0.025,po(3)*0.6,po(4)*0.6],'FontSize',11,'fontweight','bold');
    annotation('line',[0.9511 0.9618],[0.3617 0.3617],'color','k','linewidth',1);
  end
%---- plot aleatory variability of the SSM model
  bar_long=0.05;
  for ik=1:length(target_M)
    for jk=1:3
      index_Mtar_loc=find(abs(Mw-target_M_plot(jk,ik))==min(abs(Mw-target_M_plot(jk,ik))));
      eval(['target_AV_range(1)=log10(pre_S2_PSR_add2sigmaal',num2str(jk),'(index_Mtar_loc));']);
      eval(['target_AV_range(2)=log10(pre_S2_PSR_minus2sigmaal',num2str(jk),'(index_Mtar_loc));']);
      eval(['[nu]=sub_pt_AV_range(target_M_plot(jk,ik),target_AV_range,colorcode(',num2str(jk),',:),bar_long);']);
    end
  end
%----
  xlabel('M','fontsize',13,'fontweight','bold');
  ylabel('Probability (S2 PSR)','fontsize',13,'fontweight','bold');
  xlim([4.9 8.2]);
  ylim([-2 0]);
  grid on;
H6=axes('position',[0.15,0.15,0.195,0.195]);
  plot(RWR_used1(index_M55:end),pre_S2_PSR1(index_M55:end),'k- ','linewidth',2);hold on;
  plot(RWR_used2(index_M55:end),pre_S2_PSR2(index_M55:end),'linestyle','-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(RWR_used3(index_M55:end),pre_S2_PSR3(index_M55:end),'m- ','linewidth',1.5);
  plot(RWR_used1_dashed,pre_S2_PSR1_dashed,'k: ','linewidth',2);hold on;
  plot(RWR_used2_dashed,pre_S2_PSR2_dashed,'linestyle',':','color',[0/255 140/255 255/255],'linewidth',1);
  plot(RWR_used3_dashed,pre_S2_PSR3_dashed,'m: ','linewidth',1.5);
  plot(RWR_used1(index_M55:end),pre_S2_PSR1(index_M55:end),'k- ','linewidth',2);hold on;
  plot(RWR_used2(index_M55:end),pre_S2_PSR2(index_M55:end),'linestyle','-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(RWR_used3(index_M55:end),pre_S2_PSR3(index_M55:end),'m- ','linewidth',1.5);
  if(strcmp(region,'CA')==1&strcmp(Comptype,'ZBSZ')==1)
  else
    plot(RWR_used1(index_tar_LargeM),pre_S2_PSR1(index_tar_LargeM),'ko ','markersize',markersize_tarRWR-1,'linewidth',1);
    plot(RWR_used2(index_tar_LargeM),pre_S2_PSR2(index_tar_LargeM),'o','color',[0/255 140/255 255/255],'markersize',markersize_tarRWR-1,'linewidth',1);
    plot(RWR_used3(index_tar_LargeM),pre_S2_PSR3(index_tar_LargeM),'mo ','markersize',markersize_tarRWR-1,'linewidth',1);
  end
  plot(RWR_used1(index_tar_SmallM),pre_S2_PSR1(index_tar_SmallM),'k^ ','markersize',markersize_tarRWR,'linewidth',1);
  plot(RWR_used2(index_tar_SmallM),pre_S2_PSR2(index_tar_SmallM),'^','color',[0/255 140/255 255/255],'markersize',markersize_tarRWR,'linewidth',1);
  plot(RWR_used3(index_tar_SmallM),pre_S2_PSR3(index_tar_SmallM),'m^ ','markersize',markersize_tarRWR,'linewidth',1);
  xlabel('RWR','fontsize',13,'fontweight','bold');
  ylabel('Probability (S2 PSR)','fontsize',13,'fontweight','bold');
  xlim([0 1]);ylim([0 1]);
  grid on;
  if(strcmp(region,'CA')==1)
    NMZHypo_print='0.5122';
  elseif(strcmp(region,'TW')==1)
    NMZHypo_print='0.5198';
  elseif(strcmp(region,'SHI')==1)
    NMZHypo_print='0.4055';
  elseif(strcmp(region,'NWCN')==1)
    NMZHypo_print='0.3504';
  elseif(strcmp(region,'SEU')==1)
    NMZHypo_print='0.6116';
  end
  text(0.1,0.85,['NZ_H_y_p_o: ',NMZHypo_print],'fontsize',10,'fontweight','bold','backgroundcolor','white','edgecolor','none','margin',1);
  annotation('line',[0.1826 0.192],[0.337 0.337],'color','k','linewidth',1.2);
labelfontsize=12;
set(H1,'xtick',[5:0.5:8.5],'fontsize',labelfontsize,'fontweight','bold','linewidth',1);
set(H2,'ytick',[0:5:60],'yticklabel',[{0},{''},{10},{''},{20},{''},{30},{''},{40},{''},{50},{''},{60}],'xtick',[5:0.5:8.5],'fontsize',labelfontsize,'fontweight','bold','linewidth',1);
set(H4,'xtick',[5:0.5:8.5],'ytick',[0:6],'yticklabel',[{0},{''},{2},{''},{4},{''},{6}],'fontsize',labelfontsize,'fontweight','bold','linewidth',1);
ytick_lin=[0.01:0.01:0.1,0.2:0.1:1];ytick_log=log10(ytick_lin);
set(H5,'fontsize',labelfontsize,'fontweight','bold','xtick',[5:0.5:8],'ytick',ytick_log,'yticklabel',[{0.01},{''},{''},{''},{0.05},{''},{''},{''},{''},{0.1},{''},{''},{''},{0.5},{''},{''},{''},{''},{1}],'linewidth',1);
set(H6,'ytick',[0:0.2:1],'xtick',[0:0.2:1],'fontsize',labelfontsize,'fontweight','bold','linewidth',1);
print(gcf,outnm,'-deps2c');unix(['convert -density 200 ',outnm,'.eps ',outnm,'.png']);
close(1);










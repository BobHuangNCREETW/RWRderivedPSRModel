function [nu]=sub_pt_RWRPSR_w_other_models(PSR_RWR_P1,PSR_RWR_P_Psigmaal1,PSR_RWR_P_Msigmaal1,PSR_RWR_P2,PSR_RWR_P_Psigmaal2,PSR_RWR_P_Msigmaal2,PSR_RWR_P3,PSR_RWR_P_Psigmaal3,PSR_RWR_P_Msigmaal3,PSR_RWR_M1,PSR_RWR_M2,PSR_RWR_M3,Mea24_CPSR_P1,Mea24_CPSR_P2,Mea24_CPSR_P3,Mea24_CPSR_M1,Mea24_CPSR_M2,Mea24_CPSR_M3,legend_in1,legend_in2,legend_in3,legend_inCPSR1,legend_inCPSR2,legend_inCPSR3,outnm,Comp_SoF);
%                                     by Bob J.Y. Huang in March 2026
nu='';
%---- Generate PSR from Pizza et al. (2023) with various fault types.
example_M_Pea23=4.7:0.1:7.9;
c0_ALL_M_Pea23=-14.47;c1_ALL_M_Pea23=2.177; % from Pizza et al. 2023 Tab 2 all type
pre_ztmp_model_Mw_Pea23_ALL=c0_ALL_M_Pea23+(example_M_Pea23.*c1_ALL_M_Pea23);
logitfit_model_Mw_Pea23_ALL=(exp(pre_ztmp_model_Mw_Pea23_ALL))./(1+(exp(pre_ztmp_model_Mw_Pea23_ALL)));
c0_NML_M_Pea23=-13.5;c1_NML_M_Pea23=2.159; % NML
pre_ztmp_model_Mw_Pea23_NML=c0_NML_M_Pea23+(example_M_Pea23.*c1_NML_M_Pea23);
logitfit_model_Mw_Pea23_NML=(exp(pre_ztmp_model_Mw_Pea23_NML))./(1+(exp(pre_ztmp_model_Mw_Pea23_NML)));
c0_RV_M_Pea23=-10.75;c1_RV_M_Pea23=1.427; % RV
pre_ztmp_model_Mw_Pea23_RV=c0_RV_M_Pea23+(example_M_Pea23.*c1_RV_M_Pea23);
logitfit_model_Mw_Pea23_RV=(exp(pre_ztmp_model_Mw_Pea23_RV))./(1+(exp(pre_ztmp_model_Mw_Pea23_RV)));
c0_SS_M_Pea23=-28.56;c1_SS_M_Pea23=4.436; % SS
pre_ztmp_model_Mw_Pea23_SS=c0_SS_M_Pea23+(example_M_Pea23.*c1_SS_M_Pea23);
logitfit_model_Mw_Pea23_SS=(exp(pre_ztmp_model_Mw_Pea23_SS))./(1+(exp(pre_ztmp_model_Mw_Pea23_SS)));
%---- Generate PSR from Youngs et al. (2003) for all types.
example_M_Yea03=4.7:0.01:8.13;
c0_S1_M_Yea03=-12.51;c1_S1_M_Yea03=2.053; % description of Tab 1 in Pizza et al. 2023; dataset from Wells and Coppersmith 1993
pre_ztmp_model_Mw_Yea03=c0_S1_M_Yea03+(example_M_Yea03.*c1_S1_M_Yea03);
logitfit_model_Mw_Yea03=(exp(pre_ztmp_model_Mw_Yea03))./(1+(exp(pre_ztmp_model_Mw_Yea03)));
%---- Generate PSR based on Moss and Ross (2011) RV types.
example_M_MR11=5.5:0.1:8;
c0_S1_M_MR11=7.3;c1_S1_M_MR11=-1.03; % description of Eq. (5) of Moss and Ross (2011)
pre_ztmp_model_Mw_MR11=c0_S1_M_MR11+(example_M_MR11.*c1_S1_M_MR11);
%logitfit_model_Mw_MR11=(exp(pre_ztmp_model_Mw_MR11))./(1+(exp(pre_ztmp_model_Mw_MR11)));
logitfit_model_Mw_MR11=(1./(1+(exp(pre_ztmp_model_Mw_MR11))));
%---- 
if(strcmp(Comp_SoF,'SS')==1)
  Comp_P1=logitfit_model_Mw_Pea23_SS;Comp_M1=example_M_Pea23;
  Comp_P2=logitfit_model_Mw_Yea03;Comp_M2=example_M_Yea03;
  Leg_Comp1='P(M, Pea23, SS type)';Leg_Comp2='P(M, Yea03, all type, data from WC93)';
elseif(strcmp(Comp_SoF,'RV')==1)
  Comp_P1=logitfit_model_Mw_Pea23_RV;Comp_M1=example_M_Pea23;
  Comp_P2=logitfit_model_Mw_Yea03;Comp_M2=example_M_Yea03;
  Comp_P3=logitfit_model_Mw_MR11;Comp_M3=example_M_MR11;
  Leg_Comp1='P(M, Pea23, RV type)';Leg_Comp2='P(M, Yea03, all type, data from WC93)';
  Leg_Comp3='P(M, MR11, RV type)';
elseif(strcmp(Comp_SoF,'NML')==1)
  Comp_P1=logitfit_model_Mw_Pea23_NML;Comp_M1=example_M_Pea23;
  Comp_P2=logitfit_model_Mw_Yea03;Comp_M2=example_M_Yea03;
  Leg_Comp1='P(M, Pea23, NML type)';Leg_Comp2='P(M, Yea03, all type, data from WC93)';
  Leg_Comp1='P(M, Pea23, NML type)';Leg_Comp2='P(M, Yea03, all type, data from WC93)';
end  
%%% evaluate regions
index_dash=find(outnm=='_');
region=outnm(index_dash(1)+1:index_dash(2)-1);
Comptype=outnm(index_dash(2)+1:index_dash(3)-1);
%%%
%---- grab Mw=5.5 because Hea24 started from here
index_M55=find(abs(PSR_RWR_M1-5.5)==min(abs(PSR_RWR_M1-5.5)));
%----
figure(1);
H1=subplot(2,1,1);
  plot(PSR_RWR_M1(index_M55:end),log10(PSR_RWR_P1(index_M55:end)),'-k ','linewidth',2);hold on;
  plot(PSR_RWR_M2(index_M55:end),log10(PSR_RWR_P2(index_M55:end)),'-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(PSR_RWR_M3(index_M55:end),log10(PSR_RWR_P3(index_M55:end)),'m- ','linewidth',1.5);
  plot(Mea24_CPSR_M1,log10(Mea24_CPSR_P1),'-o','color',[225/255 180/255 10/255],'linewidth',1);
  plot(Mea24_CPSR_M2,log10(Mea24_CPSR_P2),'-s ','color',[225/255 180/255 10/255],'linewidth',1);
  plot(Mea24_CPSR_M3,log10(Mea24_CPSR_P3),'-^ ','color',[225/255 180/255 10/255],'linewidth',1);
  if(strcmp(Comp_SoF,'RV')==1)
    plot(Comp_M1,log10(Comp_P1),'linestyle','-','color','g','linewidth',1.5);
    plot(Comp_M3,log10(Comp_P3),'linestyle','-','color',[0.7,0.7,0.7],'linewidth',1.5);
    plot(Comp_M2,log10(Comp_P2),'linestyle','--','color',[0.4,0.4,0.4],'linewidth',1.5);
    LH=legend(legend_in1,legend_in2,legend_in3,legend_inCPSR1,legend_inCPSR2,legend_inCPSR3,Leg_Comp1,Leg_Comp3,Leg_Comp2,4);
    if(strcmp(Comptype,'dip')==1)
      po=get(LH,'position');set(LH,'position',[po(1)+0.105,po(2)-0.03,po(3)*0.7,po(4)*0.7],'FontSize',12,'fontweight','bold');
    else
      po=get(LH,'position');set(LH,'position',[po(1)+0.09,po(2)-0.03,po(3)*0.7,po(4)*0.7],'FontSize',12,'fontweight','bold');
    end
  else
    plot(Comp_M1,log10(Comp_P1),'linestyle','-','color','g','linewidth',1.5);
    plot(Comp_M2,log10(Comp_P2),'linestyle','--','color',[0.4,0.4,0.4],'linewidth',1.5);
    LH=legend(legend_in1,legend_in2,legend_in3,legend_inCPSR1,legend_inCPSR2,legend_inCPSR3,Leg_Comp1,Leg_Comp2,4);
    if(strcmp(region,'SHI')==1)
      if(strcmp(Comptype,'dip')==1)
        po=get(LH,'position');set(LH,'position',[po(1)+0.115,po(2)-0.03,po(3)*0.7,po(4)*0.7],'FontSize',12,'fontweight','bold');
      elseif(strcmp(Comptype,'ZBSZ')==1)
        po=get(LH,'position');set(LH,'position',[po(1)+0.105,po(2)-0.03,po(3)*0.7,po(4)*0.7],'FontSize',12,'fontweight','bold');
      end
    else
      if(strcmp(Comptype,'dip')==1)
        po=get(LH,'position');set(LH,'position',[po(1)+0.10,po(2)-0.03,po(3)*0.7,po(4)*0.7],'FontSize',12,'fontweight','bold');
      elseif(strcmp(Comptype,'ZBSZ')==1)
        po=get(LH,'position');set(LH,'position',[po(1)+0.085,po(2)-0.03,po(3)*0.7,po(4)*0.7],'FontSize',12,'fontweight','bold');
      end
    end
  end
  plot(PSR_RWR_M1,log10(PSR_RWR_P1),':k ','linewidth',2);hold on;
  plot(PSR_RWR_M2,log10(PSR_RWR_P2),':','color',[0/255 140/255 255/255],'linewidth',1);
  plot(PSR_RWR_M3,log10(PSR_RWR_P3),'m: ','linewidth',1.5);
  plot(PSR_RWR_M1(index_M55:end),log10(PSR_RWR_P1(index_M55:end)),'-k ','linewidth',2);
  plot(PSR_RWR_M2(index_M55:end),log10(PSR_RWR_P2(index_M55:end)),'-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(PSR_RWR_M3(index_M55:end),log10(PSR_RWR_P3(index_M55:end)),'m- ','linewidth',1.5);
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
  if(strcmp(region,'CA')==1&strcmp(Comptype,'dip')==1)
    MBP_RWReq1_3=6.39;
  elseif(strcmp(region,'CA')==1&strcmp(Comptype,'ZBSZ')==1)  
    MBP_RWReq1_3=6.39;
  elseif(strcmp(region,'TW')==1&strcmp(Comptype,'dip')==1)
    MBP_RWReq1_3=6.68;
  elseif(strcmp(region,'TW')==1&strcmp(Comptype,'ZBSZ')==1)
    MBP_RWReq1_3=7.01;
  elseif(strcmp(region,'SHI')==1&strcmp(Comptype,'dip')==1)
    MBP_RWReq1_3=6.67;
  elseif(strcmp(region,'SHI')==1&strcmp(Comptype,'ZBSZ')==1)
    MBP_RWReq1_3=7.08;
  else
    MBP_RWReq1_3=7;
  end
  for kk=1:length(target_M)
      target_M_plot(1,kk)=target_M(kk)-M_interval;target_M_plot(2,kk)=target_M(kk);
      target_M_plot(3,kk)=target_M(kk)+M_interval;
      if(target_M_plot(3,kk)>8)
        target_M_plot(3,kk)=8;
      end
  end
  for ik=1:length(target_M)
    for jk=1:3
      index_Mtar_loc=find(abs(PSR_RWR_M1-target_M_plot(jk,ik))==min(abs(PSR_RWR_M1-target_M_plot(jk,ik))));
      eval(['target_AV_range(1)=log10(PSR_RWR_P_Psigmaal',num2str(jk),'(index_Mtar_loc));']); 
      eval(['target_AV_range(2)=log10(PSR_RWR_P_Msigmaal',num2str(jk),'(index_Mtar_loc));']);
      eval(['[nu]=sub_pt_AV_range(target_M_plot(jk,ik),target_AV_range,colorcode(',num2str(jk),',:),bar_long);']);
    end
  end 
%----
  grid on;
  xlabel('M','FontSize',16,'fontweight','bold');
  ylabel('Probability','FontSize',16,'fontweight','bold');
  xlim([4.8,8.2]);
  ylim([-2 0]);
H2=subplot(2,1,2);
  plot(PSR_RWR_M1(index_M55:end),PSR_RWR_P1(index_M55:end),'-k ','linewidth',2);hold on;
  plot(PSR_RWR_M2(index_M55:end),PSR_RWR_P2(index_M55:end),'-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(PSR_RWR_M3(index_M55:end),PSR_RWR_P3(index_M55:end),'m- ','linewidth',1.5);
  plot(Mea24_CPSR_M1,Mea24_CPSR_P1,'-o','color',[225/255 180/255 10/255],'linewidth',1);
  plot(Mea24_CPSR_M2,Mea24_CPSR_P2,'-s ','color',[225/255 180/255 10/255],'linewidth',1);
  plot(Mea24_CPSR_M3,Mea24_CPSR_P3,'-^ ','color',[225/255 180/255 10/255],'linewidth',1);
  if(strcmp(Comp_SoF,'RV')==1)
    plot(Comp_M1,Comp_P1,'linestyle','-','color','g','linewidth',1.5);
    plot(Comp_M3,Comp_P3,'linestyle','-','color',[0.7,0.7,0.7],'linewidth',1.5);
    plot(Comp_M2,Comp_P2,'linestyle','--','color',[0.4,0.4,0.4],'linewidth',1.5);
  else
    plot(Comp_M1,Comp_P1,'linestyle','-','color','g','linewidth',1.5);
    plot(Comp_M2,Comp_P2,'linestyle','--','color',[0.4,0.4,0.4],'linewidth',1.5);
  end
  plot(PSR_RWR_M1,PSR_RWR_P1,':k ','linewidth',2);hold on;
  plot(PSR_RWR_M2,PSR_RWR_P2,':','color',[0/255 140/255 255/255],'linewidth',1);
  plot(PSR_RWR_M3,PSR_RWR_P3,'m: ','linewidth',1.5);
  plot(PSR_RWR_M1(index_M55:end),PSR_RWR_P1(index_M55:end),'-k ','linewidth',2);
  plot(PSR_RWR_M2(index_M55:end),PSR_RWR_P2(index_M55:end),'-','color',[0/255 140/255 255/255],'linewidth',1);
  plot(PSR_RWR_M3(index_M55:end),PSR_RWR_P3(index_M55:end),'m- ','linewidth',1.5);
%---- plot aleatory variability of the SSM model
  for ik=1:length(target_M)
    for jk=1:3
      index_Mtar_loc=find(abs(PSR_RWR_M1-target_M_plot(jk,ik))==min(abs(PSR_RWR_M1-target_M_plot(jk,ik))));
      eval(['target_AV_range(1)=PSR_RWR_P_Psigmaal',num2str(jk),'(index_Mtar_loc);']);
      eval(['target_AV_range(2)=PSR_RWR_P_Msigmaal',num2str(jk),'(index_Mtar_loc);']);
      eval(['[nu]=sub_pt_AV_range(target_M_plot(jk,ik),target_AV_range,colorcode(',num2str(jk),',:),bar_long);']);
    end
  end
%----
  grid on;
  xlabel('M','FontSize',16,'fontweight','bold');
  ylabel('Probability','FontSize',16,'fontweight','bold');
  xlim([4.8,8.2]);
  ylim([0 1]);
%-----
ytick_lin=[0.01:0.01:0.1,0.2:0.1:1];ytick_log=log10(ytick_lin);
set(H1,'position',[0.15,0.09,0.33,0.39],'fontsize',12,'fontweight','bold','xtick',[5:0.5:8],'ytick',ytick_log,'yticklabel',[{0.01},{''},{''},{''},{0.05},{''},{''},{''},{''},{0.1},{''},{''},{''},{0.5},{''},{''},{''},{''},{1}]);
set(H2,'position',[0.56,0.09,0.33,0.39],'fontsize',12,'fontweight','bold','xtick',[5:0.5:8]);
print(gcf,outnm,'-deps2c');unix(['convert -density 200 ',outnm,'.eps ',outnm,'.png']);
close all;





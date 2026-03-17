%-----
%      New empirical RWR-derived model for the probability of surface rupture (PSR)
%      RWR:        Rupture-width ratio, defined as W_Rup divided by W_Flt.
%      W_Rup:      Rupture width.
%      W_Flt:      Fault width, considering fault geometry such as the fault dip angle 
%                    and seismogenic thickness.
%      NZbar_Hypo: The regional (fault-specific) average normalized hypocenter depth.
%                                              by Bob J.Y. Huang in March 2026
%-----
assummed_RWR_case1=0:0.01:1;
assummed_NZbar_Hypo_case1=0.4733; % the average NZbar_Hypo in Table 1
for i=1:length(assummed_RWR_case1)
  [PSR_pred_case1(i),PSR_lo_case1(i),PSR_hi_case1(i)]=sub_RWR_derived_PSRmodel(assummed_RWR_case1(i),assummed_NZbar_Hypo_case1);
end
figure(1);
H1=subplot(1,1,1);
  plot(assummed_RWR_case1,PSR_pred_case1,'k- ','linewidth',4);hold on;
  plot(assummed_RWR_case1,PSR_lo_case1,'k- ','linewidth',1);
  plot(assummed_RWR_case1,PSR_hi_case1,'k- ','linewidth',1);
  xlabel('RWR (W_R_u_p/W_F_l_t)','FontSize',18,'fontweight','bold');
  ylabel('Probability (P(RWR, NZ_H_y_p_o; c_0, c_1, \beta))','FontSize',18,'fontweight','bold');
  annotation('line',[0.032 0.032],[0.58 0.607],'color','k','linewidth',2);
  legend('NZ_H_y_p_o=0.4733',2);
  annotation('line',[0.227 0.241],[0.91 0.91],'color','k','linewidth',1.5);
  grid on;set(H1,'FontSize',15,'fontweight','bold');
  outnm='RWR_PSR_S2_prediction_backbone';
  print(gcf,outnm,'-deps2c');unix(['convert -density 200 ',outnm,'.eps ',outnm,'.png']);

assummed_RWR_case2=0:0.01:1;
assummed_NZbar_Hypo_case2=0.3504; % minimum value of the empirical data 
for i=1:length(assummed_RWR_case2)
  [PSR_pred_case2(i),PSR_lo_case2(i),PSR_hi_case2(i)]=sub_RWR_derived_PSRmodel(assummed_RWR_case2(i),assummed_NZbar_Hypo_case2);
end

assummed_RWR_case3=0:0.01:1;
assummed_NZbar_Hypo_case3=0.6116; % maximum value of the empirical data
for i=1:length(assummed_RWR_case3)
  [PSR_pred_case3(i),PSR_lo_case3(i),PSR_hi_case3(i)]=sub_RWR_derived_PSRmodel(assummed_RWR_case3(i),assummed_NZbar_Hypo_case3);
end

figure(2);
H2=subplot(1,1,1);
  plot(assummed_RWR_case2,PSR_pred_case2,'b- ','linewidth',4);hold on;
  plot(assummed_RWR_case3,PSR_pred_case3,'m- ','linewidth',4);
  plot(assummed_RWR_case2,PSR_lo_case2,'b- ','linewidth',1);
  plot(assummed_RWR_case2,PSR_hi_case2,'b- ','linewidth',1);
  plot(assummed_RWR_case3,PSR_lo_case3,'m- ','linewidth',1);
  plot(assummed_RWR_case3,PSR_hi_case3,'m- ','linewidth',1);
  xlabel('RWR (W_R_u_p/W_F_l_t)','FontSize',18,'fontweight','bold');
  ylabel('Probability (P(RWR, NZ_H_y_p_o; c_0, c_1, \beta))','FontSize',18,'fontweight','bold');
  annotation('line',[0.032 0.032],[0.58 0.607],'color','k','linewidth',2);
  HL100=legend('NZ_H_y_p_o=0.3504','NZ_H_y_p_o=0.6116',2);
  annotation('line',[0.227 0.241],[0.851 0.851],'color','k','linewidth',1.5);
  annotation('line',[0.227 0.241],[0.91 0.91],'color','k','linewidth',1.5);
  grid on;set(H2,'FontSize',15,'fontweight','bold');
  outnm='RWR_PSR_S2_prediction_examples';
  print(gcf,outnm,'-deps2c');unix(['convert -density 200 ',outnm,'.eps ',outnm,'.png']);





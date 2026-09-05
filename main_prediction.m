%-----
%      New empirical RWR-derived model for the probability of surface rupture (PSR)
%      RWR:        Rupture-width ratio, defined as W_Rup divided by W_Flt.
%      W_Rup:      Rupture width.
%      W_Flt:      Fault width, considering fault geometry such as the fault dip angle
%                    and seismogenic thickness.
%      Rake:       Fault rake angles, for specify style of faulting (SoF).
%                                              by Bob J.Y. Huang in September 2026
%-----
assummed_RWR_case1=0:0.01:1;
assummed_Rake_case1=-999; % for using S1 global model, does not consider SoF
for i=1:length(assummed_RWR_case1)
  [PSR_pred_case1(i),PSR_lo_case1(i),PSR_hi_case1(i)]=sub_RWR_derived_PSRmodel(assummed_RWR_case1(i),assummed_Rake_case1);
end
figure(1);
H1=subplot(1,1,1);
  plot(assummed_RWR_case1,PSR_pred_case1,'k- ','linewidth',4);hold on;
  plot(assummed_RWR_case1,PSR_lo_case1,'k- ','linewidth',1);
  plot(assummed_RWR_case1,PSR_hi_case1,'k- ','linewidth',1);
  xlabel('RWR (W_R_u_p/W_F_l_t)','FontSize',18,'fontweight','bold');
  ylabel('PSR','FontSize',18,'fontweight','bold');
  legend('S1 global model',2);
  grid on;set(H1,'FontSize',15,'fontweight','bold','linewidth',2);
  outnm='RWR_PSR_S1_global_model';
  print(gcf,outnm,'-deps2c');unix(['convert -density 200 ',outnm,'.eps ',outnm,'.png']);

assummed_RWR_case2=0:0.01:1;
assummed_Rake_case2=0; % for using S2 SoF model, SS
for i=1:length(assummed_RWR_case2)
  [PSR_pred_case2(i),PSR_lo_case2(i),PSR_hi_case2(i)]=sub_RWR_derived_PSRmodel(assummed_RWR_case2(i),assummed_Rake_case2);
end
assummed_RWR_case3=0:0.01:1;
assummed_Rake_case3=90; % for using S2 SoF model, RV
for i=1:length(assummed_RWR_case3)
  [PSR_pred_case3(i),PSR_lo_case3(i),PSR_hi_case3(i)]=sub_RWR_derived_PSRmodel(assummed_RWR_case3(i),assummed_Rake_case3);
end
assummed_RWR_case4=0:0.01:1;
assummed_Rake_case4=-90; % for using S2 SoF model, NML
for i=1:length(assummed_RWR_case4)
  [PSR_pred_case4(i),PSR_lo_case4(i),PSR_hi_case4(i)]=sub_RWR_derived_PSRmodel(assummed_RWR_case4(i),assummed_Rake_case4);
end

figure(2);
H2=subplot(1,1,1);
  plot(assummed_RWR_case2,PSR_pred_case2,'b- ','linewidth',4);hold on;
  plot(assummed_RWR_case3,PSR_pred_case3,'m- ','linewidth',4);
  plot(assummed_RWR_case4,PSR_pred_case4,'g- ','linewidth',4);
  plot(assummed_RWR_case2,PSR_lo_case2,'b- ','linewidth',1);
  plot(assummed_RWR_case2,PSR_hi_case2,'b- ','linewidth',1);
  plot(assummed_RWR_case3,PSR_lo_case3,'m- ','linewidth',1);
  plot(assummed_RWR_case3,PSR_hi_case3,'m- ','linewidth',1);
  plot(assummed_RWR_case4,PSR_lo_case4,'g- ','linewidth',1);
  plot(assummed_RWR_case4,PSR_hi_case4,'g- ','linewidth',1);
  xlabel('RWR (W_R_u_p/W_F_l_t)','FontSize',18,'fontweight','bold');
  ylabel('PSR','FontSize',18,'fontweight','bold');
  legend('S2 SoF model (SS)','S2 SoF model (RV)','S2 SoF model (SoF)',2);
  grid on;set(H2,'FontSize',15,'fontweight','bold','linewidth',2);
  outnm='RWR_PSR_S2_prediction_examples';
  print(gcf,outnm,'-deps2c');unix(['convert -density 200 ',outnm,'.eps ',outnm,'.png']);


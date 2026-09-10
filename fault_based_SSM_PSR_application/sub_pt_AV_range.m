function [nu]=sub_pt_AV_range(target_M,target_AV_range,colorcode,bar_long);
%                                     by Bob J.Y. Huang in March 2026
nu='';
plot([target_M,target_M],[target_AV_range(1),target_AV_range(2)],'linestyle','-','color',colorcode,'linewidth',1);
plot([target_M-bar_long,target_M+bar_long],[target_AV_range(1),target_AV_range(1)],'linestyle','-','color',colorcode,'linewidth',1);
plot([target_M-bar_long,target_M+bar_long],[target_AV_range(2),target_AV_range(2)],'linestyle','-','color',colorcode,'linewidth',1);


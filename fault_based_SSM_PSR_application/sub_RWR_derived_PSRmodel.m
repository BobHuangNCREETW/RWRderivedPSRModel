function [PSR_pred, PSR_EU_lo, PSR_EU_hi] = sub_RWR_derived_PSRmodel(RWR, NZbar_Hypo)
%-------------------------------------------------------
%      New empirical RWR-derived model for the probability of surface rupture (PSR),
%        derived from the Beta cumulative distribution function (CDF). 
%      *** Inputs: ***
%      RWR:              Rupture-width ratio, calculated as the rupture width divided 
%                          by the fault width.
%      NZbar_Hypo:       The regional (fault-specific) normalized mean hypocenter depth.
%      ***
%      *** Outputs: ***
%      PSR_pred:         The mean prediction PSR.
%      PSR_EU_[lo]&[hi]: The epistemic uncertainty of the lower and higher bound, 
%                          calculated from 5% to 95% confidence interval (1.96 times
%                          sigma_{epi}). 
%      ***
%      cov_theta_S2_RWR: The covariance matrix from the maximum likelihood estimation
%                          (MLE) for the prediction parameters c1 and c1.
%      alpha:            The shape parameter for the Beta-CDF, which is replaced by 
%                          c0 and c1 in Eq. 8
%      beta:             The second shape parameter for the Beta-CDF, which fixed to 1 
%                          (rounded from coefficient derived from all dataset in 
%                          Table 1).
%      c0, c1:           Coefficients obtained from functional form in Eq. 8, calculated
%                          from the MLE.
%                                                 by Bob J.Y. Huang on March 2026
%-----
%--- The parameters used are derived on all earthquake data in Table 1, and are listed in Table 2
c0=1.5638;
c1=1.8159;
cov_mat_c0c1=[0.0201,0.0566;0.0566,3.1077];
beta = 1; % rounding from the all model (0.9914)
%---
% alpha calculated from Eq. 8
NZbar_Hypo_c=NZbar_Hypo-0.4733; % center at mean normalized Z[bar]Hypo (=0.4733)
alpha = c0 + c1*NZbar_Hypo_c;
%
if(RWR<0|RWR>1)
  error('Invalid input for RWR, which should be within 0 and 1.');
end
if(alpha<=0)
  error('Predicted alpha <= 0; recheck the input NZ[bar]_{Hypo} value');
end

% mean prediction
PSR_pred(:,1) = betacdf(RWR, alpha, beta);

% derivative wrt alpha (numerical)
delta = 1e-5;
df_dalpha = zeros(length(RWR),1);

for i = 1:length(RWR)
    f0 = betacdf(RWR(i), alpha, beta);
    f1 = betacdf(RWR(i), alpha + delta, beta);
    df_dalpha(i) = (f1 - f0)/delta;
end

% gradient wrt [c1,c2]
grad_c = [df_dalpha, df_dalpha*NZbar_Hypo_c];   % Nx2

% variance propagation
var_pred = sum((grad_c * cov_mat_c0c1) .* grad_c, 2);
std_pred(:,1) = sqrt(var_pred);

% 95% CI
zval = 1.96;
PSR_EU_lo = max(0, PSR_pred - zval*std_pred);
PSR_EU_hi = min(1, PSR_pred + zval*std_pred);

end


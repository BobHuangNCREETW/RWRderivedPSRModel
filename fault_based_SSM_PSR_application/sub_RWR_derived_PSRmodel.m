function [PSR_pred, PSR_EU_lo, PSR_EU_hi] = sub_RWR_derived_PSRmodel(RWR, Rake)
%-------------------------------------------------------
%      New empirical RWR-derived model for the probability of surface rupture (PSR),
%        derived from the Beta cumulative distribution function (CDF).
%      *** Inputs:  ***
%      RWR:              Rupture-width ratio, calculated as the rupture width divided
%                          by the fault width.
%      Rake:             For S2-Model use includes SoF consideration. 
%                        If rake = -999, meaning do not consider SoF and use the S1 
%                          global model by considering RWR term only.
%      ****************
%      *** Outputs: ***
%      PSR_pred:         The mean prediction PSR.
%      PSR_EU_[lo]&[hi]: The epistemic uncertainty of the lower and higher bound,
%                          calculated from the 5% to 95% confidence interval (1.96 times
%                          sigma_{epi}).
%      ****************
%      var_alpha:        The variance of the alpha estimation from Maximum likelihood 
%                          (MLE).
%      se_pred:          The standard error (epistemic uncertainty) of the alpha 
%                          estimation.
%      alpha:            The shape parameter for the Beta-CDF
%      beta:             The second shape parameter for the Beta-CDF, which is fixed to 1
%                          (rounded from coefficient derived from all dataset in
%                          Table S1).
%                                                 by Bob J.Y. Huang in September 2026
%-------------------------------------------------------
%--- The parameters used are derived from all earthquake data in Table S1, and are listed in Table 1
if(Rake==-999)
  alpha = 1.555;
  var_alpha=0.0194;
elseif((Rake>=-180&Rake<-150)|(Rake>=-30&Rake<30)|(Rake>=150&Rake<=180))
  alpha = 1.474;
  var_alpha=0.0497;
elseif(Rake>=30&Rake<150);
  alpha = 2.180;
  var_alpha=0.0847;
elseif(Rake>=-150&Rake<-30)
  alpha = 0.719;
  var_alpha=0.0357;
else
  error('Invalid input for Rake angle, which should be within -180 and 180 degree. Or, use -999 for the S1 global model not specific the SoF.');
end

beta = 1; % rounding from the global dataset (0.9914)

if(RWR<0|RWR>1)
  error('Invalid input for RWR, which should be within 0 and 1.');
end

% Prediction for mean probability
PSR_pred(:,1) = betacdf(RWR, alpha, beta);

% Delta method gradient, derivative wrt alpha (numerical)
delta = 1e-5;
df_dalpha  = zeros(length(RWR), 1);

for i = 1:length(RWR)
    f0 = betacdf(RWR(i), alpha, beta);
    f1 = betacdf(RWR(i), alpha + delta, beta);
    df_dalpha(i) = (f1 - f0)/delta;
end

% Prediction variance
var_pred = df_dalpha.^2 * var_alpha;
se_pred(:,1) = sqrt(var_pred);

% 95% confidence interval
zval = 1.96;
PSR_EU_lo = max(0, PSR_pred - zval*se_pred);
PSR_EU_hi = min(1, PSR_pred + zval*se_pred);

end





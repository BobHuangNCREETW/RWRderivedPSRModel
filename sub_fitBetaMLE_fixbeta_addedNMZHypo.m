function [c_hat, cov_theta, ci] = ...
    sub_fitBetaMLE_fixbeta_addedNMZHypo(x_data, y_data, z_data)
%-------------------------------------------------------
% P(x) = betacdf(x_data, alpha, beta); 
%   RWR is used for x, Surface rupture status (binary data 0 and 1 only) is for y, 
%   Normalized average hypocenter depth (NZbar_Hypo) is for z.
%   alpha > 0, beta > 0
% Fits Beta CDF Bernoulli model with:
%   alpha_i = c1 + c2*z_i
%   beta fixed 
% Uses MLE and Hessian-based covariance
%                                          by Bob J.Y. Huang in March 2026
%-------------------------------------------------------

beta = 1; % rounding from the global model (0.9914)

% Ensure valid RWR domain
min_x_data=min(x_data);max_x_data=max(x_data);
if(min_x_data<0|max_x_data>1)
  warning(['****** invalid RWR input; the RWR should be >=0 and <=1 ******',sprintf('\n'),'         ****** check carefully for the used RWR dataset         ******']);
  c_hat=[-999,-999];nll_fun=' ';cov_theta=[-999,-999;-999,-999];ci=[-999,-999;-999,-999];
  return
end
min_z_data=min(z_data);max_z_data=max(z_data);
if(min_z_data<0|max_z_data>1)
  warning(['****** invalid NZ[bar]Hypo input: the normalized average',sprintf('\n'),'hypocenter depth should be >=0 and <=1 ******',sprintf('\n'),'         ****** check carefully for the used NZ[bar]Hypo dataset         ******']);
  c_hat=[-999,-999];nll_fun=' ';cov_theta=[-999,-999;-999,-999];ci=[-999,-999;-999,-999];
  return
end

% initial guess
ini = [2; 0];   % [c1; c2]

% Centering normalized Z[bar]Hypo (as z_c) to ensure c0=alpha in alpha=c0+c1x(normalized Z[bar]Hypo-z_c)
z_mean=0.4733;
z_c=z_data-z_mean; % centering z

% negative log-likelihood handle
nll_fun = @(p) nll_beta_cdf_NMZHypo_alpha(p, x_data, y_data, z_c, beta);

%opts = optimset('Display','iter','TolX',1e-6,'TolFun',1e-6);
options = optimset('fmincon');
options.Display = 'iter';
options.TolX = 1e-6;
options.TolFun = 1e-6;

% MLE
[c_hat, ~] = fminsearch(nll_fun, ini, options);

% Hessian and covariance
H = estimateHessian(nll_fun, c_hat);
cov_theta = inv(H);

% standard errors and CI
se = sqrt(diag(cov_theta));
ci = [c_hat(:) - 1.96*se, c_hat(:) + 1.96*se];

end


function nll = nll_beta_cdf_NMZHypo_alpha(params, x, y, z_c, beta)

c0 = params(1);
c1 = params(2);

alpha = c0 + c1 .* z_c;

% enforce positivity of alpha
if any(alpha <= 0) || beta <= 0
    nll = Inf;
    return
end

p = betacdf(x, alpha, beta);

% numerical safety
epsVal = 1e-10;
p = min(max(p, epsVal), 1-epsVal);

nll = -sum(y .* log(p) + (1-y) .* log(1-p));

end

function H = estimateHessian(fun, theta)
    n = length(theta);
    h = 1e-5;
    H = zeros(n);
    for i = 1:n
        for j = 1:n
            e_i = zeros(n,1); e_i(i) = 1;
            e_j = zeros(n,1); e_j(j) = 1;
            f1 = fun(theta + h*e_i + h*e_j);
            f2 = fun(theta + h*e_i - h*e_j);
            f3 = fun(theta - h*e_i + h*e_j);
            f4 = fun(theta - h*e_i - h*e_j);
            H(i,j) = (f1 - f2 - f3 + f4) / (4*h^2);
        end
    end
end



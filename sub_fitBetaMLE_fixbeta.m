function [alpha_fit,var_alpha,ci] = sub_fitBetaMLE_fixbeta(x_data,y_data)
%-----
% Fits Beta CDF via Bernoulli Maximum Likelihood Estimation (MLE)
%   P(x) = betacdf(x, alpha, beta); RWR is used for x here
% alpha > 0, beta > 0
% Confidence intervals from Hessian (asymptotic MLE)
%                                       by Bob J.Y. Huang in March 2026
%-----

% Initial guess
ini_param = 2; % alpha; use fixed beta from global model

% Ensure valid RWR domain
min_x_data=min(x_data);max_x_data=max(x_data);
if(min_x_data<0|max_x_data>1)
  warning(['****** invalid RWR input; the RWR should be >=0 and <=1 ******',sprintf('\n'),'         ****** check carefully for the used RWR dataset         ******']);
  alpha_fit=[-999];nll_fun=' ';var_alpha=[-999];ci=[-999,-999];
  return
end

% Optimization options 
options = optimset('fmincon');
options.Display = 'iter';
options.TolX = 1e-6;
options.TolFun = 1e-6;

% Negative log-likelihood
nll_fun = @(params) nll_beta_cdf(params, x_data, y_data);

% Optimization
    [alpha_fit, ~] = fminsearch(nll_fun, ini_param, options);

% Hessian & covariance
H = estimateHessian(nll_fun, alpha_fit);
var_alpha = 1/H;

se = sqrt(var_alpha);
ci = [alpha_fit - 1.96*se, alpha_fit + 1.96*se];

% Enforce positivity in CI
ci = max(ci, 1e-3);

end

function nll = nll_beta_cdf(params, x, y)

alpha = params(1);
beta = 1; % rounding from the global model (0.9914)

% Enforce positivity
if alpha <= 0 || beta <= 0
    nll = Inf;
    return;
end


% Beta CDF
f = betacdf(x, alpha, beta);

% Numerical safety
epsVal = 1e-10;
f = min(max(f, epsVal), 1 - epsVal);

% Bernoulli negative log-likelihood
nll = -sum(y .* log(f) + (1 - y) .* log(1 - f));

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





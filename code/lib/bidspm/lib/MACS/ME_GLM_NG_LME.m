function [LME] = ME_GLM_NG_LME(P, L0, a0, b0, Ln, an, bn)
% _
% Log Model Evidence for General Linear Model with Normal-Gamma Priors
% FORMAT [LME] = ME_GLM_NG_LME(P, L0, a0, b0, Ln, an, bn)
% 
%     P   - an n x n precision matrix embodying covariance assumptions
%     L0  - a  p x p matrix (prior precision for regression coefficients)
%     a0  - a  1 x 1 scalar (prior shape for residual variance)
%     b0  - a  1 x v vector (prior rate for residual variance)
%     Ln  - a  p x p matrix (posterior precision for regression coefficients)
%     an  - a  1 x 1 scalar (posterior shape for residual variance)
%     bn  - a  1 x v vector (posterior rate for residual variance)
% 
%     LME - a  1 x v vector of log model evidences
% 
% FORMAT [LME] = ME_GLM_NG_LME(P, L0, a0, b0, Ln, an, bn) returns the log
% model evidence [1] for a general linear model with precision matrix P and
% normal-gamma distributed priors / posteriors for regression coefficients
% (L0 / Ln) and residual variance (a0, b0 / an, bn).
% 
% Further information:
%     help ME_GLM_NG
%     help ME_GLM_NG_AnC
%
% References:
% [1] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469-489, eq. 9.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 07/11/2014, 13:10 (V0.2/V7)
%  Last edit: 26/04/2019, 17:55 (V1.4/V20)


% Get model dimensions
%-------------------------------------------------------------------------%
v = size(bn,2);                 % number of time series
n = size(P, 1);                 % number of data points
p = size(Ln,1);                 % number of parameters

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ME_GLM_NG_LME: estimate');

% Calculate log model evidence
%-------------------------------------------------------------------------%
LME = 1/2*MD_mvn_logdet(P,true)  - n/2*log(2*pi) + ...
      1/2*MD_mvn_logdet(L0,true) - 1/2*MD_mvn_logdet(Ln,true) + ...
      gammaln(an)                - gammaln(a0) + ...
      a0*log(b0)                 - an*log(bn);
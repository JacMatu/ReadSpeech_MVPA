
function [Pool,Part,QV] = check_pooled_error(S);

% Pedagogical script to illustrate assumptions about error term(s) in 
% repeated-measures multi-factorial ANOVAs performed on many observations 
% (eg voxels in fMRI). Requires SPM8+ on path.
%
% Based on Henson & Penny (2003) Tech Report ("H&P03"): 
% http://www.mrc-cbu.cam.ac.uk/people/rik.henson/personal/HensonPenny_ANOVA_03.pdf
%
% ...and used for figures in Henson (in press) "ANOVA". Brain Mapping: A
% comprehensive reference.
%
% rik.henson@mrc-cbu.cam.ac.uk, Jan 2011
%
% Estimates p-values for a number of contrasts (ANOVA effects) as a
% function of a (small) number of subjects and a (large) number of voxels
% (experiments), either pooling the error over all conditions (effects), 
% or partitioning the error for each effect. Data can be generated 
% according to either a single or partitioned error, which is either white 
% (IID/spherical) or coloured (non-IID; nonspherical)
%
% Outputs:
%    Pool = subjects X effects X voxel matrix of p-values, assuming a
%           pooled error
%    Part = subjects X effects X voxel matrix of p-values, assuming
%           partitioned errors
%    AQV  = estimated error covariance of pooled model (from highest 
%           number of subjects provided)
%
% Also produces some figures - explained below!
%
% Written for ease of understanding, in relation to what SPM does, rather
% than speed or elegance!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Demonstration 1: If the data are generated by a single IID error 
% (albeit unlikely), then pooling is more efficient than partitioning error
%
% Run by:
%
%     S.ftit = 'Single white error';           % Just title for plot
%     [Pool,Part,QV] = check_pooled_error(S);     
%
%     ie using defaults below, ie a 2x2 ANOVA design with factors A and B
%
% Results: 
%
% 1. When there is no true effect (B=0 for main effect of B; red lines), 
% approx 5 percent of voxels survive p<.05 in both types of model - ie 
% type I error controlled with both pooled and partitioned error
%
% 2. When there is an effect however (eg, big main effect of A, or smaller
% interaction), a greater percentage of voxels survive p<.05 for pooled
% than partitioned error model (ie more power, because greater df's with 
% which to estimate error)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Demonstration 2: If the data are generated by partitioned but IID error 
% terms, such that fixed correlation between errors induced, then pooling 
% is still more efficient, with no bias (unless those error terms are
% non-IID, eg different (nonspherical) errors, as in demos below).
%
% Run by:
%   
%     S.GenPart = 1;
%     S.ftit = 'Partitioned errors';          
%     [Pool,Part,QV] = check_pooled_error(S);     
%
% Results: 
%
% 1. Similar to above, noting that when there is no true effect 
% (B=0 for main effect of B; red lines), pooled model still unbiased
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Demonstration 3: If there is greater correlation in the error (nonsphericity), 
% owing to similar scores across conditions from the same subject (which is 
% likely), then pooled error becomes biased, but partitioning error is less 
% biased (indeed, UNbiased when only one numerator df per effect, as here)
%
% Run by:
%   
%     S.GenPart = 0;
%     S.V = [1 0 0.5 0; 0 1 0 0; 0.5 0 1 0; 0 0 0 1];
%     S.ftit = 'Correlated error';           
%     % Eg because only true effects of A and AxB, make these correlated
%     [Pool,Part,AQV] = check_pooled_error(S);     
%
% Results: 
%
% 1. Now note that >5% of voxels have p<.05 under pooled error model (thin
% red line) - ie invalid - but still correct 5% for partitioned error model 
% (thick red line) - ie still valid.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Demonstration 4: If you try to estimate error covariance for each voxel
% separately (using ReML), that estimate is poor (unless large number of 
% subjects), resulting in "poor" p-values.
%
% Run by:
%   
%     S.GenPart = 0;
%     S.PreWhiten = 1;                      % Separate ReML for each voxel
%     S.V = [1 0 0.5 0; 0 1 0 0; 0.5 0 1 0; 0 0 0 1];
%     S.ftit = 'Voxel-wise estimation of correlated error';  
%     [Pool,Part,AQV] = check_pooled_error(S);     
%
% Results: 
%
% 1. Warnings as ReML unable to estimate error covariance for some random
% data (some voxels), particularly when few subjects
% 
% 2. Note that >5% of voxels have p<.05 under pooled error model (thin
% red line) - ie invalid.
%
% 3. Note that power to detect real effects also reduced for pooled versus
% partitioned error (green and blue lines lowered).
%
% This danger is probably why partitioning error is preferred in
% behavioural statistics when typically only one variate (one voxel) - any
% method (eg Greenhouse-Geisser) for correcting for correlated errors, 
% based on one sample of data, will entail an inefficient estimate of that 
% correlation (with few df's), so a poor correction (eg too conservative
% in case of G-G).
%
% (This is the same general problem that estimating covariance matrices
% efficiently requires large numbers of observations!)
%
% In imaging data however, we have many samples (voxels), so if we assume
% that error correlation across conditions is same across voxels, we can 
% pool the samples to get a more precise estimate of that error 
% correlation (ie effectively increase number of observations so as to 
% increase accuracy of estimate of error covariance).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Demonstration 5: If you try to estimate the error covariance by pooling
% data across voxels, in a two-step procedure like that used in SPM, then
% you recover greater efficiency of pooled than partitioned model (when
% there is a true effect), together with valid results when there is no
% effect. Note however that this is only true if the underlying error
% correlation IS the same across voxels (see next Demo)
%
% Run by:
%   
%     S.GenPart = 0;
%     S.PreWhiten = 2;            % One ReML from data pooled across voxels
%     S.V = [1 0 0.5 0; 0 1 0 0; 0.5 0 1 0; 0 0 0 1];
%     S.ftit = 'Voxel-wide estimation of correlated error';  
%     [Pool,Part,AQV] = check_pooled_error(S);     
%
% Results: 
%
% 1. We recover the situation in Demonstration 1, where pooling error is
% more powerful when effects exist, but still valid when they do not.
%
% 2. [ To visualise nonsphericity: figure,imagesc(AQV),colormap('gray') ]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Demonstration 6: If you try to estimate the error covariance by pooling
% data across voxels, in a two-step procedure like that used in SPM, BUT
% the true error correlation structure varies across voxels, the type I 
% error is again inflated (though not by much, at least when elements of
% error covariance vary randomly between 0 and N^2, where N is maximum
% element of common covariance term above
%
% Run by:
%   
%     S.GenPart = 0;
%     S.PreWhiten = 2;            % One ReML from data pooled across voxels
%     S,VarVox = 1;
%     S.V = [1 0 0.5 0; 0 1 0 0; 0.5 0 1 0; 0 0 0 1];
%     S.ftit = 'Voxel-wide estimation of correlated but nonstationary error';  
%     [Pool,Part,AQV] = check_pooled_error(S);     
%
% 1. Thin red line shows increased false positive rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    f    = S.f;        % Factorial design
    flab = S.flab;     % Factor labels for plotting
catch
    % Unless otherwise specified, assume a 2x2 design with factors A/B:
    
    f = [1 1 0 0;   % Main effect of A
         1 0 1 0;   % Main effect of B
         1 0 0 1;   % AxB Interaction
         1 1 1 1]'; % Constant
     
    flab = {'meA'; 'meB'; 'int'};   % Factor labels for plotting
end
C = detrend(f,0);    % Convert to contrasts

try 
    B = S.B;         % True Betas for 4 effects above
catch
    % Unless otherwise specified, assume main effect of A, no main effect 
    % of B, weaker interaction (and zero grand mean)
   
    B = [1 0 0.5 0]';    
end
CB = C*B;            % Convert to Betas for individual conditions

try 
    V = S.V;         % Covariance of true error
catch
    % Unless otherwise specified, assume IID error. Note that construction
    % of partitioned error below can also induce nonspherical errors
    V = speye(size(C,2));    
end

try
    subs = S.subs;   % Range of subject numbers to explore
catch
    subs = [8:2:24];
end

try
    Nvox = S.Nvox;      % Number of voxels to explore
catch
    Nvox = 1000;
end

try
    VarVox = S.VarVox;      % Whether true correlation constant across voxels
catch
    VarVox = 0;
end

% Explained above...
try  GenPart   = S.GenPart;       catch  GenPart   = 0;         end
try  PreWhiten = S.PreWhiten;     catch  PreWhiten = 0;         end
% Just figure title for plotting
try  ftit      = S.ftit;          catch  ftit      = [];        end

rand('state',0);

Ngrp = length(subs);
Ncon = size(f,1);
Neff = size(C,2) - 1;  % Ignore constant term (grand mean) 

Pool = zeros(Ngrp,Neff);      % Suprathreshold P-values for Pooled Error
Part = Pool;                  % Suprathreshold P-values for Partitioned Error

QV=[];

for s = 1:Ngrp;
    
    Nsub  = subs(s);
    N     = Nsub * Ncon;
    oS    = ones(Nsub,1);
    oC    = ones(Ncon,1);
    zS    = zeros(Nsub,1);
    zC    = zeros(Ncon,1);
    iS    = eye(Nsub);
    iC    = eye(Ncon);
    
    Bs    = randn(Nsub,1);    % Subject effects

    if PreWhiten > 0
        YY    = zeros(N);
        Q     = varcomps(Ncon,Nsub);
    end
     
    % Basic (familiar) design matrix
    X  = [kron(iC,oS) kron(oC,iS)];         %figure,imagesc(X),colormap('gray'); 
     
    fprintf('%d',Nsub);
    for v = 1:Nvox
        
        % Generic way of generating data (with subject effects)
        % y = mvnrnd(repmat(CB',Nsub,1),V) + repmat(mvnrnd(zS',iS),Ncon,1)'; 
        % y = y(:);
        
        if ~GenPart         % Generate data assuming single error term
            
            if VarVox
                 vV = rand(4) * max(V(:)); vV=vV'*vV;
            else
                 vV = V;
            end
            
            try
                e = mvnrnd(repmat(zC',Nsub,1),vV);  % If have Matlab's stats toolbox...
            catch
                e = spm_normrnd(zC',vV,Nsub)';      % ...else use SPM
            end
            
            y  = X*[CB; Bs] + e(:);
            
        else                % Generate data assuming partitioned errors
            
            % Full design matrix for partitioned model (Fig 9a H&P03)
            Xf = [kron(f,oS) kron(f,iS)];   %figure,imagesc(Xf),colormap('gray');
            e  = mvnrnd(repmat(zC',Nsub,1),V);
            y  = Xf*[B; e(:)];
        end
        
                   
        % Estimate using pooled error, with different types of pre-whitening        
        if PreWhiten == 0
            for e = 1:Neff
                Pool(s,e,v) = Fstat(X,y,C(:,e));
            end
        elseif PreWhiten == 1
            QV    = rik_reml(y*y',X,Q);
            W     = PreWfilter(QV);
            WX    = W*X;
            for e = 1:Neff
                Pool(s,e,v) = Fstat(WX,y,C(:,e));
            end
        elseif PreWhiten == 2
            Ay{s}(:,v) = y';
            YY         = YY + y*y';
        end
        
        % Estimate using partitioned errors (Section 5.1 of H&P03)     
        for e = 1:Neff              
            yc          = kron(C(:,e),iS)'*y; 
            Part(s,e,v) = Fstat(oS,yc,[1]);

            %faster...?
            %T = mean(y)*sN/std(yc);
            %p = (1 - spm_Tcdf(T,N-1))*2;
        end
    end
    fprintf('.')
  
    if PreWhiten == 2        % Estimate nonsphericity pooled across voxels
 
        QV    = rik_reml(YY/Nvox,X,Q);
        W     = PreWfilter(QV);
        WX    = W*X;
        
        for v = 1:Nvox

            y  = W * Ay{s}(:,v);  % Pre-whiten (same) data
            
            % Estimate using pooled error
            for e = 1:Neff 
                Pool(s,e,v) = Fstat(WX,y,C(:,e));         
            end
        end
    end
    fprintf('.')
end
fprintf('\n')


cols = {'b-';'r-';'g-'};
figure,hold on
for e=1:Neff
    for s=1:Ngrp
        pPool(s) = length(find(Pool(s,e,:)<.05))/Nvox;
        pPart(s) = length(find(Part(s,e,:)<.05))/Nvox;
    end
    plot(subs,pPool,cols{e},'LineWidth',1);
    plot(subs,pPart,cols{e},'LineWidth',2);
end
ylabel('Prop. Voxels with p<.05');
xlabel('Number of Subjects');
if ~isempty(flab)
    llab = {};
    for e=1:Neff
        llab{end+1} = [flab{e} 'pool'];
        llab{end+1} = [flab{e} 'part'];
    end
    legend(llab,'Location','Best')
end
if ~isempty(ftit)
    title(ftit)
end
CovEst = full(V-diag(diag(V))); CovEst = max(CovEst(:));  % Silly filename label!
eval(sprintf('print -f%d -djpeg75 GenPart%d_PreWhiten%d_Cov%3.2f.jpg',gcf,GenPart,PreWhiten,CovEst));

return

%-----------------------------------------------------------------


function [p,r,B] = Fstat(X,y,c);

% F-statistic for GLM X, data y and F-contrast matrix c.
% See Appendix A.1 of:
% http://www.mrc-cbu.cam.ac.uk/people/rik.henson/personal/HensonPenny_ANOVA_03.pdf    
%------------------------------------------------------------------

  if size(c,1) < size(X,2)
      c = [c; zeros(size(X,2)-size(c,1),1)];
  end
 
  B = pinv(X)*y;
  Y = X*B;
  r = y - Y;
  
  c_0 = eye(size(X,2)) - c*pinv(c);
  X_0 = X*c_0;
  R   = eye(size(X,1)) - X*pinv(X);
  R_0 = eye(size(X,1)) - X_0*pinv(X_0);
  M = R_0 - R;
  df = [rank(X)-rank(X_0) size(X,1)-rank(X)];
  F  = ((B'*X'*M*X*B)/df(1)) / ((y'*R*y)/df(2)); 
  p  = 1-spm_Fcdf(F,df);

return



function Q = varcomps(Ncon,Nsub)
    
% make set of variance components
%------------------------------------------------------------------
  Q={}; nv=0; z=zeros(Ncon*Nsub); id=[1:Nsub];
  
  for c1 = 1:Ncon
        nv = nv+1;
        v = z;
        v((c1-1)*Nsub + id, (c1-1)*Nsub + id)=eye(Nsub);
        Q{nv} = sparse(v);
  end

  for c1 = 1:Ncon
       for c2 = (c1+1):Ncon
            nv = nv+1;
            v = z;
            v( (c1-1)*Nsub + id, (c2-1)*Nsub + id )=eye(Nsub);
            v( (c2-1)*Nsub + id, (c1-1)*Nsub + id )=eye(Nsub);
            Q{nv} = sparse(v);   
        end
  end

return
  

    
function W = PreWfilter(V);

% make W a whitening filter W*W' = inv(V)
%------------------------------------------------------------------
   [u s] = spm_svd(V);
   s     = spdiags(1./sqrt(diag(s)),0,length(s),length(s));
   W     = u*s*u';
   W     = sparse(W.*(abs(W) > 1e-6));

return



function [V,h,Ph,F,Fa,Fc] = rik_reml(YY,X,Q,N,D,t)

% estimate covariance components (from spm_reml)
%------------------------------------------------------------------

% ReML estimation of [improper] covariance components from y*y'
% FORMAT [C,h,Ph,F,Fa,Fc] = spm_reml(YY,X,Q,N,D,t);
%
% YY  - (m x m) sample covariance matrix Y*Y'  {Y = (m x N) data matrix}
% X   - (m x p) design matrix
% Q   - {1 x q} covariance components
% N   - number of samples
% D   - Flag for positive-definite scheme
% t   - regularisation (default 4)
%
% C   - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) ReML hyperparameters h
% Ph  - (q x q) conditional precision of h
%
% F   - [-ve] free energy F = log evidence = p(Y|X,Q) = ReML objective
%
% Fa  - accuracy
% Fc  - complexity (F = Fa - Fc)
%
% Performs a Fisher-Scoring ascent on F to find ReML variance parameter
% estimates.
%
%__________________________________________________________________________
%
% John Ashburner & Karl Friston
% $Id: spm_reml.m 3791 2010-03-19 17:52:12Z karl $
 
 
% check defaults
%--------------------------------------------------------------------------
try, N; catch, N  = 1;  end       % assume a single sample if not specified
try, K; catch, K  = 32; end       % default number of iterations
try, D; catch, D  = 0;  end       % default checking
try, t; catch, t  = 4;  end       % default regularisation
 
% catch NaNs
%--------------------------------------------------------------------------
W     = Q;
q     = find(all(isfinite(YY)));
YY    = YY(q,q);
for i = 1:length(Q)
    Q{i} = Q{i}(q,q);
end
 
% dimensions
%--------------------------------------------------------------------------
n     = length(Q{1});
m     = length(Q);
 
% ortho-normalise X
%--------------------------------------------------------------------------
if isempty(X)
    X = sparse(n,0);
else
    X = spm_svd(X(q,:),0);
end
 
% initialise h and specify hyperpriors
%==========================================================================
h   = zeros(m,1);
for i = 1:m
    h(i,1) = any(diag(Q{i}));
end
hE  = sparse(m,1);
hP  = speye(m,m)/exp(32);
dF  = Inf;
D   = 8*(D > 0);
 
 
% ReML (EM/VB)
%--------------------------------------------------------------------------
for k = 1:K

    % compute current estimate of covariance
    %----------------------------------------------------------------------
    C     = sparse(n,n);
    for i = 1:m
        C = C + Q{i}*h(i);
    end
 
    % positive [semi]-definite check
    %----------------------------------------------------------------------
    for i = 1:D
        if min(real(eig(full(C)))) < 0

            % increase regularisation and re-evaluate C
            %--------------------------------------------------------------
            t     = t - 1;
            h     = h - dh;
            dh    = spm_dx(dFdhh,dFdh,{t});
            h     = h + dh;
            C     = sparse(n,n);
            for i = 1:m
                C = C + Q{i}*h(i);
            end
        else
            break
        end
    end


    % E-step: conditional covariance cov(B|y) {Cq}
    %======================================================================
    iC     = spm_inv(C);
    iCX    = iC*X;
    if ~isempty(X)
        Cq = spm_inv(X'*iCX);
    else
        Cq = sparse(0);
    end

    % M-step: ReML estimate of hyperparameters
    %======================================================================

    % Gradient dF/dh (first derivatives)
    %----------------------------------------------------------------------
    P     = iC - iCX*Cq*iCX';
    U     = speye(n) - P*YY/N;
    for i = 1:m

        % dF/dh = -trace(dF/diC*iC*Q{i}*iC)
        %------------------------------------------------------------------
        PQ{i}     = P*Q{i};
        dFdh(i,1) = -sum(sum(PQ{i}'.*U))*N/2;

    end

    % Expected curvature E{dF/dhh} (second derivatives)
    %----------------------------------------------------------------------
    for i = 1:m
        for j = i:m

            % dF/dhh = -trace{P*Q{i}*P*Q{j}}
            %--------------------------------------------------------------
            dFdhh(i,j) = -sum(sum(PQ{i}'.*PQ{j}))*N/2;
            dFdhh(j,i) =  dFdhh(i,j);

        end
    end
 
    % add hyperpriors
    %----------------------------------------------------------------------
    e     = h     - hE;
    dFdh  = dFdh  - hP*e;
    dFdhh = dFdhh - hP;

    % Fisher scoring: update dh = -inv(ddF/dhh)*dF/dh
    %----------------------------------------------------------------------
    dh    = spm_dx(dFdhh,dFdh,{t});
    h     = h + dh;

    % predicted change in F - increase regularisation if increasing
    %----------------------------------------------------------------------
    pF    = dFdh'*dh;
    if pF > dF
        t = t - 1;
    else
        t = t + 1/4;
    end
    
    % revert to SPD checking, if near phase-transition
    %----------------------------------------------------------------------
    if abs(pF) > 1e6
        [V,h,Ph,F,Fa,Fc] = rik_reml(YY,X,Q,N,1,t - 2);
        return
    else
        dF = pF;
    end
    
    % Convergence (1% change in log-evidence)
    %======================================================================
    %fprintf('%s %-23d: %10s%e [%+3.2f]\n','  ReML Iteration',k,'...',full(pF),t);
 
    % final estimate of covariance (with missing data points)
    %----------------------------------------------------------------------
    if dF < 1e-1, break, end
end

 
% re-build predicted covariance
%==========================================================================
V     = 0;
for i = 1:m
    V = V + W{i}*h(i);
end
 
% check V is positive semi-definite (if not already checked)
%==========================================================================
if ~D
    if min(eig(V)) < 0
        [V,h,Ph,F,Fa,Fc] = rik_reml(YY,X,Q,N,1,2);
        return
    end
end
 
% log evidence = ln p(y|X,Q) = ReML objective = F = trace(R'*iC*R*YY)/2 ...
%--------------------------------------------------------------------------
Ph    = -dFdhh;
if nargout > 3
 
    % tr(hP*inv(Ph)) - nh + tr(pP*inv(Pp)) - np (pP = 0)
    %----------------------------------------------------------------------
    Ft = trace(hP*inv(Ph)) - length(Ph) - length(Cq);
 
    % complexity - KL(Ph,hP)
    %----------------------------------------------------------------------
    Fc = Ft/2 + e'*hP*e/2 + spm_logdet(Ph*inv(hP))/2 - N*spm_logdet(Cq)/2;
 
    % Accuracy - ln p(Y|h)
    %----------------------------------------------------------------------
    Fa = Ft/2 - trace(C*P*YY*P)/2 - N*n*log(2*pi)/2 - N*spm_logdet(C)/2;
 
    % Free-energy
    %----------------------------------------------------------------------
    F  = Fa - Fc;
 
end

   

function CircleProg

%%%%%%%%%%%%%%%%%%%%%%%%%%%    OPTIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DOACCURATE = 1==1;

% 0 = test case
% 1 = high frequency
% 2 = low frequency
problemOption = 2;

%$$ nPts = 120;

if DOACCURATE
  if problemOption == 2 % low frequency
    nPts = 80;
  elseif problemOption == 0
    nPts = 60;
  else
    nPts = 60;
  end
else
  if problemOption == 2 % low frequency
    nPts = 80;  %70; % 80;
  else
    nPts = 60;
  end
end

% finite difference order
fdOrder = 7;

% which states to compute overlaps for
if problemOption == 0
  iSel = [1:4, 6,7, 10,11];
else
  iSel = 1:15 ;
end
iStart = 1;   % initial state

%  0 = one run
%  1 = optimization
% -1 = domcke method
WHICH_DO = 0;

% DOMCKE METHOD, and optimization
DOAVERAGE  = false;
nKp        = 3;  %4;

% DOMCKE METHOD
nPhase     = 3; %7;
DOFTLAST   = 1==1;
DOFTK      = 1==0;
DOCOMPLEXD = 1;

if WHICH_DO == 1
  longerFac = 1.1;
else
  % propagate for extra time at the end
  %
  %longerFac = 6.8541;
  %longerFac = 4.2361;
  %$$ longerFac = 2.618;
  % longerFac = 1.1;
  
  longerFac = 1.618;
end

% %$$ timeRes = 1;
% %$$ modPlot = 20;
% % timeRes = 0.5;
% % modPlot = 40;
% %
% timeRes = 0.3;
% modPlot = 70;

if DOACCURATE
  if problemOption == 2
    timeRes = 0.05;     % for slow field..  need this for high frequency
    modPlot = 200;
  elseif problemOption == 1
    timeRes = 0.1;  %0.12;  % 0.2;
    modPlot = 100;
  else
    timeRes = 0.06; %0.1;
    modPlot = 100;
  end
else
  if problemOption == 2
    timeRes = 0.07;  %0.1;     % for slow field
    modPlot = 100;
  else
    timeRes = 0.5;
    modPlot = 20;
  end
end
clear DOACCURATE

propMode = 1;         % 0 = Crank-Nicholson otherwise expm

oStr = 1e0;           % pulse strength multiplication factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[vMax,omega,duration,strXS,strXC,strYS,strYC,ENV_OPT] = GetPulseParams(problemOption);

switch WHICH_DO
  case -1   % DOMCKE METHOD
    CircleDomcke(DOAVERAGE, nPhase, nKp, ...
      problemOption,nPts,fdOrder,vMax,omega,duration,...
      strXS,strXC,strYS,strYC,oStr,ENV_OPT, ...
      iStart,iSel,        longerFac,timeRes,propMode,modPlot, ...
      DOFTLAST, DOFTK, DOCOMPLEXD);
  case 0    % JUST ONE RUN
    iP = 1;  iPh = 2*pi*(iP-1)/nPhase;
    jP = 1;  jPh = 2*pi*(jP-1)/nPhase;
    kP = 1;  kPh = 2*pi*(kP-1)/nKp;
    CircleProp(problemOption,nPts,fdOrder,vMax,omega,duration,...
      strXS,strXC,strYS,strYC,oStr,iPh,jPh,kPh,ENV_OPT,...
      iStart,iSel,   true,longerFac,timeRes,propMode,modPlot);
    disp('...DONE')
  case 1   % OPTIMIZATION
    if problemOption == 0
      iOpt = 3;
    elseif problemOption == 2  % low frequency
      iOpt = 3;  %2;
    else
      error('not supported')
    end
    CircleOpt(DOAVERAGE,nKp, ...
      problemOption,nPts,fdOrder,vMax,omega,duration,...
      strXS,strXC,strYS,strYC,oStr,ENV_OPT,...
      iStart, iOpt, longerFac,timeRes,propMode,modPlot);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%     OPTIMIZATION     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CircleOpt(DOAVERAGE,nKp, ...
  problemOption,nPts,fdOrder,vMax,omega,duration,...
  strXS,strXC,strYS,strYC,oStr,ENV_OPT,...
  iStart, iOpt,  longerFac,timeRes,propMode,modPlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optVars = [2;4;5;6;7];
% optVars = 2;
% optVars = [2,8];
% optVars = [4;5;6;7];
% optVars = [5;6;7];
% optVars = [2;5;6;7];
% optVars = [3;5;6;7];
optVars = [4;5];
% optVars = [2;4;5;6;7];
whichOpt = 1;   % 0 = fminunc otherwise cmaes
stepTolCmaes = 5e-5; %$$ 2e-4;          %5e-4;  %2e-3;  %0.02;
stepTolFminunc = 1e-6;
funTol = 1e-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tempPlot = false;

OptCore  = @(vMax,omega,duration,strXS,strXC,strYS,strYC,oStr,pShift) ...
  CircleProp(problemOption,nPts,fdOrder,vMax,omega,duration,...
  strXS,strXC,strYS,strYC,oStr,pShift,pShift,pShift,ENV_OPT,...
  iStart,iOpt, ...
  tempPlot,longerFac,timeRes,propMode,modPlot);

%%%%%%%%%%

addpath('./CMAES/');

OptFun00 = @(vMax,omega,duration,strXS,strXC,strYS,strYC,oStr) ...
  OptimizationFun(problemOption,DOAVERAGE,nKp,OptCore,vMax,omega,duration,strXS,strXC,strYS,strYC,oStr);

vars  = [vMax;omega;duration;strXS;strXC;strYS;strYC;oStr];

nVars = 8;

%$$ BB    = 0.002;
%$$ BB    = BB * sqrt(mean([strXS,strXC,strYS,strYC].^2));

%$$ stds  = [0.002; 0.002; 0.005; BB*ones(4,1); 0.01]    ;

%AA =0.002;     BB = 0.05;
%AA =0.00002;   BB = 0.0002;
%AA =0.00002;   BB = 0.004;
%
% AA =0.0005;   BB = 0.01;  CC = 0.001;
AA =0.0005;   BB = 0.1;  CC = 0.001;
stds  = [Inf; AA; CC; BB*ones(4,1); 0.01]  ;

if numel(stds)~=nVars
  error('ack')
end

% %$$ varFun = @(x) [12+exp(x(1));exp(x(2:3));x(4:7)];
% %$$ varInv = @(x) [log(x(1)-12);log(x(2:3));x(4:7)];

% constrain 4 < omega < 10
% duration positive
%varFun = @(x) [Inf;10 + 5*sin(x(2)); exp(x(3));x(4:7)];
%varInv = @(x) [Inf;asin((x(2)-10)/5);log(x(3));x(4:7)];

% varFun = @(x) [Inf;x(2); exp(x(3));x(4:8)];
% varInv = @(x) [Inf;x(2); log(x(3));x(4:8)];

%varFun = @(x) [Inf;exp(x(2:3));x(4:8)];
%varInv = @(x) [Inf;log(x(2:3));x(4:8)];
%%%%%%%%%%%

% %$$ varFun = @(x) [Inf;10 + 5*sin(x(2)); exp(x(3));x(4:8)];
% %$$ varInv = @(x) [Inf;asin((x(2)-10)/5);log(x(3));x(4:8)];


%$$ varFun = @(x) [Inf;6 + 3*sin(x(2)); exp(x(3));x(4:8)];
%$$ varInv = @(x) [Inf;asin((x(2)-6)/3);log(x(3));x(4:8)];

varFun = @(x) [Inf;omega + omega/2*sin(x(2)); exp(x(3));x(4:8)];
varInv = @(x) [Inf;asin((x(2)-omega)/omega*2);log(x(3));x(4:8)];

OptFun0 = @(x) OptFun00(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8));

nOpt = numel(optVars);

isOpt = false(1,nVars);
isOpt(optVars) = true;

optMat = full(sparse(optVars,1:nOpt,ones(1,nOpt),nVars,nOpt));

vars = varInv(vars);
varsMask = vars;
varsMask(isOpt) = 0;
OptFun = @(x) OptFun0(varFun(varsMask + optMat * x));

guess = vars(optVars);
stds  = stds(optVars);

error0 = OptFun(guess);
fprintf([ ...
  '***********************************\n' ...
  'Error at beginning of fit: %20.10f \n' ...
  '***********************************\n' ...
  ],error0);

[cmaesOptions,fminuncOptions] = EngineOptions(stepTolCmaes,stepTolFminunc, ...
  funTol*abs(error0));

if whichOpt == 0
  xFit0 = fminunc(OptFun,guess,fminuncOptions);
else
  [~,~,~,~,~,bestever] = cmaes_djh(OptFun,guess,stds,cmaesOptions);
  xFit0 = bestever.x;
  % fFit = bestever.f;
end

xFit = varsMask + optMat * xFit0;

xFit = varFun(xFit);

disp(xFit)

end


function optVal = OptimizationFun(problemOption,DOAVERAGE,nKp,OptCore,vMax,omega,duration,strXS,strXC,strYS,strYC,oStr)
if DOAVERAGE == 0
  [wmax,pmin] = OptimizationFun0(OptCore,vMax,omega,duration,strXS,strXC,strYS,strYC,oStr,0);
else
  phases = (0:nKp-1)/nKp * 2 * pi;
  wmax=0;
  pmin=0;
  for ip = 1:nKp
    [wmaxi,pmini] = OptimizationFun0(OptCore,vMax,omega,duration,strXS,strXC,strYS,strYC,oStr,phases(ip));
    wmax = wmax + wmaxi / nKp;
    pmin = pmin + pmini / nKp;
  end
end

optVal = MyOptVal(wmax,pmin,problemOption);

% fprintf('VARS    ::   vMax %f  om %f  dur %f  str %f  %f  %f  %f  ostr %f\n',...
%  vMax,omega,duration,strXS,strXC,strYS,strYC,oStr)

totStr = sqrt(strXS^2+strXC^2+strYS^2+strYC^2);

fprintf('VARS: om %f dur %f str %f %f %f %f ostr %f totstr %f \n',...
  omega,duration,strXS,strXC,strYS,strYC,oStr,totStr)

%fprintf('VARS: om %f dur %f str %f %f %f %f ostr %f totstr %f \n',...
%  omega,duration,strXS,strXC,strYS,strYC,oStr,totStr)

fprintf('                               OPTVAL  ::   %20.12f \n',optVal);
pause(1)
end


function [wmax,pmin] = OptimizationFun0(OptCore,vMax,omega,duration,strXS,strXC,strYS,strYC,oStr,pShift)
[popVel,popLen,workVel,workLen] = OptCore(vMax,omega,duration,strXS,strXC,strYS,strYC,oStr,pShift);

popVel = real(popVel);
popLen = real(popLen);
workVel = real(workVel);
workLen = real(workLen);

% wmin            = min(workVel,workLen);
wmax            = max(workVel,workLen);

if wmax < 0
  disp('WHAT')
  wmax = eps;
end

pmin            = min(popVel,popLen);

pmin = prod(pmin);

if pmin>=1 || pmin <= 0 || wmax <=0
  pause
  error('ack')
end

end


function optVal = MyOptVal(wmax,pmin,problemOption)

if problemOption == 2 % low frequency
  
  optVal = -log(pmin);
  
elseif problemOption == 0 % test case
  
  optVal = log(wmax)*0.7 - log(pmin);
  
  % optVal = log(wmax)*0.5 - log(pmin);
  optVal = log(wmax)*0.6 - log(pmin);
  
  % optVal = wmax*10 - log(pmin);
  % optVal = wmax*3 - log(pmin);
  optVal = wmax - log(pmin);
  
  %{
%popOpt          = pmin/wmax^0.67;

% popOpt          = pmin/sqrt(wmax);

% popOpt          = pmin/wmax^2;

popOpt          = pmin/wmax;

%%%%%%%

% popOpt          = pmin^0.67/wmax;

% popOpt          = pmin^0.5/wmax;

% popOpt          = pmin^0.35/wmax;

%%%%%%%%
% optVal          = 1./popOpt;
%
%$$ optVal = -popOpt;
%
optVal = -log(max(eps,popOpt));
  %}
  
else
  error('not supported')
end  % problemOption

end  % MyOptVal


function [cmaesOptions,fminuncOptions] = EngineOptions(stepTolCmaes,stepTolFminunc,funTol)
maxEval = 1000;

cmaesOptions = struct( ...
  'MaxFunEvals',maxEval ...
  ,'MaxIter',Inf ...
  ,'TolX',stepTolCmaes ...
  ,'DispModulo',1 ...
  ...%,'Noise',struct('on',1) ...
  ...%,'DiagonalOnly',1 ...
  ...%,'CMA',struct('active',1) ...
  ...%,'PopSize',7 ...
  );

myDisplay = 'iter';
fminuncOptions = optimoptions('fminunc'...
  ,'MaxFunctionEvaluations',maxEval ...
  ,'MaxIterations',Inf ...
  ,'Display',myDisplay ...
  ,'StepTolerance',stepTolFminunc ...
  ,'FunctionTolerance',funTol ...
  ,'OptimalityTolerance',0 ...
  ,'FiniteDifferenceType','central' ...
  ...%$$,'FiniteDifferenceStepSize',1e-4 ...
  ,'FiniteDifferenceStepSize',1e-6 ...
  );

% fminuncOptions = optimoptions('fminunc'...
%   ,'MaxFunctionEvaluations',maxEval ...
%   ,'MaxIterations',Inf ...
%   ,'Display',myDisplay ...
%   ,'StepTolerance',stepTolFminunc ...
%   ,'FunctionTolerance',1e-5 ...
%   ,'FiniteDifferenceType','central' ...
%   ,'FiniteDifferenceStepSize',1e-4 ...
%   );

end

%%%%%%%%%%%%%%%%%%%%%%%   END OPTIMIZATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%



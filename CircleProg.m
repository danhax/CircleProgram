
function CircleProg
%$$ nPts = 120;
% nPts = 80;
nPts = 60;

fdOrder = 7;
% iSel = 1:15 ;
iSel = [1:4, 6,7, 10,11];
iStart = 1;

problemOption = 0;

%longerFac = 6.8541;  

%$$
longerFac = 4.2361;  
%$$longerFac = 2.618;   

% longerFac = 1.0;  

%$$ timeRes = 1;
%$$ modPlot = 20; 
timeRes = 0.5;
modPlot = 40;
% timeRes = 0.3;
% modPlot = 70;

propMode = 1;   % 0 = Crank-Nicholson otherwise expm

[vMax,omega,duration,strXS,strXC,strYS,strYC] = GetPulseParams(problemOption);

WHICH_DO = -1;

oStr = 1e0;  %sqrt(2);

DOAVERAGE = true;
nPhase = 7;  %11;  % 10; %7; %7;
nKp    = 9; %7; %5; %3; %3;

switch WHICH_DO
  case -1
    DOFTLAST = 1==1;
    DOFTK    = 1==0;
    DOCOMPLEXD = 1;
    CircleDomcke(DOAVERAGE, nPhase, nKp, ...
      problemOption,nPts,fdOrder,vMax,omega,duration,...
      strXS,strXC,strYS,strYC,oStr, ...
      iStart,iSel,        longerFac,timeRes,propMode,modPlot, ...
      DOFTLAST, DOFTK, DOCOMPLEXD);
  case 0
    iP = 1;  iPh = 2*pi*(iP-1)/nPhase;
    jP = 1;  jPh = 2*pi*(jP-1)/nPhase;
    kP = 1;  kPh = 2*pi*(kP-1)/nKp;
    CircleProp(problemOption,nPts,fdOrder,vMax,omega,duration,...
      strXS,strXC,strYS,strYC,oStr,iPh,jPh,kPh,...
      iStart,iSel,   true,longerFac,timeRes,propMode,modPlot);
    disp('...DONE')
  case 1
    if problemOption == 0
      iOpt = 3;
    else
      error('not supported')
    end
    CircleOpt(DOAVERAGE,nKp, ...
      problemOption,nPts,fdOrder,vMax,omega,duration,...
      strXS,strXC,strYS,strYC,oStr,...
      iStart, iOpt, longerFac,timeRes,propMode,modPlot);
end

end


function CircleOpt(DOAVERAGE,nKp, ...
  problemOption,nPts,fdOrder,vMax,omega,duration,...
  strXS,strXC,strYS,strYC,oStr,...
  iStart, iOpt,  longerFac,timeRes,propMode,modPlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optVars = [2;4;5;6;7];
% optVars = 2;
% optVars = [2,8];
% optVars = [4;5;6;7];
% optVars = [5;6;7];
% optVars = [2;5;6;7];
optVars = [3;5;6;7];
% optVars = [2;3];
% optVars = [2;4;5;6;7];
whichOpt = 1;   % 0 = fminunc otherwise cmaes
stepTolCmaes = 5e-5; %$$ 2e-4;          %5e-4;  %2e-3;  %0.02;
stepTolFminunc = 1e-6;
funTol = 1e-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tempPlot = false;

OptCore  = @(vMax,omega,duration,strXS,strXC,strYS,strYC,oStr,pShift) ...
  CircleProp(problemOption,nPts,fdOrder,vMax,omega,duration,...
  strXS,strXC,strYS,strYC,oStr,pShift,pShift,pShift,...
  iStart,iOpt, ...
  tempPlot,longerFac,timeRes,propMode,modPlot);

%%%%%%%%%%

addpath('./CMAES/');

OptFun00 = @(vMax,omega,duration,strXS,strXC,strYS,strYC,oStr) ...
  OptimizationFun(DOAVERAGE,nKp,OptCore,vMax,omega,duration,strXS,strXC,strYS,strYC,oStr);

vars  = [vMax;omega;duration;strXS;strXC;strYS;strYC;oStr];

nVars = 8;

%$$ BB    = 0.002;
%$$ BB    = BB * sqrt(mean([strXS,strXC,strYS,strYC].^2));

%$$ stds  = [0.002; 0.002; 0.005; BB*ones(4,1); 0.01]    ;

%AA =0.002;     BB = 0.05;
%AA =0.00002;   BB = 0.0002;
%AA =0.00002;   BB = 0.004;
AA =0.0005;   BB = 0.01;  CC = 0.001;
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

%$$ varFun = @(x) [Inf;10 + 5*sin(x(2)); exp(x(3));x(4:8)];
%$$ varInv = @(x) [Inf;asin((x(2)-10)/5);log(x(3));x(4:8)];

varFun = @(x) [Inf;6 + 3*sin(x(2)); exp(x(3));x(4:8)];
varInv = @(x) [Inf;asin((x(2)-6)/3);log(x(3));x(4:8)];

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


function optVal = OptimizationFun(DOAVERAGE,nKp,OptCore,vMax,omega,duration,strXS,strXC,strYS,strYC,oStr)
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

optVal = MyOptVal(wmax,pmin);

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


function optVal = MyOptVal(wmax,pmin)

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

end



function CircleDomcke(DOAVERAGE, nPhase, nKp, ...
  problemOption,nPts,fdOrder,vMax,omega,duration,...
  strXS,strXC,strYS,strYC,oStr, ...
  iStart,iSel,       longerFac,timeRes,propMode,modPlot, ...
  DOFTLAST, DOFTK, DOCOMPLEXD)

nWay = 2;
nGauge = 2;
nPol = 2;

[dTime,nTime] = TimeParms(omega,timeRes,duration,fdOrder,longerFac);

nTimes = zeros(1,nWay);
nTimes(1) = nTime - fdOrder;
nTimes(2) = nTime ;

  function [eFieldCell,ddtDipCell] = DoProp(phasePos,phaseNeg,phaseMed)
    doPlot = false;
    eFieldCell = cell(nWay,1);
    ddtDipCell = cell(nWay,1);
    for iiway = 1:nWay
      eFieldCell{iiway} = zeros(nTimes(iiway),nPol);
      ddtDipCell{iiway} = zeros(nTimes(iiway),nGauge,nPol);
    end
    [~,~,~,~, ...
      eFieldCell{1}(:,1),   eFieldCell{1}(:,2), ...
      ddtDipCell{1}(:,1,1), ddtDipCell{1}(:,1,2), ddtDipCell{1}(:,2,1), ddtDipCell{1}(:,2,2),...
      eFieldCell{2}(:,1),   eFieldCell{2}(:,2), ...
      ddtDipCell{2}(:,1,1), ddtDipCell{2}(:,1,2), ddtDipCell{2}(:,2,1), ddtDipCell{2}(:,2,2),...
      ] = CircleProp(...
      problemOption,nPts,fdOrder,vMax,omega,duration,...
      strXS,strXC,strYS,strYC,oStr,phasePos,phaseNeg,phaseMed,...
      iStart,iSel,     doPlot,longerFac,timeRes,propMode,modPlot);
  end

phases = (0:nPhase-1)*2*pi/nPhase;
kPhases = (0:nKp-1)*2*pi/nKp;

% REGULAR DOMCKE

eFieldCell = cell(nWay,1);
ddtDip0Cell = cell(nWay,1);
for iway = 1:nWay
  eFieldCell{iway} = zeros(nTimes(iway),nPol,nKp);
  ddtDip0Cell{iway} = zeros(nTimes(iway),nGauge,nPol,nKp);
end

for kphase = 1:nKp
  phase_k = kPhases(kphase);
  [efc,ddc] = DoProp(phase_k,phase_k,phase_k);
  for iway = 1:nWay
    eFieldCell{iway}(:,:,kphase) = efc{iway};
    ddtDip0Cell{iway}(:,:,:,kphase) = ddc{iway};
  end
end

iWay = 1;
plotNum = (iWay-1)*1000;

DomckeCore(DOAVERAGE,nTimes(iWay),dTime,omega,nGauge,nPol,nKp,eFieldCell{iWay},ddtDip0Cell{iWay}, plotNum ) ;

SUBSTUFF = 0;

if DOFTLAST
  nPhase2 = nPhase;
else
  nPhase2 = 1;
end

disp('pause at end of real val domcke.')
if DOCOMPLEXD ~= 0
  disp('about to do complex domcke.')
else
  disp('about to quit.')
end
pause


if DOCOMPLEXD ~= 0
  
  % COMPLEX DOMCKE
  
  ddtDipCell = cell(nWay,1);
  for iway = 1:nWay
    ddtDipCell{iway} = zeros(nTimes(iway),nGauge,nPol,nKp,nPhase,nPhase2);
  end
  
  for kphase = 1:nKp
    for jphase = 1:nPhase2
      for iphase = 1:nPhase
        fprintf('************ DOING i,j,kphase %i %i %i of %i %i %i ************ \n',...
          iphase,jphase,kphase,nPhase,nPhase2,nKp);
        
        phase_i = phases(iphase) + kPhases(kphase);
        phase_j = phases(jphase) + kPhases(kphase);
        
        % phase_i = phases(iphase) ;
        % phase_j = phases(jphase) ;
        
        phase_k = kPhases(kphase);
        [~,ddc] = DoProp(phase_i,phase_j,phase_k);
        
        for iway = 1:nWay
          if SUBSTUFF~=0
            ddtDipCell{iway}(:,:,:,kphase,iphase,jphase) = ddc{iway} - ddtDip0Cell{iway}(:,:,:,kphase);
          else
            ddtDipCell{iway}(:,:,:,kphase,iphase,jphase) = ddc{iway};
          end
        end
      end
    end
  end
  
  ComplexDomckeCore(DOAVERAGE, DOFTK, DOFTLAST, ...
    nTimes(iWay),dTime,omega,nGauge,nPol,nPhase,nKp,...
    eFieldCell{iWay}, ddtDipCell{iWay}, plotNum );
  
end

end



function ftDdtDip0ft = DomckeCore(DOAVERAGE, nTime,dTime,omega,nGauge,nPol,nKp, eField, ddtDip0In,plotNum)

eFieldX = reshape(eField(:,1,:),nTime,nKp);
eFieldY = reshape(eField(:,2,:),nTime,nKp);

ddtDip0 = ddtDip0In;

ftDdtDip0 = FtDip(dTime,ddtDip0);

if DOAVERAGE
  kList = 1:nKp;
else
  kList = 1;
end
nKl = numel(kList);

for ispec=1:2
  ftLenX0{ispec} = reshape(ftDdtDip0{ispec}(:,1,1,kList),[],nKl);
  ftLenY0{ispec} = reshape(ftDdtDip0{ispec}(:,1,2,kList),[],nKl);
  ftVelX0{ispec} = reshape(ftDdtDip0{ispec}(:,2,1,kList),[],nKl);
  ftVelY0{ispec} = reshape(ftDdtDip0{ispec}(:,2,2,kList),[],nKl);
end

ActualResult = @() ...
  DoFtStuffCore(eFieldX(:,kList),eFieldY(:,kList),...
  ftVelX0,ftVelY0,ftLenX0,ftLenY0,dTime,omega,true,true,plotNum);
save('ActualResult','ActualResult')
ActualResult();

if DOAVERAGE
  
  ddtDip0_2 = zeros(nTime,nGauge,nPol,nKp,nKp);
  for ik = 1:nKp
    ddtDip0_2(:,:,:,:,ik) = ddtDip0(:,:,:,mod(ik+(1:nKp)-2,nKp)+1);
  end
  
  ftDdtDip0_2 = FtDip(dTime,ddtDip0_2);
  ftDdtDip0ft = domckeFt(ftDdtDip0_2,5);
  
  for ispec=1:2
    ftLenX0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,1,1,:,:),[],nKp,nKp);
    ftLenY0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,1,2,:,:),[],nKp,nKp);
    ftVelX0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,2,1,:,:),[],nKp,nKp);
    ftVelY0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,2,2,:,:),[],nKp,nKp);
  end
  
  GET = @(a,ik) {a{1}(:,:,ik),a{2}(:,:,ik)};
  
else
  
  ftDdtDip0ft = domckeFt(ftDdtDip0,4);
  
  for ispec=1:2
    ftLenX0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,1,1,:),[],nKp);
    ftLenY0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,1,2,:),[],nKp);
    ftVelX0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,2,1,:),[],nKp);
    ftVelY0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,2,2,:),[],nKp);
  end
  
  GET = @(a,ik) {a{1}(:,ik),a{2}(:,ik)};
  
end

DomckeResult = @(ik) ...
  DoFtStuffCore(eFieldX(:,kList),eFieldY(:,kList),...
  GET(ftVelX0ft,ik),GET(ftVelY0ft,ik),GET(ftLenX0ft,ik),GET(ftLenY0ft,ik),...
  dTime,omega,true,true,plotNum);

save('DomckeResult','DomckeResult')
DomckeResult(2);

end


function x = domckeFt(x,ind)
ispec = 1;   % pos frequency
xsize = size(x{ispec});
x{ispec} = fft(x{ispec},[],ind)/xsize(ind);
ispec = 2;   % neg freqency
xsize = size(x{ispec});
x{ispec} = fft(FtFlip(x{ispec},ind),[],ind)/xsize(ind);
end


function x = FtFlip(x,ind)
if ind == 2
  x(:,2:end,:)           = flip(x(:,2:end,:),ind);
elseif ind == 4
  x(:,:,:,2:end,:)       = flip(x(:,:,:,2:end,:),ind);
elseif ind == 5
  x(:,:,:,:,2:end,:)     = flip(x(:,:,:,:,2:end,:),ind);
elseif ind == 6
  x(:,:,:,:,:,2:end,:)   = flip(x(:,:,:,:,:,2:end,:),ind);
elseif ind == 7
  x(:,:,:,:,:,:,2:end,:) = flip(x(:,:,:,:,:,:,2:end,:),ind);
else
  error('not supported')
end
end


function ComplexDomckeCore(DOAVERAGE, DOFTK, DOFTLAST, nTime,dTime,omega,nGauge,nPol,nPhase,nKp,...
  eField,  ddtDipIn, plotNum )

if DOFTK
  nFt = nKp;
else
  nFt = 1;
end
if DOFTLAST
  nPhase2 = nPhase;
else
  nPhase2 = 1;
end
if DOAVERAGE
  nKl = nKp;
else
  nKl = 1;
end

kList = 1:nKl;
eFieldX = reshape(eField(:,1,kList),nTime,nKl);
eFieldY = reshape(eField(:,2,kList),nTime,nKl);
clear kList;

if DOFTK
  if DOAVERAGE
    ddtDip = zeros(nTime,nGauge,nPol,nKp,nKp,nPhase,nPhase2);
    for ik = 1:nKp
      pk = mod(ik+(1:nKp)-2,nKp)+1 ;
      ddtDip(:,:,:,:,ik,:,:) = ddtDipIn(:,:,:,pk,:,:);
    end
  else
    ddtDip = reshape(ddtDipIn,nTime,nGauge,nPol,1,nKp,nPhase,nPhase2);
  end
else
  if DOAVERAGE
    ddtDip = reshape(ddtDipIn,nTime,nGauge,nPol,nKp,1,nPhase,nPhase2);
  else
    if nKp ~= 1
      error('what?')
    end
    ddtDip = reshape(ddtDipIn,nTime,nGauge,nPol,1,1,nPhase,nPhase2);
  end
end

ftDdtDip = FtDip(dTime,ddtDip);

if DOFTLAST
  if DOFTK
    ftDdtDipft = domckeFt(domckeFt(domckeFt(ftDdtDip,5),6),7);
  else
    ftDdtDipft = domckeFt(domckeFt(ftDdtDip,6),7);
  end
else
  if DOFTK
    ftDdtDipft = domckeFt(domckeFt(ftDdtDip,5),6);
  else
    ftDdtDipft = domckeFt(ftDdtDip,6);
  end
end

for ispec=1:2
  ftLenXft{ispec} = reshape(ftDdtDipft{ispec}(:,1,1,:),[],nKl,nFt,nPhase,nPhase2);
  ftLenYft{ispec} = reshape(ftDdtDipft{ispec}(:,1,2,:),[],nKl,nFt,nPhase,nPhase2);
  ftVelXft{ispec} = reshape(ftDdtDipft{ispec}(:,2,1,:),[],nKl,nFt,nPhase,nPhase2);
  ftVelYft{ispec} = reshape(ftDdtDipft{ispec}(:,2,2,:),[],nKl,nFt,nPhase,nPhase2);
end

GET = @(a,kp,ip,jp) {a{1}(:,:,kp,ip,jp),a{2}(:,:,kp,ip,jp)};

CDomckeResult = @(kp,ip,jp) ...
  DoFtStuffCore(eFieldX,eFieldY, ...
  GET(ftVelXft,kp,ip,jp), GET(ftVelYft,kp,ip,jp), GET(ftLenXft,kp,ip,jp), GET(ftLenYft,kp,ip,jp), ...
  dTime,omega,true,true,plotNum);

save('CDomckeResult','CDomckeResult');

DoAllCD = @()  DoAllCDResults(nFt,nPhase,nPhase2,CDomckeResult);

save('DoAllCD','DoAllCD');

DoAllCD();

CDomckeResult(1,2,1);

end


function DoAllCDResults(nFt,nPh,nPh2,CDomckeResult)
allResults = zeros(nFt,nPh,nPh2);
for ift = 1:nFt
  for ip = 1:nPh
    % for jp = ip:nPh2
    for jp = 1:nPh2
      allResults(ift,ip,jp) = CDomckeResult(ift,ip,jp);
    end
  end
end
[~,sord] = sort(-abs(allResults));
maxResult = abs(allResults(sord(1)));
iBig = abs(allResults) > 1e-5 * maxResult;
allResults = allResults(iBig);
[IFT,IPH,JPH] = ndgrid(1:nFt,1:nPh,1:nPh2);
IFT = IFT(iBig);
IPH = IPH(iBig);
JPH = JPH(iBig);
nBig = sum(iBig(:));
for ii=1:nBig
  fprintf(' %i  %i  %i  %e  \n ',IFT(ii),IPH(ii),JPH(ii),allResults(ii));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   CircleProp    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [popVel,popLen,workVel,workLen,...
  eFieldXdip, eFieldYdip, dipLenX, dipLenY, dipVelX, dipVelY,...
  eFieldXder, eFieldYder, derLenX, derLenY, derVelX, derVelY ] = ...
  CircleProp( ...
  problemOption,nPts,fdOrder,vMax,omega,duration,...
  strXS,strXC,strYS,strYC,oStr,phasePos,phaseNeg,phaseMed,...
  iStart,iSel,     doPlot,longerFac,timeRes,propMode,modPlot)

DOPROJ = 1==0 ;
nQuad1       = 5;

DOPARITY = 1==0 ;

DOCONJ = false;
if phasePos ~= phaseNeg
  DOCONJ = true;
end

startVec = CircleEigs(problemOption,nPts,fdOrder,vMax,[iStart,iSel],doPlot,false,omega);
eigVects = startVec(:,2:end);
startVec = startVec(:,1);

startVecp = [];
eigVectsp = [];
if DOPARITY
  startVecp = CircleEigs(problemOption,nPts,fdOrder,vMax,[iStart,iSel],false,true,omega);
  eigVectsp = startVecp(:,2:end);
  startVecp = startVecp(:,1);
end

clear iStart iSel

nEig = size(eigVects,2);

stateProj = [];
stateProjp = [];
if DOPROJ
  stateProj = eigVects * eigVects';
  stateProjp = eigVectsp * eigVectsp';
end

[fdmat,~,ham0,~,dipoleOpX,dipoleOpY,thetaVals] = CircleParts(problemOption,nPts,fdOrder,vMax,false);

ham0p = [];
if DOPARITY
  [~,~,ham0p] = CircleParts(problemOption,nPts,fdOrder,vMax,true);
end

if propMode==0   % do sparse for Crank-Nicholson...
  dipoleOpX = sparse( dipoleOpX );
  dipoleOpY = sparse( dipoleOpY );
  fdmat        = sparse( fdmat );
  ham0         = sparse( ham0 );
  ham0p        = sparse( ham0p );
  EYE          = speye( nPts );
else
  dipoleOpX = full( dipoleOpX );
  dipoleOpY = full( dipoleOpY );
  fdmat        = full( fdmat );
  ham0         = full( ham0 );
  ham0p        = full( ham0p );
  EYE          = eye( nPts );
end

  function  [   sovl, popLen, el ...
      ,         sovv, popVel, ev] = OvlStuffCore(...
      vl,vv,vlc,vvc,vlp,vvp,vlcp,vvcp,    startVec, startVecp, ...
      eigVects, eigVectsp, ham0, ham0p,  DOPARITY)
    
    sovl   = ((startVec'*vl).*conj(startVec'*vlc));
    sovv   = ((startVec'*vv).*conj(startVec'*vvc));
    popLen = ((eigVects'*vl).*conj(eigVects'*vlc));
    popVel = ((eigVects'*vv).*conj(eigVects'*vvc));
    ev     = vvc' * ham0 * vv ;
    el     = vlc' * ham0 * vl ;
    
    if DOPARITY
      sovlp   = ((startVecp'*vlp).*conj(startVecp'*vlcp));
      sovvp   = ((startVecp'*vvp).*conj(startVecp'*vvcp));
      popLenp = ((eigVectsp'*vlp).*conj(eigVectsp'*vlcp));
      popVelp = ((eigVectsp'*vvp).*conj(eigVectsp'*vvcp));
      evp     = vvcp' * ham0p * vvp ;
      elp     = vlcp' * ham0p * vlp ;
      sovl    =    (   sovl   +       sovlp   ) / 2;
      sovv    =    (   sovv   +       sovvp   ) / 2;
      popLen  =    (   popLen +       popLenp ) / 2;
      popVel  =    (   popVel +       popVelp ) / 2;
      ev      =    (   ev     +       evp     ) / 2;
      el      =    (   el     +       elp     ) / 2;
    end
    
  end

OvlStuff = @(vl,vv,vlc,vvc,vlp,vvp,vlcp,vvcp) OvlStuffCore(...
  vl,vv,vlc,vvc, ...
  vlp,vvp,vlcp,vvcp, ...
  startVec, startVecp, eigVects, eigVectsp, ham0, ham0p, DOPARITY);

clear OvlStuffCore eigVects eigVectsp

%%%%%%%%%%%%%  TIME   %%%%%%%%%%%

[dTime,nTime] = TimeParms(omega,timeRes,duration,fdOrder,longerFac);

disp(['duration ' num2str(duration)])
disp(['dTime    ' num2str(dTime)])
disp(['nTime    ' num2str(nTime)])

startTimes = (0:(nTime-2))' * dTime;
% halfTimes = (0.5:(nTime-1.5)) * dTime;
endTimes = (1:(nTime-1))' * dTime;
dataTimes = (0:(nTime-1))' * dTime;

if 1==0
  GetFieldsHere = @(intimes,doplot,shiftPos,shiftNeg,shiftMed) ...
    GetFields3(intimes,strXS,strXC,strYS,strYC,oStr,omega,duration,shiftPos,shiftNeg,shiftMed, doplot);
else
  GetFieldsHere = @(intimes,doplot,shiftPos,shiftNeg,~) ...
    GetFields3B(intimes,strXS,strXC,strYS,strYC,oStr,omega,duration,shiftPos,shiftNeg, doplot);
end

% [~,~,eFieldXder,eFieldYder] = GetFieldsHere(dataTimes,doPlot,phasePos,phaseNeg,phaseMed);

[aFieldXder,aFieldYder,eFieldXder,eFieldYder] = ...
  GetFieldsHere(dataTimes,doPlot,phasePos,phaseNeg,phaseMed)...
  ; %#ok<ASGLU> ok unused; for debug

%%%%%%%%%%%%%%%%%    Dipole Operators     %%%%%%%%%%%%%%%%%%%%%%%%%%

dipoleOpLenX = -dipoleOpX ;
dipoleOpVelX =  dipoleOpY;
dipoleOpLenY =  dipoleOpY;
dipoleOpVelY =  dipoleOpX;

%$$ elseif EXPECTMODE == 1
%$$   expectOpX = 1i/2 * (fdmat * dipoleOpY + dipoleOpY * fdmat);
%$$   expectOpY = 1i/2 * (fdmat * dipoleOpX + dipoleOpX * fdmat);

expectOpX = dipoleOpX;
expectOpY = -dipoleOpY;

if DOPARITY
  GetExpectX = @(vc,v,vcp,vp) 1/2 * ( ...
    vc' * expectOpX * v + ...
    vcp' * expectOpX * vp );
  GetExpectY = @(vc,v,vcp,vp) 1/2 * (  ...
    vc' * expectOpY * v + ...
    vcp' * expectOpY * vp );
  GetDerivX = @(vc,v,ham,vcp,vp,hamp) 1/2 * ((-1i)*( ...
    (vc'  * ham  * dipoleOpX * v  - vc'  * dipoleOpX * ham  * v ) + ...
    (vcp' * hamp * dipoleOpX * vp - vcp' * dipoleOpX * hamp * vp ) ));
  GetDerivY = @(vc,v,ham,vcp,vp,hamp) 1/2 * ((+1i)*( ...
    (vc'  * ham  * dipoleOpY * v  - vc'  * dipoleOpY * ham  * v ) + ...
    (vcp' * hamp * dipoleOpY * vp - vcp' * dipoleOpY * hamp * vp ) ));
else
  GetExpectX = @(vc,v,~,~) vc' * expectOpX * v;
  GetExpectY = @(vc,v,~,~) vc' * expectOpY * v;
  GetDerivX = @(vc,v,ham,~,~,~) (-1i) * (vc' * ham * dipoleOpX * v - vc' * dipoleOpX * ham * v );
  GetDerivY = @(vc,v,ham,~,~,~) (+1i) * (vc' * ham * dipoleOpY * v - vc' * dipoleOpY * ham * v );
end

clear dipoleOpX dipoleOpY expectOpX expectOpY

%%%%%%%%%%%%%%%%%%  Hamiltonians   %%%%%%%%%%%%%%%%%

% HamVel = @(aval) (1/2) * (1i*fdmat - aval*dipoleOpVel)' * (1i*fdmat - aval*dipoleOpVel) + diag(potential);

HamLen = @(evalx,evaly) ham0 ...
  + evalx * dipoleOpLenX ...
  + evaly * dipoleOpLenY ;
HamLenp = @(evalx,evaly) ham0p ...
  + evalx * dipoleOpLenX ...
  + evaly * dipoleOpLenY ;

% HamVel = @(avalx,avaly) ham0 ...
%   - 1i/2 * avalx * (fdmat*dipoleOpVelX + dipoleOpVelX*fdmat) ...
%   - 1i/2 * avaly * (fdmat*dipoleOpVelY + dipoleOpVelY*fdmat) ...
%   + 1/2 * (avalx * dipoleOpVelX + avaly * dipoleOpVelY)^2;

hamVelFacX = (fdmat*dipoleOpVelX + dipoleOpVelX*fdmat) ;
hamVelFacY = (fdmat*dipoleOpVelY + dipoleOpVelY*fdmat) ;
dOVXX = dipoleOpVelX^2;
dOVYY = dipoleOpVelY^2;
dOVcross = dipoleOpVelX * dipoleOpVelY + dipoleOpVelY * dipoleOpVelX;

%$$ HamVel = @(avalx,avaly) ham0 ...
%$$   - 1i/2 * avalx * hamVelFacX ...
%$$   - 1i/2 * avaly * hamVelFacY ...
%$$   + 1/2 * (avalx * dipoleOpVelX + avaly * dipoleOpVelY)^2;

HamVel = @(avalx,avaly) ham0 ...
  - 1i/2 * avalx * hamVelFacX ...
  - 1i/2 * avaly * hamVelFacY ;
HamVel = @(avalx,avaly) HamVel(avalx,avaly) ...
  + 1/2 * (avalx^2 * dOVXX + avaly^2 * dOVYY + avalx*avaly * dOVcross);
HamVelp = @(avalx,avaly) ham0p ...
  - 1i/2 * avalx * hamVelFacX ...
  - 1i/2 * avaly * hamVelFacY ;
HamVelp = @(avalx,avaly) HamVelp(avalx,avaly) ...
  + 1/2 * (avalx^2 * dOVXX + avaly^2 * dOVYY + avalx*avaly * dOVcross);

clear dipoleOpLenX dipoleOpLenY dipoleOpVelX dipoleOpVelY ...
  hamVelFacX hamVelFacY dOVXX dOVYY dOVcross

  function hamvel = VelHam0(intime,GetFieldsHere,phasePos,phaseNeg,phaseMed,HamVel)
    [afx,afy,~,~] =  GetFieldsHere(intime,false,phasePos,phaseNeg,phaseMed);
    hamvel = HamVel(afx,afy);
  end
  function hamvel = LenHam0(intime,GetFieldsHere,phasePos,phaseNeg,phaseMed,HamLen)
    [~,~,efx,efy] =  GetFieldsHere(intime,false,phasePos,phaseNeg,phaseMed);
    hamvel = HamLen(efx,efy);
  end

VelHam = @(t) VelHam0(t,GetFieldsHere,phasePos,phaseNeg,phaseMed,HamVel);
LenHam = @(t) LenHam0(t,GetFieldsHere,phasePos,phaseNeg,phaseMed,HamLen);
VelHamp = @(t) VelHam0(t,GetFieldsHere,phasePos,phaseNeg,phaseMed,HamVelp);
LenHamp = @(t) LenHam0(t,GetFieldsHere,phasePos,phaseNeg,phaseMed,HamLenp);
clear VelHam0 LenHam0 HamVel HamLen HamVelp HamLenp

  function [hl,hv] = TDHams0(intime,dTime,LenHam,VelHam,nQuad1,DOPROJ,stateProj)
    hv = Magnus5(intime,dTime,VelHam,nQuad1);
    hl = Magnus5(intime,dTime,LenHam,nQuad1);
    if DOPROJ
      hv = stateProj * hv * stateProj;
      hl = stateProj * hl * stateProj;
    end
  end

TDHams  = @(intime) TDHams0(intime,dTime,LenHam, VelHam, nQuad1,DOPROJ,stateProj);
TDHamsp = @(intime) TDHams0(intime,dTime,LenHamp,VelHamp,nQuad1,DOPROJ,stateProjp);
clear TDHams0 stateProj stateProjp nQuad1

  function  [vl,vv,vlc,vvc] = PropAllCore(hl,hv,  vl,vv,vlc,vvc,  propMode, DOCONJ, EYE, dTime, PropIt0, usePropIt0)
    if usePropIt0
      LPropIt  = PropIt0;
      VPropIt  = PropIt0;
      LPropItc = PropIt0;
      VPropItc = PropIt0;
    else
      if propMode == 0
        ml = (EYE + 1i*dTime/2 * hl);
        LPropIt = @(v) ml \ ( (v - 1i*dTime/2 * hl * v ) );
        mv = (EYE + 1i*dTime/2 * hv);
        VPropIt = @(v) mv \ ( (v - 1i*dTime/2 * hv * v ) );
      else
        lU = expm(-1i*dTime * hl) ;
        vU = expm(-1i*dTime * hv) ;
        LPropIt = @(v) lU * v;
        VPropIt = @(v) vU * v;
      end
      if DOCONJ
        if propMode == 0
          LPropItc = @(v) (EYE + 1i*dTime/2 * hl') \ ( (v - 1i*dTime/2 * hl' * v ) );
          VPropItc = @(v) (EYE + 1i*dTime/2 * hv') \ ( (v - 1i*dTime/2 * hv' * v ) );
        else
          lU = expm(-1i*dTime * hl') ;
          vU = expm(-1i*dTime * hv') ;
          LPropItc = @(v) lU * v;
          VPropItc = @(v) vU * v;
        end
      end
    end
    vv = VPropIt(vv);
    vl = LPropIt(vl);
    if DOCONJ
      vvc = VPropItc(vvc);
      vlc = LPropItc(vlc);
    else
      vvc = vv;
      vlc = vl;
    end
  end

if propMode == 0
  U0 = (EYE + 1i*dTime/2 * ham0) \ (EYE - 1i*dTime/2 * ham0) ;
else
  U0 = expm(-1i*dTime * ham0);
end
U0p = [];
if DOPARITY
  if propMode == 0
    U0p = (EYE + 1i*dTime/2 * ham0p) \ (EYE - 1i*dTime/2 * ham0p) ;
  else
    U0p = expm(-1i*dTime * ham0p);
  end
end

PropIt0 = @(v) U0 * v ;
PropAll = @(hl,hv,  vl,vv,vlc,vvc)  PropAllCore(hl,hv,  vl,vv,vlc,vvc,  propMode, DOCONJ, EYE, dTime, PropIt0, false);
PropAll0 = @(hl,hv,  vl,vv,vlc,vvc)  PropAllCore(hl,hv,  vl,vv,vlc,vvc,  propMode, DOCONJ, EYE, dTime, PropIt0, true);
PropIt0p = @(v) U0p * v ;
PropAllp = @(hl,hv,  vl,vv,vlc,vvc)  PropAllCore(hl,hv,  vl,vv,vlc,vvc,  propMode, DOCONJ, EYE, dTime, PropIt0p, false);
PropAll0p = @(hl,hv,  vl,vv,vlc,vvc)  PropAllCore(hl,hv,  vl,vv,vlc,vvc,  propMode, DOCONJ, EYE, dTime, PropIt0p, true);

clear PropIt0 PropIt0p U0 U0p PropAllCore propMode

  function  [nsql,nsqv] = VecNorm(vl,vv,vlc,vvc,DOPARITY,vlp,vvp,vlcp,vvcp)
    nsqv = vvc'*vv ;
    nsql = vlc'*vl ;
    if DOPARITY
      nsqvp = vvcp'*vvp ;
      nsqlp = vlcp'*vlp ;
      % nsqv = sqrt(nsqv * nsqvp);
      % nsql = sqrt(nsql * nsqlp);
      nsqv = 1/2 * (nsqv + nsqvp);
      nsql = 1/2 * (nsql + nsqlp);
    end
  end

  function [ expLenX,expLenY,expVelX,expVelY, ...
      derLenX,derLenY,derVelX,derVelY ] = ...
      AllMatelsCore(intime,vl,vv,vlc,vvc, vlp,vvp,vlcp,vvcp, ham0, ham0p, ...
      LenHam, VelHam, LenHamp, VelHamp, GetExpectX, GetExpectY, GetDerivX, GetDerivY, DOPARITY)
    
    if intime < 0
      hve = ham0;
      hle = ham0;
      hlep = ham0p;
      hvep = ham0p;
    else
      hve = VelHam(intime);
      hle = LenHam(intime);
      hvep = [];
      hlep = [];
      if DOPARITY
        hvep = VelHamp(intime);
        hlep = LenHamp(intime);
      end
    end
    expVelX = GetExpectX(vvc , vv, vvcp , vvp ) ;
    expLenX = GetExpectX(vlc , vl, vlcp , vlp ) ;
    expVelY = GetExpectY(vvc , vv, vvcp , vvp ) ;
    expLenY = GetExpectY(vlc , vl, vlcp , vlp) ;
    
    derVelX = GetDerivX(vvc , vv, hve, vvcp , vvp, hvep ) ;
    derLenX = GetDerivX(vlc , vl, hle, vlcp , vlp, hlep) ;
    derVelY = GetDerivY(vvc , vv, hve, vvcp , vvp, hvep ) ;
    derLenY = GetDerivY(vlc , vl, hle, vlcp , vlp, hlep) ;
  end

AllMatels = @(intime,vl,vv,vlc,vvc,vlp,vvp,vlcp,vvcp) ...
  AllMatelsCore(intime,vl,vv,vlc,vvc,vlp,vvp,vlcp,vvcp, ham0, ham0p, ...
  LenHam, VelHam, LenHamp, VelHamp, GetExpectX, GetExpectY, GetDerivX, GetDerivY, DOPARITY);

clear AllMatelsCore LenHam VelHam LenHamp VelHamp
clear  GetExpectX GetExpectY GetDerivX GetDerivY ham0 ham0p;

%   function [vl,vv,vlc,vvc,vlp,vvp,vlcp,vvcp] = DealV(vCell)
%     % nGauge,nConj,nParity
%     vl = vCell{1,1,1};
%     vv = vCell{2,1,1};
%     vlc = vCell{1,2,1};
%     vvc = vCell{2,2,1};
%     if DOPARITY
%       vlp = vCell{1,1,2};
%       vvp = vCell{2,1,2};
%       vlcp = vCell{1,2,2};
%       vvcp = vCell{2,2,2};
%     end
%   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[   ~,~, eStart ]                       = OvlStuff(...
  startVec,startVec,startVec,startVec,  startVecp,startVecp,startVecp,startVecp);
[ exStart,eyStart,~,~, dxStart,dyStart] =   AllMatels(-99, ...
  startVec,startVec,startVec,startVec,  startVecp,startVecp,startVecp,startVecp);

%   eStart = startVec' * ham0 * startVec ;
%
%   exStart = GetExpectX(startVec,startVec);
%   eyStart = GetExpectY(startVec,startVec);
%
%   dxStart = GetDerivX(startVec,startVec,ham0);
%   dyStart = GetDerivY(startVec,startVec,ham0);
%   clear GetExpectX GetExpectY GetDerivX GetDerivY ham0 ham0p;
%
%   if DOPARITY
%     error('notdone')
%   end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

enVel = zeros(nTime,1);
enLen = zeros(nTime,1);
eOvlVel = zeros(nEig,nTime);
eOvlLen = zeros(nEig,nTime);

enVel(1) = eStart;
enLen(1) = eStart;

expVelX = zeros(nTime,1);
expLenX = zeros(nTime,1);
expVelY = zeros(nTime,1);
expLenY = zeros(nTime,1);
expVelX(1) = exStart;
expLenX(1) = exStart;
expVelY(1) = eyStart;
expLenY(1) = eyStart;

derVelX = zeros(nTime,1);
derLenX = zeros(nTime,1);
derVelY = zeros(nTime,1);
derLenY = zeros(nTime,1);
derVelX(1) = dxStart;
derLenX(1) = dxStart;
derVelY(1) = dyStart;
derLenY(1) = dyStart;

vv   = startVec;
vl   = vv;
vvc  = vv;
vlc  = vv;
vvp  = startVecp;
vlp  = vvp;
vvcp = vvp;
vlcp = vvp;

clear startVec startVecp dxStart dyStart exStart eyStart

  function NicePrint(dispArr)
    mmyisreal = @(x) norm(imag(x)) < 1e-10 * norm(real(x)) ;
    fprintf('%7.2f ',real(dispArr(1)))
    fprintf(' %9.5f ',real(dispArr(2:end)))
    fprintf('\n');
    if ~mmyisreal(dispArr)
      fprintf('        ')
      fprintf(' %9.5f ',imag(dispArr(2:end)))
      fprintf('\n');
    end
  end

NicePrint([0,0,0,eStart,eStart,1,1])

workLenX = 0;
workLenY = 0;
workVelX = 0;
workVelY = 0;

pulseDone = false;
pulseDone2 = false;
pulseDone3 = false;

for itime = 2:nTime
  
  starttime = startTimes(itime-1);
  endtime = endTimes(itime-1);
  
  if starttime > duration  && 1==1
    pulseDone = true;
  end
  if pulseDone
    PropAll = PropAll0;
    PropAllp = PropAll0p;
  else
    [hl,hv] = TDHams(starttime);
    if DOPARITY
      [hlp,hvp] = TDHamsp(starttime);
    end
  end
  
  [vl,vv,vlc,vvc] = PropAll(hl,hv,  vl,vv,vlc,vvc);
  if DOPARITY
    [vlp,vvp,vlcp,vvcp] = PropAllp(hlp,hvp,  vlp,vvp,vlcp,vvcp);
  end
  
  [nsql,nsqv] = VecNorm(vl,vv,vlc,vvc,DOPARITY,vlp,vvp,vlcp,vvcp);
  
  [   sovl, popLen, el ...
    , sovv, popVel, ev] = OvlStuff(vl,vv,vlc,vvc,vlp,vvp,vlcp,vvcp);
  
  eOvlLen(:,itime) = popLen;
  eOvlVel(:,itime) = popVel;
  
  enVel(itime) = ev;
  enLen(itime) = el;
  
  [ expLenX(itime),expLenY(itime),expVelX(itime),expVelY(itime), ...
    derLenX(itime),derLenY(itime),derVelX(itime),derVelY(itime) ] = ...
    AllMatels(endtime,vl,vv,vlc,vvc,vlp,vvp,vlcp,vvcp);
  
  if ~pulseDone
    [~,~,efxtemp,efytemp] =  GetFieldsHere(endtime,false,phasePos,phaseNeg,phaseMed);
    ff = -1;
    if itime == nTime
      ff = -0.5;
    end
    % efxtemp = real(efxtemp);
    % efytemp = real(efytemp);
    workVelX = workVelX + derVelX(itime) * efxtemp * dTime * ff;
    workVelY = workVelY + derVelY(itime) * efytemp * dTime * ff;
    workLenX = workLenX + derLenX(itime) * efxtemp * dTime * ff;
    workLenY = workLenY + derLenY(itime) * efytemp * dTime * ff;
  end
  
  if itime == nTime || ( mod(itime,modPlot) == 1 && ~pulseDone2 )
    dispArr = [endtime,nsqv-1,nsql-1,ev,el,sovv,sovl];
    NicePrint(dispArr)
    pulseDone2 = pulseDone;
  end
  FS = 22;
  
  if doPlot  %%%% || itime == nTime
    if mod(itime,modPlot) == 1 || itime == nTime
      figure(103)
      fac = 10;
      %  thetaVals,fac*abs(vv.^2) ,'o',...
      %  thetaVals,fac*abs(vl.^2) ,'.','LineWidth',2)
      plot(...
        thetaVals,fac*real(vv.*conj(vvc)) ,'o',...
        thetaVals,fac*real(vl.*conj(vlc)) ,'.','LineWidth',2)
      title(['Density t= ' num2str(itime*dTime)]);
    end
    if ~pulseDone3
      if ( mod(itime,modPlot*100) == 1 || itime == nTime || pulseDone )
        figure(203)
        plot((1:itime)*dTime,real(enLen(1:itime)),(1:itime)*dTime,real(enVel(1:itime)),'LineWidth',2)
        legend('Len','Vel');
        xlabel('Time (atomic units)')
        title(['<Psi(t)|H_0|Psi(t)> t= ' num2str(itime*dTime)]);

        figure(204)
        semilogy((1:itime)*dTime,real(eOvlLen(:,1:itime)),'LineWidth',2)
        title(['State Overlaps Length t= ' num2str(itime*dTime)]);
        set(gca,'FontSize',FS)
        ylim([1e-8 inf])
        xlabel('Time (atomic units)')

        figure(205)
        semilogy((1:itime)*dTime,real(eOvlVel(:,1:itime)),'LineWidth',2)
        title(['State Overlaps Velocity t= ' num2str(itime*dTime)]);
        set(gca,'FontSize',FS)
        ylim([1e-8 inf])
        xlabel('Time (atomic units)')
        
        pulseDone3 = pulseDone;
      end
    end
    drawnow
  end
end

disp('POPS: log(v) log(l) v l')
disp([(1:nEig)',log([popVel,popLen]),[popVel,popLen]])

workVel = ev - eStart;
workLen = el - eStart;

workVelTot = workVelX + workVelY;
workLenTot = workLenX + workLenY;

disp('WORK from change in energy')

if any(abs(imag([workVel,workLen,workVelTot,workLenTot])) > 1e-8)
  fprintf('WORK:     Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(workVel),real(workLen),imag(workVel),imag(workLen));
  disp('Work integral dt from commutator')
  fprintf('WORK X:   Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(workVelX),real(workLenX),imag(workVelX),imag(workLenX));
  fprintf('WORK Y:   Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(workVelY),real(workLenY),imag(workVelY),imag(workLenY));
  % fprintf('WORK P: Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(workVelP),real(workLenP),imag(workVelP),imag(workLenP));
  % fprintf('WORK M: Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(workVelM),real(workLenM),imag(workVelM),imag(workLenM));
  fprintf('WORK Tot: Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(workVelTot),real(workLenTot),imag(workVelTot),imag(workLenTot));
else
  fprintf('WORK:     Vel, Len   %13.9f  %13.9f    \n',real(workVel),real(workLen));
  disp('Work integral dt from commutator')
  fprintf('WORK X:   Vel, Len   %13.9f  %13.9f    \n',real(workVelX),real(workLenX));
  fprintf('WORK Y:   Vel, Len   %13.9f  %13.9f    \n',real(workVelY),real(workLenY));
  % fprintf('WORK P: Vel, Len   %13.9f  %13.9f    \n',real(workVelP),real(workLenP));
  % fprintf('WORK M: Vel, Len   %13.9f  %13.9f    \n',real(workVelM),real(workLenM));
  fprintf('WORK Tot: Vel, Len   %13.9f  %13.9f    \n',real(workVelTot),real(workLenTot));
end

if nargout < 5 && ~doPlot
  return
end

PadStart = @(x) [ones(fdOrder,1)*x(1);x];
ChopPad = @(x) x(fdOrder+1:end-fdOrder);
fdvec = FDVec(fdOrder);
DoDiff = @(x) ChopPad(FDMult(PadStart(x),fdvec,fdOrder)) / dTime * (-1);
dipVelX = DoDiff(expVelX);  % dip is time derivative of exp
dipLenX = DoDiff(expLenX);  % should be equal to der
dipVelY = DoDiff(expVelY);
dipLenY = DoDiff(expLenY);

eFieldXdip = eFieldXder(1:end-fdOrder);
eFieldYdip = eFieldYder(1:end-fdOrder);

clear om;

eFieldX = eFieldXdip;
eFieldY = eFieldYdip;
velX = dipVelX;
velY = dipVelY;
lenX = dipLenX;
lenY = dipLenY;

% if ~doPlot
%   return;              %%%% RETURN if no plot
% end

plotNum = 0;
[cumVelT,cumVelX,cumVelY,cumVelP,cumVelM, ...
  cumLenT,cumLenX,cumLenY,cumLenP,cumLenM ] = ...
  DoFtStuff(eFieldX,eFieldY,velX,velY,lenX,lenY,dTime,omega,doPlot,doPlot,plotNum);

disp('Work integral domega:  With dipole operator, then differentiate:')
DispCumulative();

eFieldX = eFieldXder;
eFieldY = eFieldYder;
velX = derVelX;
velY = derVelY;
lenX = derLenX;
lenY = derLenY;

plotNum = 1000;
[cumVelT,cumVelX,cumVelY,cumVelP,cumVelM, ...
  cumLenT,cumLenX,cumLenY,cumLenP,cumLenM ] = ...
  DoFtStuff(eFieldX,eFieldY,velX,velY,lenX,lenY,dTime,omega,doPlot,doPlot,plotNum);

disp('Work integral domega:  from commmutator:')
DispCumulative();

if 1==0
  emissVelX = DoDiff(derVelX);  % d/dt of (d/dt D computed from commutator)
  emissLenX = DoDiff(derLenX);
  emissVelY = DoDiff(derVelY);
  emissLenY = DoDiff(derLenY);
  
  eFieldX = eFieldXdip;
  eFieldY = eFieldYdip;
  velX = emissVelX;
  velY = emissVelY;
  lenX = emissLenX;
  lenY = emissLenY;
  
  plotNum = 2000;
  DoFtStuff(eFieldX,eFieldY,velX,velY,lenX,lenY,dTime,omega,0,doPlot,plotNum);
end


disp(' ')
disp('OKAY done with run.');
%  hit enter')
% pause

return;  % END CircleProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function DispCumulative()
    if any(abs(imag([cumVelT,cumLenT])) > 1e-8)
      fprintf('WORK Tot: Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(cumVelT),real(cumLenT),imag(cumVelT),imag(cumLenT));
      fprintf('WORK X:   Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(cumVelX),real(cumLenX),imag(cumVelX),imag(cumLenX));
      fprintf('WORK Y:   Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(cumVelY),real(cumLenY),imag(cumVelY),imag(cumLenY));
      fprintf('WORK +:   Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(cumVelP),real(cumLenP),imag(cumVelP),imag(cumLenP));
      fprintf('WORK -:   Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(cumVelM),real(cumLenM),imag(cumVelM),imag(cumLenM));
    else
      fprintf('WORK Tot: Vel, Len   %13.9f  %13.9f \n',real(cumVelT),real(cumLenT));
      fprintf('WORK X:   Vel, Len   %13.9f  %13.9f \n',real(cumVelX),real(cumLenX));
      fprintf('WORK Y:   Vel, Len   %13.9f  %13.9f \n',real(cumVelY),real(cumLenY));
      fprintf('WORK +:   Vel, Len   %13.9f  %13.9f \n',real(cumVelP),real(cumLenP));
      fprintf('WORK -:   Vel, Len   %13.9f  %13.9f \n',real(cumVelM),real(cumLenM));
    end
  end

end


function varargout = FtDip(dTime,varargin)
if nargin-1 ~= nargout
  error('oops')
end
varargout = cell(nargout,1);
for iarg = 1:nargout
  [varargout{iarg}{1},varargout{iarg}{2}] = TimeFT(varargin{iarg},dTime);
end
end


function [cumVelT,cumVelX,cumVelY,cumVelP,cumVelM, ...
  cumLenT,cumLenX,cumLenY,cumLenP,cumLenM ] = ...
  DoFtStuff(...
  eFieldX,eFieldY,dipVelX,dipVelY,dipLenX,dipLenY,...
  dTime,omega,doPlot,doPlot2,plotNum)

% efield and TIME DERIVATIVE of dipole moments input
[ftVelX,ftVelY,ftLenX,ftLenY] = FtDip(dTime,dipVelX,dipVelY,dipLenX,dipLenY);

[cumVelT,cumVelX,cumVelY,cumVelP,cumVelM, ...
  cumLenT,cumLenX,cumLenY,cumLenP,cumLenM ] = ...
  DoFtStuffCore(...
  eFieldX, eFieldY, ftVelX, ftVelY, ftLenX, ftLenY, ...
  dTime,omega,doPlot,doPlot2,plotNum);

end


function [cumVelTout,cumVelXout,cumVelYout,cumVelPout,cumVelMout, ...
  cumLenTout,cumLenXout,cumLenYout,cumLenPout,cumLenMout ] = ...
  DoFtStuffCore(...
  eFieldX, eFieldY, ftVelX, ftVelY, ftLenX, ftLenY, ...
  dTime,omega,doPlot,doPlot2,plotNum)

[eFTX{1},eFTX{2},om] = TimeFT(eFieldX,dTime) ;
[eFTY{1},eFTY{2}]    = TimeFT(eFieldY,dTime) ;

%#ok<*AGROW>
for ii = 1:2
  ff = (-1)^(ii-1);
  eFTP{ii}   = 1/sqrt(2) * (eFTX{ii}   + ff * 1i * eFTY{ii});
  ftLenP{ii} = 1/sqrt(2) * (ftLenX{ii} + ff * 1i * ftLenY{ii});
  ftVelP{ii} = 1/sqrt(2) * (ftVelX{ii} + ff * 1i * ftVelY{ii});
  ftLenM{ii} = 1/sqrt(2) * (ftLenX{ii} - ff * 1i * ftLenY{ii});
  ftVelM{ii} = 1/sqrt(2) * (ftVelX{ii} - ff * 1i * ftVelY{ii});
  eFTM{ii}   = 1/sqrt(2) * (eFTX{ii}   - ff * 1i * eFTY{ii});
end

specFac = 1 / pi * dTime^2;
dom = om(2)-om(1);

% % yes, this one..  sum over positive and negative frequencies;
%  for either positive or negative, take real part of dot product
% SpecFun = @(eft,dipft) ...
%   -1/2 * real(conj(eft{1}) .* dipft{1}) .* specFac ...
%   -1/2 * real(conj(eft{2}) .* dipft{2}) .* specFac;

% % no, not this one
% SpecFun = @(eft,dipft) -1/2 * specFac * real( ...
SpecFun = @(eft,dipft) -1/2 * specFac * ( ...
  mean(eft{2} .* dipft{1} + eft{1} .* dipft{2}, 2) ...
  );

% SpecFun = @(eft,dipft) ...
%   -1 * real(eft{2} .* dipft{1}) .* specFac ;

% SpecFun = @(eft,dipft) ...
%  -1 * mean(conj(eft{1}) .* dipft{1},2) .* specFac;

%$$  EmitFun = @(x) 1/2 * specFac * mean( abs(x{1}.^2) + abs(x{2}.^2), 2);

EmitFun = @(x) specFac * mean( x{1}.*x{2}, 2);

specVelX = SpecFun(eFTX, ftVelX);
specLenX = SpecFun(eFTX, ftLenX);
specVelY = SpecFun(eFTY, ftVelY);
specLenY = SpecFun(eFTY, ftLenY);
specVelP = SpecFun(eFTP, ftVelP);
specLenP = SpecFun(eFTP, ftLenP);
specVelM = SpecFun(eFTM, ftVelM);
specLenM = SpecFun(eFTM, ftLenM);

emitVelX = EmitFun(ftVelX);
emitLenX = EmitFun(ftLenX);
emitVelY = EmitFun(ftVelY);
emitLenY = EmitFun(ftLenY);
emitVelP = EmitFun(ftVelP);
emitLenP = EmitFun(ftLenP);
emitVelM = EmitFun(ftVelM);
emitLenM = EmitFun(ftLenM);

CumSpecFun = @(x) cumsum(x) * dom ;

cumVelX = CumSpecFun(specVelX);
cumVelY = CumSpecFun(specVelY);
cumVelP = CumSpecFun(specVelP);
cumVelM = CumSpecFun(specVelM);
cumLenX = CumSpecFun(specLenX);
cumLenY = CumSpecFun(specLenY);
cumLenP = CumSpecFun(specLenP);
cumLenM = CumSpecFun(specLenM);

cumVelXout = cumVelX(end);
cumVelYout = cumVelY(end);
cumVelPout = cumVelP(end);
cumVelMout = cumVelM(end);
cumVelTout = cumVelXout + cumVelYout;
cumLenXout = cumLenX(end);
cumLenYout = cumLenY(end);
cumLenPout = cumLenP(end);
cumLenMout = cumLenM(end);
cumLenTout = cumLenXout + cumLenYout;

totEnd = real(cumLenX(end))+real(cumLenY(end));

FS = 22;
plotArgs = {'Linewidth',4,'MarkerSize',1.5};

if doPlot
  
  figure(500+plotNum)
  for ii=1:2
    % yunit = totEnd/dom;
    % yunit = max(abs(eFTX{1}.^2 + eFTY{1}.^2))*specFac;
    yunit = 1;
    
    plot(...
      ... %     om, real(specVelX), 'o', om, real(specLenX), '-', ...
      ... %     om, real(specVelY), 'o', om, real(specLenY), '-', ...
      ... %     om, real(specVelP), 'x', om, real(specLenP), '--', ...
      ... %     om, real(specVelM), 'x', om, real(specLenM), '--', ...
      om, real(specLenX)/yunit, '-', ...
      om, real(specLenY)/yunit, '-', ...
      om, real(specLenP)/yunit, '--', ...
      om, real(specLenM)/yunit, '--', ...
      om, real(specVelX)/yunit, 'o', ...
      om, real(specVelY)/yunit, 'o', ...
      om, real(specVelP)/yunit, 'x', ...
      om, real(specVelM)/yunit, 'x', ...
      plotArgs{:} )
    set(gca,'FontSize',FS)
    % title('StimE/Absorb work per unit \omega or t S(\omega or t) = E dot d/dt D')
    title('Absorption/StimEmission S = E dot d/dt D')
    xlabel('omega')
    % legend('Vx','Lx','Vy','Ly','V+','L+','V-','L-')
    legend('Lx','Ly','L+','L-','Vx','Vy','V+','V-')
    xlim([0.7 1.4]*omega)

    % yticks((-1:1)/10)
    % yticks([0,yunit]);
    % yticklabels({'0','1'})
    % % yticklabels({'0',sprintf('%4.2e',yunit)})

    yrr = floor(log10(diff(ylim)/3));
    trr =10^yrr;
    % yticks([-trr,0,trr]);
    yticks([-10*trr,-trr,-trr/10,-trr/100,0,trr/100,trr/10,trr,10*trr]);
    yticklabels({...
      sprintf('-10^{%2.0i}',yrr+1),...
      sprintf('-10^{%2.0i}',yrr),...
      sprintf('-10^{%2.0i}',yrr-1),...
      sprintf('-10^{%2.0i}',yrr-2),...
      '0',...
      sprintf('10^{%2.0i}',yrr-2)...
      sprintf('10^{%2.0i}',yrr-1)...
      sprintf('10^{%2.0i}',yrr)...
      sprintf('10^{%2.0i}',yrr+1)...
      })
    
    ylabel('S(\omega), work per unit \omega') %, units max |E(\omega)|^2')
    xlabel('Frequency \omega (inverse Hartree)')
    if ii==1
      figure(553+plotNum)
    else
      yrange = 1.5*max(abs(specVelP)/yunit);
      ylim([-yrange yrange]);
    end
  end
  
  drawnow
  
end

if doPlot2
  
  figure(499+plotNum)
  if 1==0
    loglog(...
      ...%       om, real(emitVelX), 'o', om, real(emitLenX),'-', ...
      ...%       om, real(emitVelY), 'o', om, real(emitLenY),'-', ...
      ...%       om, real(emitVelP), 'x', om, real(emitLenP),'--', ...
      ...%       om, real(emitVelM), 'x', om, real(emitLenM),'--', ...
      om, real(emitLenX),'-', ...
      om, real(emitLenY),'-', ...
      om, real(emitLenP),'--', ...
      om, real(emitLenM),'--', ...
      om, real(emitVelX), 'o', ...
      om, real(emitVelY), 'o', ...
      om, real(emitVelP), 'x', ...
      om, real(emitVelM), 'x', ...
      plotArgs{:} )
    set(gca,'FontSize',FS)  
    xlim([0.7 10]*omega)
    xlabel('\omega (Hartree)')
  else
    omega = 7.844;
    % lom = log(om+omega)/log(omega) - 1;
    % lom = (log(om/omega + 1/2) - log(1/2)) / log(2);
    lom = om/omega;
    lom(lom>1) = log(lom(lom>1))/log(2)+1;
    semilogy(...
      lom, real(emitLenX),'-', ...
      lom, real(emitLenY),'-', ...
      lom, real(emitLenP),'--', ...
      lom, real(emitLenM),'--', ...
      lom, real(emitVelX), 'o', ...
      lom, real(emitVelY), 'o', ...
      lom, real(emitVelP), 'x', ...
      lom, real(emitVelM), 'x', ...
      plotArgs{:} )
    set(gca,'FontSize',FS)  
    xlim([0 Inf])
    xlabel('x=log_2(\omega/\omega_0)+1 (x>1) or x=\omega/\omega_0 (x<1)')
  end  
  
  if ~doPlot  % it is meant that we are plotting emission.. temporary hack for logic
    title('|d^2/dt^2 D|^2')
  else
    title('|d/dt D|^2')
  end
  % legend('Vx','Lx','Vy','Ly','V+','L+','V-','L-')  
  legend('Lx','Ly','L+','L-','Vx','Vy','V+','V-')
  %$$ % ylim([10^-10 Inf])
  %$$ ylim([10^-14 10])
  %$$ ylim([10^-19 0.01])
  drawnow
end

if doPlot
  
  figure(501+plotNum)
  if 1==0
    mystuff = {...
      ... %       om, real(cumVelX), 'o', om, real(cumLenX), '-', ...
      ... %       om, real(cumVelY), 'o', om, real(cumLenY), '-', ...
      ... %       om, real(cumVelP), 'x', om, real(cumLenP), '--', ...
      ... %       om, real(cumVelM), 'x', om, real(cumLenM), '--', ...
      om, real(cumLenX), '-', ...
      om, real(cumLenY), '-', ...
      om, real(cumLenP), '--', ...
      om, real(cumLenM), '--', ...
      om, real(cumVelX), 'o', ...
      om, real(cumVelY), 'o', ...
      om, real(cumVelP), 'x', ...
      om, real(cumVelM), 'x', ...
      };
    % mylegend = {'Vx','Lx','Vy','Ly','V+','L+','V-','L-'};
    mylegend = { 'Lx','Ly','L+','L-','Vx','Vy','V+','V-' };
  else
    mystuff = {...
      ...% om, real(cumVelX+cumVelY), 'o', om, real(cumLenX+cumLenY),'-', ...
      ...% om, real(cumVelP+cumVelM), 'x', om, real(cumLenP+cumLenM),'--', ...
      om, real(cumLenX+cumLenY),'-', ...
      ...%om, real(cumLenP+cumLenM),'--', ...
      om, real(cumVelX+cumVelY), 'o', ...
      ...%om, real(cumVelP+cumVelM), 'x', ...
      };
    % mylegend = {'Vx/y','Lx/y','V+/-','L+/-'};
    % mylegend = {'Lx,Ly','L+/-','Vx/y','V+/-'};
    mylegend = {'Lx+Ly','Vx+Vy'};
  end
  plot(mystuff{:},  plotArgs{:});
  set(gca,'FontSize',FS)  
  if totEnd > 0
  yticks([0,totEnd]);
  yticklabels({'0',sprintf('%5.3e',totEnd)})
  % ytickformat('%4.2e')
  ylim(totEnd*[-0.1,1.2])
  end
  
  legend(mylegend{:},'Location','SouthEast');
  title('Total Absorb/StimEmission integral S(\omega) d\omega')
  xlim([0.7 1.4]*omega)
  ylabel('Energy (Hartree)')
  xlabel('Frequency \omega (inverse Hartree)')

  vec_pos = get(get(gca, 'YLabel'), 'Position');
  %set(get(gca, 'YLabel'), 'Position', vec_pos + [0.7 0 0]);
  set(get(gca, 'YLabel'), 'Position', vec_pos + [1.1 0 0]);

  drawnow
end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y = FDMult(x,fdvec,ord)
m = 2*ord+1;
if numel(fdvec)~=m
  error('ack')
end
n = numel(x);
if n < m
  error('oops toosmall')
end
y         = zeros(size(x));
for i     = 1:m
  iind    = mod((0:n-1) + (i - ord - 1), n) + 1 ;
  y       = y + fdvec(i) * x(iind);
end
end

function [fdvec,sdvec] = FDVec(ord)

n = 2*ord+1;
pts = -ord:ord;
fac = ord^((ord-1)/ord);
pts = pts / fac;

eqns = zeros(n,n);

for ieqn = 1:n
  eqns(ieqn,:) = pts.^(ieqn-1);
end

fdvec = eqns \ [0;1;zeros(n-2,1)];

sdvec = eqns \ [0;0;2;zeros(n-3,1)];

fdvec = fdvec / fac;
sdvec = sdvec / fac^2;

end

function [fdmat,sdmat] = FDMat(nPts,ord)

% n = 2*ord+1;

pts = -ord:ord;

[fdvec,sdvec] = FDVec(ord);

fdmat = zeros(nPts,nPts);
sdmat = zeros(nPts,nPts);
for ipt = 1:nPts
  iind = mod(pts + ipt - 1, nPts) + 1;
  fdmat(iind,ipt) = fdvec;
  sdmat(iind,ipt) = sdvec;
end

end



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


function [dTime,nTime] = TimeParms(omega,timeRes,duration,fdOrder,longerFac)
% there are nTime-1 steps

dTime = 1/omega * timeRes;

nTime = ceil(duration/dTime) + 1;

dTime = duration/(nTime-1);

% ensure differentiability at end!
nTime = nTime + fdOrder ;

% go for longer
nTime = ceil(nTime * longerFac);

end



function [fdmat,sdmat,ham0,potential,dipoleOpX,dipoleOpY,thetaVals] = ...
  CircleParts(problemOption,nPts,fdOrder,vMax,FLIPPARITY)

dTheta = 2*pi/nPts;
thetaVals  = ( 0.5:nPts-0.5 ) * dTheta;

[fdmat,sdmat] = FDMat(nPts,fdOrder);
fdmat = fdmat / dTheta;
sdmat = sdmat / dTheta^2;

if 1==0
  eit = exp(1i * thetaVals(:));
  enit = exp(-1i * thetaVals(:));
  fd2 = -1i/4 * ( eit(:) .* sdmat .* enit(:).' - enit(:) .* sdmat .* eit(:).' ) ;
  fd2 = real(fd2);
  disp('USING FDMAT FROM COMMUTATOR')
  fdmat = fd2;
end

%%%%%%%%%

[potential] = GetPotential(problemOption,thetaVals,vMax,FLIPPARITY);

kemat = -1/2 * sdmat ;

ham0 = kemat + diag(potential);

dipoleOpX = diag(cos(thetaVals));
dipoleOpY = diag(sin(thetaVals));

end


% omega just for plotting
function [eigVects,eigVals] = CircleEigs(problemOption,nPts,fdOrder,vMax,iSel,doOut,FLIPPARITY,omega)

[~,~,ham0,potential,dipoleOpX,dipoleOpY,thetaVals] = CircleParts(problemOption,nPts,fdOrder,vMax,FLIPPARITY);

ham0 = (ham0 + ham0')/2;
[eigVects0,eigVals0] = eig(ham0);

eigVals0 = diag(eigVals0);

[~,sord] = sort(eigVals0);
eigVals0 = eigVals0(sord);
eigVects0 = eigVects0(:,sord);

% figure(100)
% plot(eigVals0(1:10),'o','LineWidth',4)
% title('eigVals')

eigVects = eigVects0(:,iSel);
eigVals = eigVals0(iSel);

if doOut
  disp('eigVals')
  disp([iSel',eigVals,eigVals-eigVals(1)])
end

eigVectsSq = eigVects.^2;

norms = sqrt(diag(eigVects'*eigVects));
eigVects = eigVects ./ norms(:)';

norms = sqrt(diag(eigVectsSq'*eigVectsSq));
eigVectsSq = eigVectsSq ./ norms(:)';

sqOvlMat = eigVectsSq'*eigVectsSq;

if doOut
  figure(101)
  fac = 10/omega;
  tfac = 1; %180/pi;
  plot(thetaVals*tfac,potential/omega - eigVals(1)/omega,'-',...
    thetaVals*tfac,fac*eigVects + eigVals'/omega - eigVals(1)/omega,'.',...
    'Linewidth',4,'MarkerSize',8)
  set(gca,'FontSize',16)
  xlim([0 2*pi*tfac])
  xlabel('theta')
  ylabel('V(theta) / omega')
  title('PES and Eigs 1-4,6,7,10,11')
  drawnow
end

if doOut
  xOpMat = eigVects' * dipoleOpX * eigVects;
  yOpMat = eigVects' * dipoleOpY * eigVects;
  
  sqOpMat = xOpMat.^2 + yOpMat.^2;
  
  disp(['Selected States: ' num2str(iSel)])
  
  disp('Sq Ovl Mat')
  disp(sqOvlMat)
  disp('Sq Dipole Mat')
  disp(sqOpMat)
  
  disp('X Dipole Mat')
  disp(xOpMat)
  disp('Y Dipole Mat')
  disp(yOpMat)
end

end




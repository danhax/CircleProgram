%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   DOMCKE WAVE MIXING   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#ok<*AGROW>

function CircleDomcke(DOAVERAGE, nPhase, nKp, ...
  problemOption,nPts,fdOrder,vMax,omega,duration,...
  strXS,strXC,strYS,strYC,oStr, ENV_OPT, ...
  iStart,iSel,       longerFac,timeRes,propMode,modPlot, ...
  DOFTLAST, DOFTK, DOCOMPLEXD)

nWay = 2;
nGauge = 2;
nPol = 2;
nDer = 3;

[dTime,nTime] = TimeParms(omega,timeRes,duration,fdOrder,longerFac);

nTimes = zeros(1,nWay);
nTimes(1) = nTime - fdOrder;
nTimes(2) = nTime ;

  function [eFieldCell,aFieldCell,ddtDipCell] = DoProp(phasePos,phaseNeg,phaseMed)
    doPlot = false;
    eFieldCell = cell(nWay,1);
    aFieldCell = cell(nWay,1);
    ddtDipCell = cell(nWay,1);
    for iiway = 1:nWay
      eFieldCell{iiway} = zeros(nTimes(iiway),nPol);
      aFieldCell{iiway} = zeros(nTimes(iiway),nPol);
      ddtDipCell{iiway} = zeros(nTimes(iiway),nDer,nGauge,nPol);
    end
    [~,~,~,~, ...
      eFieldCell{1}(:,1),   eFieldCell{1}(:,2), aFieldCell{1}(:,1),   aFieldCell{1}(:,2), ...
      ddtDipCell{1}(:,:,1,1), ddtDipCell{1}(:,:,1,2), ddtDipCell{1}(:,:,2,1), ddtDipCell{1}(:,:,2,2),...
      eFieldCell{2}(:,1),   eFieldCell{2}(:,2), aFieldCell{2}(:,1),   aFieldCell{2}(:,2), ...
      ddtDipCell{2}(:,:,1,1), ddtDipCell{2}(:,:,1,2), ddtDipCell{2}(:,:,2,1), ddtDipCell{2}(:,:,2,2),...
      ] = CircleProp(...
      problemOption,nPts,fdOrder,vMax,omega,duration,...
      strXS,strXC,strYS,strYC,oStr,phasePos,phaseNeg,phaseMed,ENV_OPT,...
      iStart,iSel,     doPlot,longerFac,timeRes,propMode,modPlot);
  end

phases = (0:nPhase-1)*2*pi/nPhase;
kPhases = (0:nKp-1)*2*pi/nKp;

% REGULAR DOMCKE

eFieldCell = cell(nWay,1);
aFieldCell = cell(nWay,1);
ddtDip0Cell = cell(nWay,1);
for iway = 1:nWay
  eFieldCell{iway} = zeros(nTimes(iway),nPol,nKp);
  aFieldCell{iway} = zeros(nTimes(iway),nPol,nKp);
  ddtDip0Cell{iway} = zeros(nTimes(iway),nDer,nGauge,nPol,nKp);
end

for kphase = 1:nKp
  phase_k = kPhases(kphase);
  [efc,afc,ddc] = DoProp(phase_k,phase_k,phase_k);
  for iway = 1:nWay
    eFieldCell{iway}(:,:,kphase) = efc{iway};
    aFieldCell{iway}(:,:,kphase) = afc{iway};
    ddtDip0Cell{iway}(:,:,:,:,kphase) = ddc{iway};
  end
end

iWay = 2;
plotNum = (iWay-1)*1000;

DomckeCore(DOAVERAGE,nTimes(iWay),dTime,omega,nDer,nGauge,nPol,nKp,eFieldCell{iWay},aFieldCell{iWay},ddtDip0Cell{iWay}, plotNum ) ;

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
    ddtDipCell{iway} = zeros(nTimes(iway),nDer,nGauge,nPol,nKp,nPhase,nPhase2);
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
        [~,~,ddc] = DoProp(phase_i,phase_j,phase_k);
        
        for iway = 1:nWay
          if SUBSTUFF~=0
            ddtDipCell{iway}(:,:,:,:,kphase,iphase,jphase) = ddc{iway} - ddtDip0Cell{iway}(:,:,:,:,kphase);
          else
            ddtDipCell{iway}(:,:,:,:,kphase,iphase,jphase) = ddc{iway};
          end
        end
      end
    end
  end
  
  ComplexDomckeCore(DOAVERAGE, DOFTK, DOFTLAST, ...
    nTimes(iWay),dTime,omega,nDer,nGauge,nPol,nPhase,nKp,...
    eFieldCell{iWay}, aFieldCell{iWay}, ddtDipCell{iWay}, plotNum );
  
end

end



function ftDdtDip0ft = DomckeCore(DOAVERAGE, nTime,dTime,omega,nDer,nGauge,nPol,nKp, eField, aField, ddtDip0In,plotNum)

eFieldX = reshape(eField(:,1,:),nTime,nKp);
eFieldY = reshape(eField(:,2,:),nTime,nKp);

aFieldX = reshape(aField(:,1,:),nTime,nKp);
aFieldY = reshape(aField(:,2,:),nTime,nKp);

ddtDip0 = ddtDip0In;

ftDdtDip0 = FtDip(dTime,ddtDip0);

if DOAVERAGE
  kList = 1:nKp;
else
  kList = 1;
end
nKl = numel(kList);

for ispec=1:2
  ftLenX0{ispec} = reshape(ftDdtDip0{ispec}(:,:,1,1,kList),[],nDer,nKl);
  ftLenY0{ispec} = reshape(ftDdtDip0{ispec}(:,:,1,2,kList),[],nDer,nKl);
  ftVelX0{ispec} = reshape(ftDdtDip0{ispec}(:,:,2,1,kList),[],nDer,nKl);
  ftVelY0{ispec} = reshape(ftDdtDip0{ispec}(:,:,2,2,kList),[],nDer,nKl);
end

ActualResult = @() ...
  DoFtStuffCore(eFieldX(:,kList),eFieldY(:,kList),aFieldX(:,kList),aFieldY(:,kList),...
  ftVelX0,ftVelY0,ftLenX0,ftLenY0,dTime,omega,true,true,plotNum);
save('ActualResult','ActualResult')
ActualResult();

if DOAVERAGE
  
  ddtDip0_2 = zeros(nTime,nDer,nGauge,nPol,nKp,nKp);
  for ik = 1:nKp
    ddtDip0_2(:,:,:,:,:,ik) = ddtDip0(:,:,:,:,mod(ik+(1:nKp)-2,nKp)+1);
  end
  
  ftDdtDip0_2 = FtDip(dTime,ddtDip0_2);
  ftDdtDip0ft = domckeFt(ftDdtDip0_2,6);
  
  for ispec=1:2
    ftLenX0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,:,1,1,:,:),[],nDer,nKp,nKp);
    ftLenY0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,:,1,2,:,:),[],nDer,nKp,nKp);
    ftVelX0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,:,2,1,:,:),[],nDer,nKp,nKp);
    ftVelY0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,:,2,2,:,:),[],nDer,nKp,nKp);
  end
  
  GET = @(a,ik) {a{1}(:,:,:,ik),a{2}(:,:,:,ik)};
  
else
  
  ftDdtDip0ft = domckeFt(ftDdtDip0,5);
  
  for ispec=1:2
    ftLenX0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,:,1,1,:),[],nDer,nKp);
    ftLenY0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,:,1,2,:),[],nDer,nKp);
    ftVelX0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,:,2,1,:),[],nDer,nKp);
    ftVelY0ft{ispec} = reshape(ftDdtDip0ft{ispec}(:,:,2,2,:),[],nDer,nKp);
  end
  
  GET = @(a,ik) {a{1}(:,:,ik),a{2}(:,:,ik)};
  
end

DomckeResult = @(ik) ...
  DoFtStuffCore(eFieldX(:,kList),eFieldY(:,kList),aFieldX(:,kList),aFieldY(:,kList),...
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
  x(:,2:end,:)             = flip(x(:,2:end,:),ind);
elseif ind == 4
  x(:,:,:,2:end,:)         = flip(x(:,:,:,2:end,:),ind);
elseif ind == 5
  x(:,:,:,:,2:end,:)       = flip(x(:,:,:,:,2:end,:),ind);
elseif ind == 6
  x(:,:,:,:,:,2:end,:)     = flip(x(:,:,:,:,:,2:end,:),ind);
elseif ind == 7
  x(:,:,:,:,:,:,2:end,:)   = flip(x(:,:,:,:,:,:,2:end,:),ind);
elseif ind == 8
  x(:,:,:,:,:,:,:,2:end,:) = flip(x(:,:,:,:,:,:,:,2:end,:),ind);
else
  error('not supported')
end
end


function ComplexDomckeCore(DOAVERAGE, DOFTK, DOFTLAST, nTime,dTime,omega,...
  nDer,nGauge,nPol,nPhase,nKp,...
  eField, aField, ddtDipIn, plotNum )

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

aFieldX = reshape(aField(:,1,kList),nTime,nKl);
aFieldY = reshape(aField(:,2,kList),nTime,nKl);

clear kList;

if DOFTK
  if DOAVERAGE
    ddtDip = zeros(nTime,nDer,nGauge,nPol,nKp,nKp,nPhase,nPhase2);
    for ik = 1:nKp
      pk = mod(ik+(1:nKp)-2,nKp)+1 ;
      ddtDip(:,:,:,:,:,ik,:,:) = ddtDipIn(:,:,:,:,pk,:,:);
    end
  else
    ddtDip = reshape(ddtDipIn,nTime,nDer,nGauge,nPol,1,nKp,nPhase,nPhase2);
  end
else
  if DOAVERAGE
    ddtDip = reshape(ddtDipIn,nTime,nDer,nGauge,nPol,nKp,1,nPhase,nPhase2);
  else
    if nKp ~= 1
      error('what?')
    end
    ddtDip = reshape(ddtDipIn,nTime,nDer,nGauge,nPol,1,1,nPhase,nPhase2);
  end
end

ftDdtDip = FtDip(dTime,ddtDip);

if DOFTLAST
  if DOFTK
    ftDdtDipft = domckeFt(domckeFt(domckeFt(ftDdtDip,6),7),8);
  else
    ftDdtDipft = domckeFt(domckeFt(ftDdtDip,7),8);
  end
else
  if DOFTK
    ftDdtDipft = domckeFt(domckeFt(ftDdtDip,6),7);
  else
    ftDdtDipft = domckeFt(ftDdtDip,7);
  end
end

for ispec=1:2
  ftLenXft{ispec} = reshape(ftDdtDipft{ispec}(:,:,1,1,:),[],nDer,nKl,nFt,nPhase,nPhase2);
  ftLenYft{ispec} = reshape(ftDdtDipft{ispec}(:,:,1,2,:),[],nDer,nKl,nFt,nPhase,nPhase2);
  ftVelXft{ispec} = reshape(ftDdtDipft{ispec}(:,:,2,1,:),[],nDer,nKl,nFt,nPhase,nPhase2);
  ftVelYft{ispec} = reshape(ftDdtDipft{ispec}(:,:,2,2,:),[],nDer,nKl,nFt,nPhase,nPhase2);
end

GET = @(a,kp,ip,jp) {a{1}(:,:,:,kp,ip,jp),a{2}(:,:,:,kp,ip,jp)};

CDomckeResult = @(kp,ip,jp) ...
  DoFtStuffCore(eFieldX,eFieldY,aFieldX,aFieldY, ...
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


% function [popVel,popLen,workVel,workLen,...
%   eFieldX2, eFieldY2, cur2LenX, cur2LenY, cur2VelX, cur2VelY,...
%   eFieldX, eFieldY,  curLenX,  curLenY,  curVelX,  curVelY ] = ...
function [popVel,popLen,workVel,workLen,...
  eFieldX2, eFieldY2, aFieldX2, aFieldY2, ...
  lenX2, lenY2, velX2, velY2,...
  eFieldX, eFieldY, aFieldX, aFieldY, ...
  lenX, lenY, velX, velY ] = ...
  CircleProp( ...
  problemOption,nPts,fdOrder,vMax,omega,duration,...
  strXS,strXC,strYS,strYC,oStr,phasePos,phaseNeg,phaseMed,ENV_OPT,...
  iStart,iSel,     doPlot,longerFac,timeRes,propMode,modPlot)

DOPROJ        = 1==0 ;
nQuad1        = 5;
DOPARITY      = 1==0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DOCONJ        = false;
if phasePos ~= phaseNeg
  DOCONJ      = true;
end

startVec      = CircleEigs(problemOption,nPts,fdOrder,vMax,[iStart,iSel],doPlot,false,omega);
eigVects      = startVec(:,2:end);
startVec      = startVec(:,1);

startVecp     = [];
eigVectsp     = [];
if DOPARITY
  startVecp   = CircleEigs(problemOption,nPts,fdOrder,vMax,[iStart,iSel],false,true,omega);
  eigVectsp   = startVecp(:,2:end);
  startVecp   = startVecp(:,1);
end

clear iStart iSel

nEig          = size(eigVects,2);

stateProj     = [];
stateProjp    = [];
if DOPROJ
  stateProj   = eigVects * eigVects';
  stateProjp  = eigVectsp * eigVectsp';
end

[fdmat,sdmat,ham0,potential,cosThetaDiag,sinThetaDiag,thetaVals,gradpot,hesspot] ...
  = CircleParts(problemOption,nPts,fdOrder,vMax,false);

ham0p         = [];
if DOPARITY
  [~,~,ham0p] = CircleParts(problemOption,nPts,fdOrder,vMax,true);
end

if propMode==0   % do sparse for Crank-Nicholson...
  fdmat       = sparse( fdmat );
  ham0        = sparse( ham0 );
  ham0p       = sparse( ham0p );
  EYE         = speye( nPts );
else
  fdmat       = full( fdmat );
  ham0        = full( ham0 );
  ham0p       = full( ham0p );
  EYE         = eye( nPts );
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

startTimes    = (0:(nTime-2))' * dTime;
% halfTimes   = (0.5:(nTime-1.5)) * dTime;
endTimes      = (1:(nTime-1))' * dTime;
dataTimes     = (0:(nTime-1))' * dTime;

%%%%%%%%%%%%%%%%%    Fields    %%%%%%%%%%%%%%%%%%%%%%%%%%

GetFieldsHere = @(intimes,doplot,shiftPos,shiftNeg,~) ...
  GetFields3B(intimes,strXS,strXC,strYS,strYC,oStr,omega,duration,shiftPos,shiftNeg, doplot,ENV_OPT);

[aFieldX,aFieldY,eFieldX,eFieldY] = ...
  GetFieldsHere(dataTimes,doPlot,phasePos,phaseNeg,phaseMed) ;

AvalX    = @(t) GetOut([1,2],t,@(x)GetFieldsHere(x,false,phasePos,phaseNeg,phaseMed));
AvalY    = @(t) GetOut(2,t,@(x)GetFieldsHere(x,false,phasePos,phaseNeg,phaseMed));
EvalX    = @(t) GetOut(3,t,@(x)GetFieldsHere(x,false,phasePos,phaseNeg,phaseMed));
EvalY    = @(t) GetOut(4,t,@(x)GetFieldsHere(x,false,phasePos,phaseNeg,phaseMed));
DdtEvalX = @(t) GetOut(5,t,@(x)GetFieldsHere(x,false,phasePos,phaseNeg,phaseMed));
DdtEvalY = @(t) GetOut(6,t,@(x)GetFieldsHere(x,false,phasePos,phaseNeg,phaseMed));

%%%%%%%%%%%%%%%%%%  Hamiltonians   %%%%%%%%%%%%%%%%%

dipoleLenX =  cosThetaDiag(:) ;
dipoleLenY =  sinThetaDiag(:);

dipoleVelX = -sinThetaDiag(:);
dipoleVelY =  cosThetaDiag(:);

%$ dipoleLenXDer = -sinThetaDiag(:) ;    % same as dipoleVelX
%$ dipoleLenYDer =  cosThetaDiag(:);     % same as dipoleVelY

LenE          = @(t) ( EvalX(t) * dipoleLenX + EvalY(t) * dipoleLenY);
LenEder       = @(t) ( EvalX(t) * dipoleVelX + EvalY(t) * dipoleVelY);  % e.g. dipoleVelY = dipoleLenYDer

DdtLenE       = @(t) ( DdtEvalX(t) * dipoleLenX + DdtEvalY(t) * dipoleLenY);
DdtLenEder    = @(t) ( DdtEvalX(t) * dipoleVelX + DdtEvalY(t) * dipoleVelY);

LenHam     = @(t) ham0  + diag(LenE(t));
LenHamp    = @(t) ham0p + diag(LenE(t));

dipoleVelXDer = -cosThetaDiag(:);
dipoleVelYDer = -sinThetaDiag(:);

%$$ HamVel = @(aval) (1/2) * (1i*fdmat + aval*dipoleOpVel)' * (1i*fdmat + aval*dipoleOpVel) + diag(potential);

VelA          = @(t) ( AvalX(t) * dipoleVelX + AvalY(t) * dipoleVelY);
DdtVelA       = @(t) ( EvalX(t) * dipoleVelX + EvalY(t) * dipoleVelY);
DdtVelAder    = @(t) ( EvalX(t) * dipoleVelXDer + EvalY(t) * dipoleVelYDer);
DdtSqVelA    = @(t)(DdtEvalX(t) * dipoleVelX + DdtEvalY(t) * dipoleVelY);

dOVXX         = dipoleVelX.^2;
dOVYY         = dipoleVelY.^2;
dOVcross      = 2 * dipoleVelX .* dipoleVelY ;

% VelASq = @(t) AvalX(t)^2 * dOVXX + AvalY(t)^2 * dOVYY + AvalX(t)*AvalY(t) * dOVcross ;
  function y = VelASq0(t,AvalX,dOVXX,dOVYY,dOVcross)
    [avx,avy] = AvalX(t);
    y = avx^2 * dOVXX + avy^2 * dOVYY + avx * avy * dOVcross;
  end
VelASq = @(t) VelASq0(t,AvalX,dOVXX,dOVYY,dOVcross);
clear VelASq0;

DdtVelASq = @(t) ...
  2 * EvalX(t) * AvalX(t) * dOVXX ...
  + 2 * EvalY(t) * AvalY(t) * dOVYY ...
  + ( EvalX(t)*AvalY(t) + AvalX(t)*EvalY(t) ) * dOVcross ;

VelHam  = @(t) ham0  + 1i/2 * ( VelA(t) .* fdmat + fdmat .* VelA(t).' ) + 1/2 * diag(VelASq(t));
VelHamp = @(t) ham0p + 1i/2 * ( VelA(t) .* fdmat + fdmat .* VelA(t).' ) + 1/2 * diag(VelASq(t));
VelSD   = @(t) sdmat - 1i   * ( VelA(t) .* fdmat + fdmat .* VelA(t).' ) -       diag(VelASq(t));

% DdtVelHam = @(t) 1i/2 * ( DdtVelA(t) .* fdmat + fdmat .* DdtVelA(t).' ) + 1/2 * diag(DdtVelASq(t));

DdtVelSD = @(t) -1i * ( DdtVelA(t) .* fdmat + fdmat .* DdtVelA(t).' )  - diag(DdtVelASq(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DIPOLE D, (ddt D) DER, (ddt^2 D) ACC, (ddt^3 D) JERK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DOPARITY
  error('cannot DOPARITY any more with acceleration, need to program')
end

if 1==1    % done with commutators, no algebra

  %VelP    = @(t) (1i*fdmat + diag(VelA(t)));

  %VelJX = @(t) 1/2*(dipoleVelX(:) .* VelP(t) + VelP(t) .* dipoleVelX(:).');
  %VelJY = @(t) 1/2*(dipoleVelY(:) .* VelP(t) + VelP(t) .* dipoleVelY(:).');

  LenJX =  -1i/2*(dipoleVelX .* fdmat + fdmat .* dipoleVelX.');
  LenJY =  -1i/2*(dipoleVelY .* fdmat + fdmat .* dipoleVelY.');
  
  VelJX = @(t) LenJX - diag(VelA(t) .* dipoleVelX);
  VelJY = @(t) LenJY - diag(VelA(t) .* dipoleVelY);

  DdtVelJX = @(t) -1 * diag(DdtVelA(t).*dipoleVelX) ;
  DdtVelJY = @(t) -1 * diag(DdtVelA(t).*dipoleVelY) ;
  DdtSqVelJX = @(t) -1 * diag(DdtSqVelA(t).*dipoleVelX) ;
  DdtSqVelJY = @(t) -1 * diag(DdtSqVelA(t).*dipoleVelY) ;

  VelCurX = @(vc,t,v) vc' * VelJX(t) * v;
  VelCurY = @(vc,t,v) vc' * VelJY(t) * v;
  
  LenFieldX = @(vc,t,v) sum( conj(vc) .* VelA(t) .* dipoleVelX .* v, 1);
  LenFieldY = @(vc,t,v) sum( conj(vc) .* VelA(t) .* dipoleVelY .* v, 1);
  VelFieldX = @(vc,t,v) sum( conj(vc) .* VelA(t) .* dipoleVelX .* v, 1);
  VelFieldY = @(vc,t,v) sum( conj(vc) .* VelA(t) .* dipoleVelY .* v, 1);

  %%%%%%   ACCELERATION, JERK:  LENGTH
  %
  % % DdtVelA(t) and LenEder(t) are the same
  
  %   LenAXop = @(t) -1i*(LenJX * LenHam(t) - LenHam(t) * LenJX) ;

  %$ LenJX__sdmat = LenJX * sdmat - sdmat * LenJX;
  %$ LenJY__sdmat = LenJY * sdmat - sdmat * LenJY;
  %$ LenJX__pot = LenJX .* potential.' - potential .* LenJX;
  %$ LenJY__pot = LenJY .* potential.' - potential .* LenJY;
  %$ LenJX__LenE = @(t) LenJX .* LenE(t).' - LenE(t) .* LenJX ;
  %$ LenJY__LenE = @(t) LenJY .* LenE(t).' - LenE(t) .* LenJY ;
  
  LenJX__potdiag = -1i/2 * ( 2 * dipoleVelX .* gradpot );
  LenJY__potdiag = -1i/2 * ( 2 * dipoleVelY .* gradpot );
  LenJX__LenEdiag = @(t) -1i/2 * ( 2 * dipoleVelX.*LenEder(t) );
  LenJY__LenEdiag = @(t) -1i/2 * ( 2 * dipoleVelY.*LenEder(t) );
  DdtLenJX__LenE = @(t) -1i/2 * diag( 2 * dipoleVelX.*DdtLenEder(t) );
  DdtLenJY__LenE = @(t) -1i/2 * diag( 2 * dipoleVelY.*DdtLenEder(t) );
  % -2 DD F^2 + DDDD - 2 F^2 DD
  LenJX__sdmat = -1i/2 * ...
    ( -2*(dipoleVelXDer .* sdmat + sdmat .* dipoleVelXDer.') + diag(dipoleLenX) );
  LenJY__sdmat = -1i/2 * ...
    ( -2*(dipoleVelYDer .* sdmat + sdmat .* dipoleVelYDer.') + diag(dipoleLenY) );
  
  LenAXop = @(t) -1i * ( -1/2 * LenJX__sdmat + diag(LenJX__LenEdiag(t) + LenJX__potdiag) );
  LenAYop = @(t) -1i * ( -1/2 * LenJY__sdmat + diag(LenJY__LenEdiag(t) + LenJY__potdiag) );

  LenAccX   = @(vc,t,v) (vc' * LenAXop(t) * v);
  LenAccY   = @(vc,t,v) (vc' * LenAYop(t) * v);
  
  %   LenJrXop = @(t) -1i*(LenAXop(t) * LenHam(t) - LenHam(t) * LenAXop(t)) ...
  %     -1i*DdtLenJX__LenE(t) ;
  %   LenJrYop = @(t) -1i*(LenAYop(t) * LenHam(t) - LenHam(t) * LenAYop(t)) ...
  %     -1i*DdtLenJY__LenE(t) ;
  %
  
  LenJrXop0 = -1i*(-1i)*( ...
    -1/2 * ( LenJX__sdmat .* potential.' - potential .* LenJX__sdmat) ...
    +1/4 * ( LenJX__sdmat    * sdmat - sdmat  * LenJX__sdmat) ...
    -1/2 * ( LenJX__potdiag     .* sdmat - sdmat .* LenJX__potdiag.' ) ...
    ) ;
  LenJrXop = @(t) LenJrXop0 -1i*(-1i)*( ...
    -1/2 * ( LenJX__sdmat .* LenE(t).' - LenE(t) .* LenJX__sdmat) ...
    -1/2 * ( LenJX__LenEdiag(t) .* sdmat - sdmat .* LenJX__LenEdiag(t).' ) ...
    ) ...
    -1i*DdtLenJX__LenE(t) ;
  
  LenJrYop0 = -1i*(-1i)*( ...
    -1/2 * ( LenJY__sdmat .* potential.' - potential .* LenJY__sdmat) ...
    +1/4 * ( LenJY__sdmat    * sdmat - sdmat  * LenJY__sdmat) ...
    -1/2 * ( LenJY__potdiag     .* sdmat - sdmat .* LenJY__potdiag.' ) ...
    ) ;
  LenJrYop = @(t) LenJrYop0 -1i*(-1i)*( ...
    -1/2 * ( LenJY__sdmat .* LenE(t).' - LenE(t) .* LenJY__sdmat) ...
    -1/2 * ( LenJY__LenEdiag(t) .* sdmat - sdmat .* LenJY__LenEdiag(t).' ) ...
    ) ...
    -1i*DdtLenJY__LenE(t) ;
  
  
  %{
  LenJrXop1 = @(t) -1i*(-1i)*( ...
    -1/2 * ( LenJX__sdmat .* LenE(t).' - LenE(t) .* LenJX__sdmat) ...
    -1/2 * ( LenJX__sdmat .* potential.' - potential .* LenJX__sdmat) ...
    +1/4 * ( LenJX__sdmat    * sdmat - sdmat  * LenJX__sdmat) ...
    -1/2 * ( LenJX__potdiag     .* sdmat - sdmat .* LenJX__potdiag.' ) ) ...
    ;
  LenJrXop2 = @(t) -1i*(-1i)*( ...
    -1/2 * ( LenJX__LenEdiag(t) .* sdmat - sdmat .* LenJX__LenEdiag(t).' ) ) ...
    -1i*DdtLenJX__LenE(t) ...
  ;

  LenJrYop1 = @(t) -1i*(-1i)*( ...
    -1/2 * ( LenJY__sdmat .* LenE(t).' - LenE(t) .* LenJY__sdmat) ...
    -1/2 * ( LenJY__sdmat .* potential.' - potential .* LenJY__sdmat) ...
    +1/4 * ( LenJY__sdmat    * sdmat - sdmat  * LenJY__sdmat) ...
    -1/2 * ( LenJY__potdiag     .* sdmat - sdmat .* LenJY__potdiag.' ) ) ...
    ;
  LenJrYop2 = @(t) -1i*(-1i)*( ...
    -1/2 * ( LenJY__LenEdiag(t) .* sdmat - sdmat .* LenJY__LenEdiag(t).' ) ) ...
    -1i*DdtLenJY__LenE(t) ...
  ;
  %}
  %{
  LenJrXop1 = @(t) -1i*(-1i)*( 0 + ...
    -1/2 * ( LenJX__sdmat .* LenE(t).' - LenE(t) .* LenJX__sdmat) ...
    -1/2 * ( LenJX__sdmat .* potential.' - potential .* LenJX__sdmat) ...
    +1/4 * ( LenJX__sdmat    * sdmat - sdmat  * LenJX__sdmat) ...
    -1/2 * ( LenJX__potdiag     .* sdmat - sdmat .* LenJX__potdiag.' ) ...
    ) ...
    ;
  LenJrXop2 = @(t) -1i*(-1i)*( 0 + ...
    -1/2 * ( LenJX__LenEdiag(t) .* sdmat - sdmat .* LenJX__LenEdiag(t).' ) ...
    ) ...
    -1i*DdtLenJX__LenE(t) ...
  ;

  LenJrYop1 = @(t) -1i*(-1i)*( 0 + ...
    -1/2 * ( LenJY__sdmat .* LenE(t).' - LenE(t) .* LenJY__sdmat) ...
    -1/2 * ( LenJY__sdmat .* potential.' - potential .* LenJY__sdmat) ...
    +1/4 * ( LenJY__sdmat    * sdmat - sdmat  * LenJY__sdmat) ...
    -1/2 * ( LenJY__potdiag     .* sdmat - sdmat .* LenJY__potdiag.' ) ...
    ) ...
    ;
  LenJrYop2 = @(t) -1i*(-1i)*( 0 + ...
    -1/2 * ( LenJY__LenEdiag(t) .* sdmat - sdmat .* LenJY__LenEdiag(t).' ) ...
    ) ...
    -1i*DdtLenJY__LenE(t) ...
    ;

  % LenJrXop = @(t) LenJrXop1(t) + LenJrXop2(t);
  % LenJrYop = @(t) LenJrYop1(t) + LenJrYop2(t);
  
  LenJrXop = @(t) LenJrXop1(t) ;
  LenJrYop = @(t) LenJrYop1(t) ;
  %}

  
  %%%%% VELOCITY
  
  % VelJX__VelHam =  @(t) VelJX(t) * VelHam(t) - VelHam(t) * VelJX(t) ;
  % VelJY__VelHam =  @(t) VelJY(t) * VelHam(t) - VelHam(t) * VelJY(t) ;
  
  % -2 DD F^2 + DDDD - 2 F^2 DD
  VelJX__VelHam = @(t) diag(LenJX__potdiag) + (-1/2) * -1i/2 * ...
    ( -2*(dipoleVelXDer .* VelSD(t) + VelSD(t) .* dipoleVelXDer.') + diag(dipoleLenX) );
  VelJY__VelHam = @(t) diag(LenJY__potdiag) + (-1/2) * -1i/2 * ...
    ( -2*(dipoleVelYDer .* VelSD(t) + VelSD(t) .* dipoleVelYDer.') + diag(dipoleLenY) );
  
  % there is cancellation which has not been accounted for.
  % note the form of LenAXop, LenAYop with diagonal time-dependent part
  VelAXop = @(t) -1i*VelJX__VelHam(t) + DdtVelJX(t) ;
  VelAYop = @(t) -1i*VelJY__VelHam(t) + DdtVelJY(t) ;
  
  DdtVelAXop = @(t) DdtSqVelJX(t) + -1/2 * ...   %-1i * (-1/2) * -1i/2 * (-2) * ...
    (dipoleVelXDer .* DdtVelSD(t) + DdtVelSD(t) .* dipoleVelXDer.');
  DdtVelAYop = @(t) DdtSqVelJY(t) + -1/2 * ...   %-1i * (-1/2) * -1i/2 * (-2) * ...
    (dipoleVelYDer .* DdtVelSD(t) + DdtVelSD(t) .* dipoleVelYDer.');
  
  VelAccX   = @(vc,t,v) (vc' * VelAXop(t) * v);
  VelAccY   = @(vc,t,v) (vc' * VelAYop(t) * v);
  
  VelJrXop = @(t) -1i*(VelAXop(t) * VelHam(t) - VelHam(t) * VelAXop(t)) ...
    + DdtVelAXop(t);
  VelJrYop = @(t) -1i*(VelAYop(t) * VelHam(t) - VelHam(t) * VelAYop(t)) ...
    + DdtVelAYop(t);

  LenJerkX   = @(vc,t,v) (vc' * LenJrXop(t) * v);
  LenJerkY   = @(vc,t,v) (vc' * LenJrYop(t) * v);
  VelJerkX   = @(vc,t,v) (vc' * VelJrXop(t) * v);
  VelJerkY   = @(vc,t,v) (vc' * VelJrYop(t) * v);

  VelCenX = @(vc,t,v) 0;
  VelCenY = @(vc,t,v) 0;

else
  error('oops')
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear dipoleLenX dipoleLenY dipoleVelX dipoleVelY VelJ VelJ2 Den VelP VelPc VelA VelASq
clear dOVXX dOVYY dOVcross
clear gradpot hesspot 
clear AccDiagXDer AccDiagYDer AccDiagX AccDiagY

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
clear LenHam VelHam LenHamp VelHamp

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
PropAll = @(hl,hv,  vl,vv,vlc,vvc)  ...
  PropAllCore(hl,hv,  vl,vv,vlc,vvc,  propMode, DOCONJ, EYE, dTime, PropIt0, false);
PropAll0 = @(hl,hv,  vl,vv,vlc,vvc)  ...
  PropAllCore(hl,hv,  vl,vv,vlc,vvc,  propMode, DOCONJ, EYE, dTime, PropIt0, true);

PropIt0p = @(v) U0p * v ;
PropAllp = @(hl,hv,  vl,vv,vlc,vvc)  ...
  PropAllCore(hl,hv,  vl,vv,vlc,vvc,  propMode, DOCONJ, EYE, dTime, PropIt0p, false);
PropAll0p = @(hl,hv,  vl,vv,vlc,vvc)  ...
  PropAllCore(hl,hv,  vl,vv,vlc,vvc,  propMode, DOCONJ, EYE, dTime, PropIt0p, true);

clear PropIt0 PropIt0p U0 U0p PropAllCore propMode ham0 ham0p

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

  function [ dipLenX,dipLenY,dipVelX,dipVelY, ...
      curLenX,curLenY,curVelX,curVelY, ...
      accLenX,accLenY,accVelX,accVelY, ...
      jerLenX,jerLenY,jerVelX,jerVelY, ...
      fldLenX,fldLenY,fldVelX,fldVelY, ...
      cenLenX,cenLenY,cenVelX,cenVelY ...
      ] = ...
      AllMatelsCore(intime,vl,vv,vlc,vvc, vlp,vvp,vlcp,vvcp,...
      VelCurX, VelCurY, ...
      LenAccX, LenAccY, VelAccX, VelAccY, ...
      LenJerkX, LenJerkY, VelJerkX, VelJerkY, ...
      LenFieldX, LenFieldY, VelFieldX, VelFieldY, ...
      VelCenX, VelCenY, ...
      cosThetaDiag, sinThetaDiag, DOPARITY)    
    
    % value
    GetExpectX = @(vc,v,~,~) vc' * diag(cosThetaDiag) * v;
    GetExpectY = @(vc,v,~,~) vc' * diag(sinThetaDiag) * v;
    if DOPARITY
      GetExpectX = @(vc,v,vcp,vp) 1/2 * ( GetExpectX(vc,v) + GetExpectX(vcp,vp) );
      GetExpectY = @(vc,v,vcp,vp) 1/2 * ( GetExpectY(vc,v) + GetExpectY(vcp,vp) );
    end
    dipVelX = GetExpectX(vvc , vv, vvcp , vvp ) ;
    dipVelY = GetExpectY(vvc , vv, vvcp , vvp ) ;
    %
    dipLenX = GetExpectX(vlc , vl, vlcp , vlp ) ;
    dipLenY = GetExpectY(vlc , vl, vlcp , vlp) ;    
    
    % derivative
    GetCurXVel = @(vc,v,~,~) sum(VelCurX(vc,intime,v),1);
    GetCurYVel = @(vc,v,~,~) sum(VelCurY(vc,intime,v),1);
    GetCurXLen = @(vc,v,~,~) sum(VelCurX(vc,-99,v),1);
    GetCurYLen = @(vc,v,~,~) sum(VelCurY(vc,-99,v),1);
    if DOPARITY
      GetCurXVel = @(vc,v,vcp,vp) 1/2 * ( GetCurXVel(vc,v) + GetCurXVel(vcp,vp) );
      GetCurYVel = @(vc,v,vcp,vp) 1/2 * ( GetCurYVel(vc,v) + GetCurYVel(vcp,vp) );
      GetCurXLen = @(vc,v,vcp,vp) 1/2 * ( GetCurXLen(vc,v) + GetCurXLen(vcp,vp) );
      GetCurYLen = @(vc,v,vcp,vp) 1/2 * ( GetCurYLen(vc,v) + GetCurYLen(vcp,vp) );
    end
    curVelX = GetCurXVel(vvc , vv,  vvcp , vvp ) ;
    curVelY = GetCurYVel(vvc , vv,  vvcp , vvp  ) ;
    %
    curLenX = GetCurXLen(vlc , vl,  vlcp , vlp ) ;
    curLenY = GetCurYLen(vlc , vl,  vlcp , vlp  ) ;
    
    if nargout > 8
      % acceleration
      GetAccXVel = @(vc,v,~,~) sum( VelAccX(vc,intime,v), 1);
      GetAccYVel = @(vc,v,~,~) sum( VelAccY(vc,intime,v), 1);
      GetAccXLen = @(vc,v,~,~) sum( LenAccX(vc,intime,v), 1);
      GetAccYLen = @(vc,v,~,~) sum( LenAccY(vc,intime,v), 1);
      if DOPARITY
        GetAccXVel = @(vc,v,vcp,vp) 1/2 * ( GetAccXVel(vc,v) + GetAccXVel(vcp,vp) );
        GetAccYVel = @(vc,v,vcp,vp) 1/2 * ( GetAccYVel(vc,v) + GetAccYVel(vcp,vp) );
        GetAccXLen = @(vc,v,vcp,vp) 1/2 * ( GetAccXLen(vc,v) + GetAccXLen(vcp,vp) );
        GetAccYLen = @(vc,v,vcp,vp) 1/2 * ( GetAccYLen(vc,v) + GetAccYLen(vcp,vp) );
      end
      accVelX = GetAccXVel(vvc , vv,  vvcp , vvp ) ;
      accVelY = GetAccYVel(vvc , vv,  vvcp , vvp  ) ;
      %
      accLenX = GetAccXLen(vlc , vl,  vlcp , vlp ) ;
      accLenY = GetAccYLen(vlc , vl,  vlcp , vlp  ) ;
    end
    
    if nargout > 12
      % jerk
      GetJerkXVel = @(vc,v,~,~) sum( VelJerkX(vc,intime,v), 1);
      GetJerkYVel = @(vc,v,~,~) sum( VelJerkY(vc,intime,v), 1);
      GetJerkXLen = @(vc,v,~,~) sum( LenJerkX(vc,intime,v), 1);
      GetJerkYLen = @(vc,v,~,~) sum( LenJerkY(vc,intime,v), 1);
      if DOPARITY
        GetJerkXVel = @(vc,v,vcp,vp) 1/2 * ( GetJerkXVel(vc,v) + GetJerkXVel(vcp,vp) );
        GetJerkYVel = @(vc,v,vcp,vp) 1/2 * ( GetJerkYVel(vc,v) + GetJerkYVel(vcp,vp) );
        GetJerkXLen = @(vc,v,vcp,vp) 1/2 * ( GetJerkXLen(vc,v) + GetJerkXLen(vcp,vp) );
        GetJerkYLen = @(vc,v,vcp,vp) 1/2 * ( GetJerkYLen(vc,v) + GetJerkYLen(vcp,vp) );
      end
      jerVelX = GetJerkXVel(vvc , vv,  vvcp , vvp ) ;
      jerVelY = GetJerkYVel(vvc , vv,  vvcp , vvp  ) ;
      %
      jerLenX = GetJerkXLen(vlc , vl,  vlcp , vlp ) ;
      jerLenY = GetJerkYLen(vlc , vl,  vlcp , vlp  ) ;
    end
    
    if nargout > 16
      % field
      GetFldXVel = @(vc,v,~,~) sum( VelFieldX(vc,intime,v), 1);
      GetFldYVel = @(vc,v,~,~) sum( VelFieldY(vc,intime,v), 1);
      GetFldXLen = @(vc,v,~,~) sum( LenFieldX(vc,intime,v), 1);
      GetFldYLen = @(vc,v,~,~) sum( LenFieldY(vc,intime,v), 1);
      if DOPARITY
        GetFldXVel = @(vc,v,vcp,vp) 1/2 * ( GetFldXVel(vc,v) + GetFldXVel(vcp,vp) );
        GetFldYVel = @(vc,v,vcp,vp) 1/2 * ( GetFldYVel(vc,v) + GetFldYVel(vcp,vp) );
        GetFldXLen = @(vc,v,vcp,vp) 1/2 * ( GetFldXLen(vc,v) + GetFldXLen(vcp,vp) );
        GetFldYLen = @(vc,v,vcp,vp) 1/2 * ( GetFldYLen(vc,v) + GetFldYLen(vcp,vp) );
      end
      fldVelX = GetFldXVel(vvc , vv,  vvcp , vvp ) ;
      fldVelY = GetFldYVel(vvc , vv,  vvcp , vvp  ) ;
      %
      fldLenX = GetFldXLen(vlc , vl,  vlcp , vlp ) ;
      fldLenY = GetFldYLen(vlc , vl,  vlcp , vlp  ) ;
    end
    
    if nargout > 20
      % centripetal
      GetCenXVel = @(vc,v,~,~) sum( VelCenX(vc,intime,v), 1);
      GetCenYVel = @(vc,v,~,~) sum( VelCenY(vc,intime,v), 1);
      GetCenXLen = @(vc,v,~,~) sum( VelCenX(vc,-99,v), 1);
      GetCenYLen = @(vc,v,~,~) sum( VelCenY(vc,-99,v), 1);
      if DOPARITY
        GetCenXVel = @(vc,v,vcp,vp) 1/2 * ( GetCenXVel(vc,v) + GetCenXVel(vcp,vp) );
        GetCenYVel = @(vc,v,vcp,vp) 1/2 * ( GetCenYVel(vc,v) + GetCenYVel(vcp,vp) );
        GetCenXLen = @(vc,v,vcp,vp) 1/2 * ( GetCenXLen(vc,v) + GetCenXLen(vcp,vp) );
        GetCenYLen = @(vc,v,vcp,vp) 1/2 * ( GetCenYLen(vc,v) + GetCenYLen(vcp,vp) );
      end
      cenVelX = GetCenXVel(vvc , vv,  vvcp , vvp ) ;
      cenVelY = GetCenYVel(vvc , vv,  vvcp , vvp  ) ;
      %
      cenLenX = GetCenXLen(vlc , vl,  vlcp , vlp ) ;
      cenLenY = GetCenYLen(vlc , vl,  vlcp , vlp  ) ;
    end
    
  end  % AllMatelsCore

%%%%   ALLMATELS WRAPPER  %%%%

AllMatels = @(intime,vl,vv,vlc,vvc,vlp,vvp,vlcp,vvcp) ...
  AllMatelsCore(intime,vl,vv,vlc,vvc,vlp,vvp,vlcp,vvcp, ...
  VelCurX, VelCurY, ...
  LenAccX, LenAccY, VelAccX, VelAccY, ...
  LenJerkX, LenJerkY, VelJerkX, VelJerkY, ...
  LenFieldX, LenFieldY, VelFieldX, VelFieldY, ...
  VelCenX, VelCenY, ...  
  cosThetaDiag, sinThetaDiag, DOPARITY);

clear AllMatelsCore VelCurX VelCurY VelCenX VelCenY cosThetaDiag sinThetaDiag
clear AccelX AccelY VelJerkX VelJerkY LenJerkX LenJerkY VelAccX VelAccY LenAccX LenAccY
clear LenFieldX LenFieldY VelFieldX VelFieldY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial values
%
[   ~,~, eStart ]                       = OvlStuff(...
  startVec,startVec,startVec,startVec,  startVecp,startVecp,startVecp,startVecp);
%
% ax,ay,jx,jy,cx,cystart should be zero certainly
[ exStart,eyStart, ~,~, ...
  dxStart,dyStart, ~,~, ...
  axStart,ayStart, ~,~, ...
  jxStart,jyStart, ~,~, ...
  fxStart,fyStart, ~,~, ...
  cxStart,cyStart   ...
  ] =   AllMatels(-99, ...
  startVec,startVec,startVec,startVec,  startVecp,startVecp,startVecp,startVecp);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

enVel = zeros(nTime,1);
enLen = zeros(nTime,1);
eOvlVel = zeros(nEig,nTime);
eOvlLen = zeros(nEig,nTime);

enVel(1) = eStart;
enLen(1) = eStart;

dipVelX = zeros(nTime,1);
dipLenX = zeros(nTime,1);
dipVelY = zeros(nTime,1);
dipLenY = zeros(nTime,1);
dipVelX(1) = exStart;
dipLenX(1) = exStart;
dipVelY(1) = eyStart;
dipLenY(1) = eyStart;

curVelX = zeros(nTime,1);
curLenX = zeros(nTime,1);
curVelY = zeros(nTime,1);
curLenY = zeros(nTime,1);
curVelX(1) = dxStart;
curLenX(1) = dxStart;
curVelY(1) = dyStart;
curLenY(1) = dyStart;

accVelX = zeros(nTime,1);
accLenX = zeros(nTime,1);
accVelY = zeros(nTime,1);
accLenY = zeros(nTime,1);
accVelX(1) = axStart;
accLenX(1) = axStart;
accVelY(1) = ayStart;
accLenY(1) = ayStart;

jerVelX = zeros(nTime,1);
jerLenX = zeros(nTime,1);
jerVelY = zeros(nTime,1);
jerLenY = zeros(nTime,1);
jerVelX(1) = jxStart;
jerLenX(1) = jxStart;
jerVelY(1) = jyStart;
jerLenY(1) = jyStart;

fldVelX = zeros(nTime,1);
fldLenX = zeros(nTime,1);
fldVelY = zeros(nTime,1);
fldLenY = zeros(nTime,1);
fldVelX(1) = fxStart;
fldLenX(1) = fxStart;
fldVelY(1) = fyStart;
fldLenY(1) = fyStart;

cenVelX = zeros(nTime,1);
cenLenX = zeros(nTime,1);
cenVelY = zeros(nTime,1);
cenLenY = zeros(nTime,1);
cenVelX(1) = cxStart;
cenLenX(1) = cxStart;
cenVelY(1) = cyStart;
cenLenY(1) = cyStart;

vv   = startVec;
vl   = vv;
vvc  = vv;
vlc  = vv;
vvp  = startVecp;
vlp  = vvp;
vvcp = vvp;
vlcp = vvp;

clear startVec startVecp dxStart dyStart exStart eyStart cxStart cyStart
clear axStart ayStart jxStart jyStart

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

% time, norm, norm, en, en, ovl, ovl   (len and vel)
NicePrint([0,0,0,eStart,eStart,1,1])

workLenX = 0;
workLenY = 0;
workVelX = 0;
workVelY = 0;

pulseDone = false;
pulseDone2 = false;
pulseDone3 = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%      PROPAGATE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DOFTan = nargout > 8 || doPlot;

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
  
  if DOFTan
    [ dipLenX(itime),dipLenY(itime),dipVelX(itime),dipVelY(itime), ...
      curLenX(itime),curLenY(itime),curVelX(itime),curVelY(itime), ...
      accLenX(itime),accLenY(itime),accVelX(itime),accVelY(itime), ...
      jerLenX(itime),jerLenY(itime),jerVelX(itime),jerVelY(itime), ...
      fldLenX(itime),fldLenY(itime),fldVelX(itime),fldVelY(itime), ...
      cenLenX(itime),cenLenY(itime),cenVelX(itime),cenVelY(itime) ] = ...
      AllMatels(endtime,vl,vv,vlc,vvc,vlp,vvp,vlcp,vvcp);
  else
    [ dipLenX(itime),dipLenY(itime),dipVelX(itime),dipVelY(itime), ...
      curLenX(itime),curLenY(itime),curVelX(itime),curVelY(itime) ] = ...
      AllMatels(endtime,vl,vv,vlc,vvc,vlp,vvp,vlcp,vvcp);
  end
  
  if ~pulseDone
    [~,~,efxtemp,efytemp] =  GetFieldsHere(endtime,false,phasePos,phaseNeg,phaseMed);
    ff = -1;
    if itime == nTime
      ff = -0.5;
    end
    % efxtemp = real(efxtemp);
    % efytemp = real(efytemp);
    workVelX = workVelX + curVelX(itime) * efxtemp * dTime * ff;
    workVelY = workVelY + curVelY(itime) * efytemp * dTime * ff;
    workLenX = workLenX + curLenX(itime) * efxtemp * dTime * ff;
    workLenY = workLenY + curLenY(itime) * efytemp * dTime * ff;
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

%%%%%%%%%%%%%%%%%%%%        DONE PROPAGATING.       %%%%%%%%%%%%%%%%%%%%%%

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
  fprintf('WORK Tot: Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(workVelTot),real(workLenTot),imag(workVelTot),imag(workLenTot));
  fprintf('WORK X:   Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(workVelX),real(workLenX),imag(workVelX),imag(workLenX));
  fprintf('WORK Y:   Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(workVelY),real(workLenY),imag(workVelY),imag(workLenY));
  % fprintf('WORK P: Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(workVelP),real(workLenP),imag(workVelP),imag(workLenP));
  % fprintf('WORK M: Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(workVelM),real(workLenM),imag(workVelM),imag(workLenM));
else
  fprintf('WORK:     Vel, Len   %13.9f  %13.9f    \n',real(workVel),real(workLen));
  disp('Work integral dt from commutator')
  fprintf('WORK Tot: Vel, Len   %13.9f  %13.9f    \n',real(workVelTot),real(workLenTot));
  fprintf('WORK X:   Vel, Len   %13.9f  %13.9f    \n',real(workVelX),real(workLenX));
  fprintf('WORK Y:   Vel, Len   %13.9f  %13.9f    \n',real(workVelY),real(workLenY));
  % fprintf('WORK P: Vel, Len   %13.9f  %13.9f    \n',real(workVelP),real(workLenP));
  % fprintf('WORK M: Vel, Len   %13.9f  %13.9f    \n',real(workVelM),real(workLenM));
end

if ~DOFTan
  return
end

PadStart = @(x) [ones(fdOrder,1)*x(1);x];
ChopPad = @(x) x(fdOrder+1:end-fdOrder);
[fdvec,sdvec,tdvec] = FDVec(fdOrder);
DoDiff = @(x) ChopPad(VecMult(PadStart(x),fdvec,fdOrder)) / dTime ;
cur2VelX = DoDiff(dipVelX);  % dip is time derivative of exp
cur2LenX = DoDiff(dipLenX);  % should be equal to der
cur2VelY = DoDiff(dipVelY);
cur2LenY = DoDiff(dipLenY);

DoDiff2 = @(x) ChopPad(VecMult(PadStart(x),sdvec,fdOrder)) / dTime^2 ;
acc3VelX = DoDiff2(dipVelX);
acc3LenX = DoDiff2(dipLenX);
acc3VelY = DoDiff2(dipVelY);
acc3LenY = DoDiff2(dipLenY);
acc2VelX = DoDiff(curVelX);
acc2LenX = DoDiff(curLenX);
acc2VelY = DoDiff(curVelY);
acc2LenY = DoDiff(curLenY);

jer2VelX = DoDiff(accVelX);
jer2LenX = DoDiff(accLenX);
jer2VelY = DoDiff(accVelY);
jer2LenY = DoDiff(accLenY);


eFieldX2 = eFieldX(1:end-fdOrder);
eFieldY2 = eFieldY(1:end-fdOrder);
aFieldX2 = aFieldX(1:end-fdOrder);
aFieldY2 = aFieldY(1:end-fdOrder);


if 1==0
  disp([acc3VelX(2:10), acc3LenX(2:10), acc2VelX(2:10), acc2LenX(2:10)])
  disp([cenVelX(2:10), cenLenX(2:10)]) 
  disp(' ')  
  disp([cur2VelY(2:10), cur2LenY(2:10), curVelY(2:10), curLenY(2:10)])
  disp(aFieldY(2:10))
end


if 1==0
  figure(123)
  plot(....
    1:nTime,-jerLenX.*curLenX,'-', ...
    1:nTime,-jerVelX.*curVelX,'.' )
  figure(124)
  plot(....
    1:nTime,-jerLenY.*curLenY,'-', ...
    1:nTime,-jerVelY.*curVelY,'.' )
  
  xxxLenX = curLenX + fldLenX;
  xxxLenY = curLenY + fldLenY;
  xxxVelX = curVelX + fldVelX;
  xxxVelY = curVelY + fldVelY;
  
  figure(125)
  plot(....
    1:nTime,-jerLenX.*xxxLenX,'-', ...
    1:nTime,-jerVelX.*xxxVelX,'.' )
  figure(126)
  plot(....
    1:nTime,-jerLenY.*xxxLenY,'-', ...
    1:nTime,-jerVelY.*xxxVelY,'.' )
  
  drawnow
end


lenX = [curLenX,accLenX,jerLenX];
lenY = [curLenY,accLenY,jerLenY];
velX = [curVelX,accVelX,jerVelX];
velY = [curVelY,accVelY,jerVelY];
description = 'From commmutator';

lenX2 = [cur2LenX,acc2LenX,jer2LenX];
lenY2 = [cur2LenY,acc2LenY,jer2LenY];
velX2 = [cur2VelX,acc2VelX,jer2VelX];
velY2 = [cur2VelY,acc2VelY,jer2VelY];
description2 = 'With dipole operator, then differentiate';

DoDiff3 = @(x) ChopPad(VecMult(PadStart(x),tdvec,fdOrder)) / dTime^3 ;
dip3VelX = DoDiff3(dipVelX);
dip3LenX = DoDiff3(dipLenX);
dip3VelY = DoDiff3(dipVelY);
dip3LenY = DoDiff3(dipLenY);

lenX3 = [cur2LenX,acc3LenX,dip3LenX];
lenY3 = [cur2LenY,acc3LenY,dip3LenY];
velX3 = [cur2VelX,acc3VelX,dip3VelX];
velY3 = [cur2VelY,acc3VelY,dip3VelY];
description3 = 'Finite difference from plain dipole';

if 1==0
  lenX(:,1) = lenX(:,1) + fldLenX;
  lenY(:,1) = lenY(:,1) + fldLenY;
  velX(:,1) = velX(:,1) + fldVelX;
  velY(:,1) = velY(:,1) + fldVelY;
end
if 1==0
  lenX(:,1) = - fldLenX;
  lenY(:,1) = - fldLenY;
  velX(:,1) = - fldVelX;
  velY(:,1) = - fldVelY;
end

AllFtStuff(dTime,eFieldX,eFieldY,aFieldX,aFieldY,...
  lenX,lenY,velX,velY,doPlot,0,description,omega);
AllFtStuff(dTime,eFieldX2,eFieldY2,aFieldX2,aFieldY2,...
  lenX2,lenY2,velX2,velY2,doPlot,1000,description2,omega);
AllFtStuff(dTime,eFieldX2,eFieldY2,aFieldX2,aFieldY2,...
  lenX3,lenY3,velX3,velY3,doPlot,2000,description3,omega);

disp(' ')
disp('OKAY done with run.');
%  hit enter')
% pause

end  % END CIRCLEPROP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    END CIRCLEPROP    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% omega just for plotting
function [eigVects,eigVals] = CircleEigs(problemOption,nPts,fdOrder,vMax,iSel,doOut,FLIPPARITY,omega)

[~,~,ham0,potential,cosThetaDiag,sinThetaDiag,thetaVals] = CircleParts(problemOption,nPts,fdOrder,vMax,FLIPPARITY);

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

maxSel = max(iSel);
eigVals0 = eigVals0(1:maxSel);

if doOut
  disp('all transition energies')
  te = eigVals0(:)-eigVals0(:).';
  for ii=1:maxSel
    fprintf(' %8.4f',te(1:ii,ii))
    fprintf('\n')
  end
  disp('selected eigVals')
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

if doOut && 1==0
  xOpMat = eigVects' * diag(cosThetaDiag) * eigVects;
  yOpMat = eigVects' * diag(sinThetaDiag) * eigVects;
  
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

end   % END CircleEigs


%%%%%%%%%%
function varargout = GetOut(iwhich,x,Fun)
iwhich             = iwhich(1:nargout);
barargout          = cell(1,max(iwhich));
[barargout{:}]     = Fun(x);
varargout          = cell(1,nargout);
[varargout{:}]     = barargout{iwhich};
end





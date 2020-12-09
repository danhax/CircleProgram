
function [aFieldX,aFieldY,eFieldX,eFieldY] = ...
  GetFields3B(calcTimes,strXS,strXC,strYS,strYC,oStr,omega,duration,phsh1,phsh2,doplot)

if 1==0
  totStr = sqrt(strXS^2+strXC^2+strYS^2+strYC^2);
  strXS = strXS / totStr;
  strXC = strXC / totStr;
  strYS = strYS / totStr;
  strYC = strYC / totStr;
end

insize = size(calcTimes);
calcTimes = calcTimes(:);

strXS = strXS * oStr;
strXC = strXC * oStr;
strYS = strYS * oStr;
strYC = strYC * oStr;
clear oStr;

doplot2 = doplot;
doplot = false;

[aField1,eField1] = GetField(calcTimes,1,omega,duration,phsh1,doplot2);
[aField2,eField2] = GetField(calcTimes,1,omega,duration,phsh2,doplot);

field1 = [aField1,eField1];
field2 = [aField2,eField2];

% wMode = 0:  regular, just divide pulse in two
%         1:  complex-valued rotating waves
%         2:  circular polarization
%

wMode = 2;

switch wMode
  case {0,1}
    field1X = (     -strXC + 1i * strXS) * field1;
    field1Y = (     -strYC + 1i * strYS) * field1;
    field2X = conj((-strXC + 1i * strXS) * field2 );
    field2Y = conj((-strYC + 1i * strYS) * field2 );
    
  case {2}
    
    % XS, XC, YS, YC
    strXS1 = strXS - strYC;     strXS2 = strYC + strXS;
    strXC1 = strXC + strYS;     strXC2 = strXC - strYS;
    strYS1 = strXC1;            strYS2 = -strXC2;
    strYC1 = -strXS1;           strYC2 = strXS2;
    
    field1X = (     -strXC1 + 1i * strXS1) * field1 / 1;
    field1Y = (     -strYC1 + 1i * strYS1) * field1 / 1;
    field2X = (     -strXC2 + 1i * strXS2) * field2 / 1;
    field2Y = (     -strYC2 + 1i * strYS2) * field2 / 1;
    
  case {3}
    
    field1X = (     -strXC + 1i * strXS) * field1     * 2 ;
    field1Y = (     -strYC + 1i * strYS) * field1     * 0 ;
    field2X = conj((-strXC + 1i * strXS) * field2 )   * 0 ;
    field2Y = conj((-strYC + 1i * strYS) * field2 )   * 2 ;
    
  otherwise
    error('not supported')
end

aFieldX = field1X(:,1) + field2X(:,1);
eFieldX = field1X(:,2) + field2X(:,2);
aFieldY = field1Y(:,1) + field2Y(:,1);
eFieldY = field1Y(:,2) + field2Y(:,2);

if wMode~=1 || phsh1 == phsh2
  aFieldX = real(aFieldX);
  eFieldX = real(eFieldX);
  aFieldY = real(aFieldY);
  eFieldY = real(eFieldY);
  %
else   %  leave as complex-valued rotating waves
  if 1==0
    phi = (1+sqrt(5))/2;
    ifac = 1/phi;
    % ifac = 1/3;
    aFieldX = real(aFieldX) + ifac * 1i * imag(aFieldX);
    aFieldY = real(aFieldY) + ifac * 1i * imag(aFieldY);
    eFieldX = real(eFieldX) + ifac * 1i * imag(eFieldX);
    eFieldY = real(eFieldY) + ifac * 1i * imag(eFieldY);
  end
end

aFieldX = aFieldX / 2;
eFieldX = eFieldX / 2;
aFieldY = aFieldY / 2;
eFieldY = eFieldY / 2;

aFieldX = reshape(aFieldX,insize);
eFieldX = reshape(eFieldX,insize);
aFieldY = reshape(aFieldY,insize);
eFieldY = reshape(eFieldY,insize);

end



function [aField,eField] = GetField(calcTimes,fieldStrength,omega,duration,phaseShift,doplot)

[~,aField,eField] = FieldFunction(calcTimes,omega,duration,phaseShift);

aField = aField * sqrt(2*pi*fieldStrength) / omega / sqrt(duration) ;
eField = eField * sqrt(2*pi*fieldStrength) / omega / sqrt(duration) ;

if doplot && numel(calcTimes)>1 && 1==1
  figure(5)
  plot(calcTimes,real(aField),'-',calcTimes,real(eField),'-','LineWidth',4)
  title('Fields')
  legend('A(t)','E(t)')
  yticks(0);
  xlabel('t (atomic units)');
  set(gca,'FontSize',24)
  drawnow
  
  dTime = calcTimes(2)-calcTimes(1);
  checkSum = sum(real(aField).^2)*dTime;
  disp(['                                    CheckSum aField :::: ' num2str(checkSum)])
end

end


function [y,yp,ypp] = FieldFunction(t,omega,duration,phaseShift)
%  differentiated once for A, twice for E.
%

%$$ N = 3;   % N=3 is minimal value ensuring E is zero at t=0.
%$$          % with N=3, d/dt E is not zero at t=0.

N = 5;

O = pi/duration;

%if any(t<0) || any(t> duration)
%  error('ack checkme')
%end

badt = t<0 | t>duration;

% NOW CENTERED
t = t - duration / 2;

if 1==0
  fun   = @sin;
  dfun  = @(x) cos(x);
  ddfun = @(x) -sin(x);
else
  fun   = @(x)    exp(1i*x);
  dfun  = @(x) 1i*exp(1i*x);
  ddfun = @(x) -1*exp(1i*x);
end

funt = fun(omega*t+phaseShift);
dfunt = omega * dfun(omega*t+phaseShift);
ddfunt = omega^2 * ddfun(omega*t+phaseShift);

env = cos(O*t).^N;
denv = -O * N * cos(O*t).^(N-1) .* sin(O*t);
ddenv =   ...
  O^2 * N * (N-1) * cos(O*t).^(N-2) .* sin(O*t).^2 ...
  - O^2 * N * cos(O*t).^(N-1) .* cos(O*t) ;

y = ...
  env .* funt ;

yp = ...
  denv .* funt + ...
  env .* dfunt ;

ypp = ...
  ddenv    .* funt + ...
  2 * denv .* dfunt + ...
  env      .* ddfunt ;

y(badt) = 0;
yp(badt) = 0;
ypp(badt) = 0;

end



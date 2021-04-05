
function [aFieldX,aFieldY,eFieldX,eFieldY,dFieldX,dFieldY] = ...
  GetFields3B(calcTimes,strXS,strXC,strYS,strYC,oStr,omega,duration,phsh1,phsh2,doplot,ENV_OPT)

ENV_PWR = 1;

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

doplot2 = doplot;
doplot = false;

[aField1,eField1,dField1] = GetField(calcTimes,1,omega,duration,phsh1,doplot2,ENV_OPT,ENV_PWR);
if phsh1 == phsh2
  aField2 = aField1;
  eField2 = eField1;
  dField2 = dField1;
else
  [aField2,eField2,dField2] = GetField(calcTimes,1,omega,duration,phsh2,doplot,ENV_OPT,ENV_PWR);
end

field1 = [aField1,eField1,dField1];
field2 = [aField2,eField2,dField2];

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
dFieldX = field1X(:,3) + field2X(:,3);

aFieldY = field1Y(:,1) + field2Y(:,1);
eFieldY = field1Y(:,2) + field2Y(:,2);
dFieldY = field1Y(:,3) + field2Y(:,3);

if wMode~=1 || phsh1 == phsh2
  aFieldX = real(aFieldX);
  eFieldX = real(eFieldX);
  dFieldX = real(dFieldX);
  aFieldY = real(aFieldY);
  eFieldY = real(eFieldY);
  dFieldY = real(dFieldY);
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
    
    dFieldX = real(dFieldX) + ifac * 1i * imag(dFieldX);
    dFieldY = real(dFieldY) + ifac * 1i * imag(dFieldY);
  end
end

aFieldX = aFieldX / 2;
eFieldX = eFieldX / 2;
dFieldX = dFieldX / 2;

aFieldY = aFieldY / 2;
eFieldY = eFieldY / 2;
dFieldY = dFieldY / 2;

aFieldX = reshape(aFieldX,insize);
eFieldX = reshape(eFieldX,insize);
dFieldX = reshape(dFieldX,insize);

aFieldY = reshape(aFieldY,insize);
eFieldY = reshape(eFieldY,insize);
dFieldY = reshape(dFieldY,insize);

end



function [aField,eField,dField] = GetField(calcTimes,fieldStrength,omega,duration,phaseShift,doplot,ENV_OPT,ENV_PWR)

%{
persistent cT_ fS_ o_ d_ p_ E_
persistent aF_ eF_ dF_
if isequal({calcTimes,fieldStrength,omega,duration,phaseShift,ENV_OPT}, {cT_, fS_, o_, d_, p_, E_})
  aField = aF_;
  eField = eF_;
  dField = dF_;
  return;
end
cT_ = calcTimes;
fS_ = fieldStrength;
o_ = omega;
d_ = duration;
p_ = phaseShift;
E_ = ENV_OPT;
%}

if ENV_PWR == 0
  [~,aField,eField,dField]   = FieldFunction(calcTimes,omega,duration,phaseShift,ENV_OPT);
elseif ENV_PWR == 1
  [~,~,aField,eField,dField] = FieldFunction(calcTimes,omega,duration,phaseShift-pi/2,ENV_OPT);
  aField = aField ./ omega ;
  eField = eField ./ omega ;
  dField = dField ./ omega ;
else
  error('not supported')
end

fac = sqrt(2*pi*fieldStrength) / omega / sqrt(duration);

aField = aField * fac ;
eField = eField * fac ;
dField = dField * fac ;

%{
aF_ = aField;
eF_ = eField;
dF_ = dField;
%}

if doplot && numel(calcTimes)>1 && 1==1
  figure(5-ENV_PWR)
  plot(calcTimes,real(aField),'-',calcTimes,real(eField),'-','LineWidth',4)
  title('Fields')
  legend('A(t)','E(t)')
  % yticks(0);
  xlabel('t (atomic units)');
  set(gca,'FontSize',24)
  drawnow
  
  dTime = calcTimes(2)-calcTimes(1);
  checkSum = sum(real(aField).^2)*dTime;
  disp(['                                    CheckSum aField :::: ' num2str(checkSum)])
end

end


function [y,yp,ypp,yp3,yp4] = FieldFunction(t,omega,duration,phaseShift,ENV_OPT)
%  differentiated once for A, twice for E.
%
TEMPPLOT = 1==0;

O = pi/duration;

%if any(t<0) || any(t> duration)
%  error('ack checkme')
%end

badt = t<=0 | t>=duration;

% NOW CENTERED
t = t - duration / 2;

if 1==0
  fun   = @sin;
  dfun  = @(x) cos(x);
  ddfun = @(x) -sin(x);
  d3fun = @(x) -cos(x);
  d4fun = @(x) sin(x);
else
  fun   = @(x)    exp(1i*x);
  dfun  = @(x) 1i*exp(1i*x);
  ddfun = @(x) -1*exp(1i*x);
  d3fun = @(x)-1i*exp(1i*x);
  d4fun = @(x)    exp(1i*x);
end

funt = fun(omega*t+phaseShift);
dfunt = omega * dfun(omega*t+phaseShift);
ddfunt = omega^2 * ddfun(omega*t+phaseShift);
d3funt = omega^3 * d3fun(omega*t+phaseShift);
d4funt = omega^4 * d4fun(omega*t+phaseShift);

% disp('TEMP!')
% ENV_OPT = 0;

if ENV_OPT == 0
  
  %$$ N = 3;   % N=3 is minimal value ensuring E is zero at t=0.
  %$$          % with N=3, d/dt E is not zero at t=0.

  N = 5;

  env = cos(O*t).^N;
  denv = -O * N * cos(O*t).^(N-1) .* sin(O*t);
  ddenv =   ...
    O^2 * N * (N-1) * cos(O*t).^(N-2) .* sin(O*t).^2 ...
    - O^2 * N * cos(O*t).^N ;
  %$- O^2 * N * cos(O*t).^(N-1) .* cos(O*t) ;
  
  d3env =   ...
    - O^3 * N * (N-1) * (N-2) * cos(O*t).^(N-3) .* sin(O*t).^3 ...
    + O^3 * N * (3*N-2)       * cos(O*t).^(N-1) .* sin(O*t) ;
  %$+ O^3 * N * (N-1) * 2     * cos(O*t).^(N-1) .* sin(O*t) ...
  %$+ O^3 * N^2               * cos(O*t).^(N-1) .* sin(O*t) ;

  d4env =   ...
    + O^4 * N * (N-1) * (N-2) * (N-3) * cos(O*t).^(N-4) .* sin(O*t).^4 ...
    - O^4 * N * (N-1) * (6*N-8)   * cos(O*t).^(N-2) .* sin(O*t).^2 ;      
    ...%- O^4 * N * (N-1) * (N-2) * 3 * cos(O*t).^(N-2) .* sin(O*t).^2 ...
    ...%- O^4 * N * (N-1) * (3*N-2)   * cos(O*t).^(N-2) .* sin(O*t).^2 ;
    + O^4 * N * (3*N-2)           * cos(O*t).^N  ;

  ppp = 0;
  
else
  
  if 1==0
    N = 9;
    
    env0     =        cos(O*t);
    denv0    = -O*    sin(O*t);
    ddenv0   = -O^2 * cos(O*t);
    d3env0   = O^3*   sin(O*t);
    d4env0   = O^4*   cos(O*t);
  else
    N = 12;
    
    env0     = 1-(O*2/pi)^2*t.^2;
    denv0    = -2*(O*2/pi)^2*t;
    ddenv0   = -2*(O*2/pi)^2*(t*0+1);
    d3env0   = t*0;
    d4env0   = t*0;
  end
  
  env1    = exp(N) * exp(-N./env0 + log(env0))           ;
  denv1   = (N./env0.^2 + 1./env0) .* env1               ;
  ddenv1  = (N^2./env0.^4) .* env1                       ;
  d3env1  = (-3 * N^2./env0.^5 + N^3./env0.^6 ) .* env1  ;

  d4env1  = ((15-3) .* N^2./env0.^6   +(-2-6) .* N^3./env0.^7    + N^4./env0.^8) .* env1 ;
  
  
  env     = env1;
  denv    = denv1 .* denv0;
  ddenv   = ddenv1 .* denv0.^2  +  denv1 .* ddenv0;
  d3env   = d3env1 .* denv0.^3  +  3 * ddenv1 .* denv0 .* ddenv0  +  denv1 .* d3env0;
  
  d4env   = d4env1 .* denv0.^4  + 3 * d3env1 .* denv0.^2 .* ddenv0 ...
    +       3 * (d3env1 .* denv0.^2 .* ddenv0 + ddenv1 .* ddenv0.^2 + ddenv1 .* denv0 .* d3env0) ...
    +       ddenv1 .* denv0 .* d3env0 + denv1 .* d4env0;

  ppp = 1000;
  
end


if TEMPPLOT
  env(badt) = 0;
  denv(badt) = 0;
  ddenv(badt) = 0;
  d3env(badt) = 0;
  d4env(badt) = 0;
  figure(ppp+997); plot(env)
  figure(ppp+996); plot(denv)
  figure(ppp+995); plot(ddenv)
  figure(ppp+994); plot(d3env)
  figure(ppp+894); plot(diff(ddenv))
  figure(ppp+794); plot(diff(diff(denv)))
  figure(ppp+993); plot(d4env)
  figure(ppp+893); plot(diff(d3env))
  figure(ppp+793); plot(diff(diff(ddenv)))
  figure(ppp+792); plot(d4env(2:end-1)./diff(diff(ddenv)))
  figure(ppp+791); plot(ddenv(2:end-1)./diff(diff(env)))
  disp('plotted')
end

y = ...
  env .* funt ;

yp = ...
  denv .* funt + ...
  env .* dfunt ;

ypp = ...
  ddenv    .* funt + ...
  2 * denv .* dfunt + ...
  env      .* ddfunt ;

yp3 = ...
  d3env     .* funt + ...
  3 * ddenv .* dfunt + ...
  3 * denv  .* ddfunt + ...
  env       .* d3funt ;

yp4 = ...
  d4env     .* funt + ...
  4 * d3env .* dfunt + ...
  6 * ddenv .* ddfunt + ...
  4 * denv  .* d3funt + ...
  env       .* d4funt ;

y(badt) = 0;
yp(badt) = 0;
ypp(badt) = 0;
yp3(badt) = 0;
yp4(badt) = 0;

end



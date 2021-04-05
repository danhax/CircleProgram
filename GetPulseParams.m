
function [vMax,omega,duration,strXS,strXC,strYS,strYC,ENV_OPT] = ...
  GetPulseParams(problemOption)

ENV_OPT = 1;  % 0 = orig  1 = new

vMax = Inf;

if problemOption == 0
  
  duration = {12,17.5};
  omega = 7.65;
  strXS = 0;
  strXC = -0.55;
  strYS = 0.60;
  strYC = 0.095;
  
elseif problemOption == 1
  
  duration = {1.2,1.75}; %1.7};
  omega = 76.5;
  strXS = 0.01;
  strXC = 0.005;
  strYS = 1;
  strYC = 0.5;
    
  
elseif problemOption == 2
  
  duration = {200,290};
  omega = 0.314; % 7th harmonic 1st excited state
  omega = 0.275;  % 8th
  omega = 0.365;  % 6th
  
  strXS = 0.1;
  strXC = 0.05;
  strYS = 10;
  strYC = 5;
    
  strXS = 0.2;
  strXC = 0.1;
  strYS = 20;
  strYC = 10;
  
  strXS = 0.32;
  strXC = 0.16;
  strYS = 32;
  strYC = 16;
  
  strXS = 0.4;
  strXC = 0.2;
  strYS = 40;
  strYC = 20;

  % % increase population state 1
  % strXS = 13;
  % strXC = -1;
  % % increase pop state 2
  strXS = 25;
  strXC = 12.5;
  
else
  error('not supported')
end

duration = duration{ENV_OPT+1};

end



function [vMax,omega,duration,strXS,strXC,strYS,strYC] = ...
  GetPulseParams(problemOption)

if problemOption == 0
  
%%%%%%%%%%%%%%%%

% optimized first with sqrt(pop)/work
%
% logPOPS:   (??)
%    -0.4071   -0.4086
%    -9.8779   -9.8904
%    -1.1780   -1.1744
%    -8.6071   -8.6051
%    -7.9603   -7.9559
%    -5.8168   -5.8298
%    -4.4850   -4.4843
%   -11.2384  -11.2338
%    -5.8605   -5.8561
%    -9.3926   -9.3872
%
vMax     =  23.7811 ;
omega    =   7.5180 ;
duration =  33.0295 ;
strXS   =  0.1967 ;
strXC   =  3.7981 ;
strYS   =  1.9374 ;
strYC   = -0.9028 ;

% obtained from that with pop^0.35/work
% logPOPS:
%    -0.0237   -0.0240
%   -11.0017  -11.0050
%    -4.1352   -4.1202
%   -12.0418  -12.0393
%   -10.9487  -10.9425
%    -6.5633   -6.5622
%    -6.4843   -6.4733
%   -13.1623  -13.1733
%   -10.1610  -10.1579
%    -5.8892   -5.8918
vMax     =     25.4330 ;
omega    =     7.4860;
duration =   47.0734;
strXS    =     0.1583;
strXC    =     3.7822;
strYS    =     2.1451;
strYC    =    -0.7892;

% from that, pop^0.47/work.. not much change but maybe bad run
vMax     =     25.4113;
omega    =     7.4867;
duration =    47.3585;
strXS    =     0.1588;
strXC    =     3.7809;
strYS    =     2.1503;
strYC    =    -0.7867;

% now pop^0.67/work
duration = 46.7295;

%{
strXC = 1;
strXS = 0;
strYC = 0;
strYS = 0;
%omega = 5;

%omega = 6.28;
%omega = 7;
%duration = 12;
%}
%duration = 20;

duration = 12;
FFF = 0.05;
strXS    =     strXS * FFF;
strXC    =     strXC * FFF;
strYS    =     strYS * FFF;
strYC    =     strYC * FFF;

vMax     =     Inf;
omega    =     7.8467;
duration =    12.0000;
strXS    =     0.0079;
strXC    =     0.7336;
strYS    =     0.7471;
strYC    =     0.1870;

omega    =     7.8;
duration =    12.0000;
strXS    =     0;
strXC    =     0.75;
strYS    =     0.75;
strYC    =     0.0;

%%%%%%%%%% thisone
omega    =     7.8441;
duration =    12.0000;
strXS    =     0.1172;
strXC    =     0.7235;
strYS    =     0.7663;
strYC    =     0.0734;
%%%%%%%%%%

% new opt with N=3 in Fields3.m, oStr = 1e-2
% omega    =  7.5229;
% now N=5 in Fields3.m, oStr=1e-2
% omega    =  7.5062;  
% oStr = 1e-1
% omega    =   7.5202;
% oStr = 1e0
% omega    =   7.6645;

%now strXS etc are proportions.  oStr=1e0
strXS = 0.0647;
strXC = 0.8231;
strYS = 0.7755;
strYC = 0.1200;

% oStr = 1e0 re-opt omega
omega =  7.6715;

strXS =    0.0568;
strXC =    0.7226;
strYS =    0.6808;
strYC =    0.1053;

% turn off DOAVERAGE 

strXS = 0.990448; 
strXC = 3.280258;
strYS = 1.679278;
strYC = 0.344517;

% optimize pop/work (roughly, stopped calc) :

strXS = -0.462791 ;
strXC = 3.046014 ;
strYS = 2.238873 ;
strYC = 0.608654 ;
            
% pop/work^2 roughly

strXS = -0.442914;
strXC = 2.684227;
strYS = 2.368409;
strYC = 0.763688;

% log(work)*0.7 - log(pop) 
omega = 7.5588;
strXS = 0;
strXC = 1;
strYS = 0;
strYC = 0;

strXS = 0;
strXC = 0.529839;
strYS = 0.579521;
strYC = 0.080951;


omega = 7.642798;

strXS = 0;
strXC = 0.550109;
strYS = 0.597652;
strYC = 0.093775;

omega = 7.646306;

%$$$ GOOD FIT log(work)*0.7 - log(pop) 
omega =  7.646519;
strXS = 0;
strXC = 0.550740;
strYS = 0.597872;
strYC = 0.094484;

%$$$ GOOD FIT log(work)*0.7 - log(pop) 
omega =  7.65;
strXS = 0;
strXC = 0.55;
strYS = 0.60;
strYC = 0.095;

strXS = 0;
strXC = 0.67;
strYS = 0.70;
strYC = 0.136;

strXS = 0;
strXC = 1.14;
strYS = 0.75;
strYC = 0.32;

omega = 7.38;
strXS = 0;
strXC = 1.61;
strYS = 2.03;
strYC = 0.72;

duration = 15;
strXS = 0;
strXC = 2;
strYS = 2.25;
strYC = 0.8;

%strXS = 2;
%strXC = 0;
%strYS = 0.8;
%strYC = 2.25;

duration = 12;
omega = 7.65;
strXS = 0;
strXC = 0.55;
strYS = 0.60;
strYC = 0.095;


else
  error('not supported')
end

end


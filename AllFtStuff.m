
function AllFtStuff(dTime,eFieldX,eFieldY,aFieldX,aFieldY,...
  lenX,lenY,velX,velY,doPlot,plotNum,description,omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%      FT      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% efield and TIME DERIVATIVEs of dipole moments input
[ftVelX,ftVelY,ftLenX,ftLenY] = FtDip(dTime,velX,velY,lenX,lenY);

% ftVelX, etc cell arrays size (nDers,2)   2 = pos/neg freq

[cumVelT,cumVelX,cumVelY,cumVelP,cumVelM, ...
  cumLenT,cumLenX,cumLenY,cumLenP,cumLenM, hhgVelT, hhgLenT  ] = ...
  DoFtStuffCore(...
  eFieldX, eFieldY, aFieldX, aFieldY, ...
  ftVelX, ftVelY, ftLenX, ftLenY, ...
  dTime,omega,doPlot,doPlot,plotNum);

disp(['Work integral domega: ' description])

if any(abs(imag([hhgVelT,hhgLenT])) > 1e-8)
  fprintf('HHG  Tot: Vel, Len   %13.9f  %13.9f    %13.9f  %13.9f \n',real(hhgVelT),real(hhgLenT),imag(hhgVelT),imag(hhgLenT));
else
  fprintf('HHG  Tot: Vel, Len   %13.9f  %13.9f \n',real(hhgVelT),real(hhgLenT));
end

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


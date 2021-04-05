%#ok<*AGROW>

function [cumVelTout,cumVelXout,cumVelYout,cumVelPout,cumVelMout, ...
  cumLenTout,cumLenXout,cumLenYout,cumLenPout,cumLenMout ] = ...
  FourierAseCore(...
  eFieldX, eFieldY, ~, ~, ...   %   eFieldX, eFieldY, aFieldX, aFieldY, ...
  ftVelX, ftVelY, ftLenX, ftLenY, ...
  dTime,omega,doPlot,plotNum)

[eFTX{1},eFTX{2},om] = TimeFT(eFieldX,dTime,0) ;
[eFTY{1},eFTY{2}]    = TimeFT(eFieldY,dTime,0) ;

for ii = 1:2
  %$$
  ff = (-1)^(ii-1);
  % ff = 1;
  eFTP{ii}   = 1/sqrt(2) * (eFTX{ii}   + ff * 1i * eFTY{ii});  
  ftLenP{ii} = 1/sqrt(2) * (ftLenX{ii} + ff * 1i * ftLenY{ii});
  ftVelP{ii} = 1/sqrt(2) * (ftVelX{ii} + ff * 1i * ftVelY{ii});
  
  eFTM{ii}   = 1/sqrt(2) * (eFTX{ii}   - ff * 1i * eFTY{ii});
  ftLenM{ii} = 1/sqrt(2) * (ftLenX{ii} - ff * 1i * ftLenY{ii});
  ftVelM{ii} = 1/sqrt(2) * (ftVelX{ii} - ff * 1i * ftVelY{ii});
end


% PFun = @(x1,x2,y1,y2) ( (x1+x2) + 1i * (y1-y2) ) / 2 ;
% MFun = @(x1,x2,y1,y2) ( (y1+y2) + 1i * (x1-x2) ) / 2 ;

%{
  function [pf,mf] = PMF(x1,x2,y1,y2)
    pf = ( (x1+x2) + 1i * (y1-y2) ) / 2 ;
    mf = ( (y1+y2) + 1i * (x1-x2) ) / 2 ;
  end
for ii = 1:2
  jj = 3-ii;
  [   eFTP{ii},    eFTM{ii}] = PMF(  eFTX{ii},   eFTX{jj},   eFTY{ii},   eFTY{jj});
  [ ftLenP{ii},  ftLenM{ii}] = PMF(ftLenX{ii}, ftLenX{jj}, ftLenY{ii}, ftLenY{jj});
  [ ftVelP{ii},  ftVelM{ii}] = PMF(ftVelX{ii}, ftVelX{jj}, ftVelY{ii}, ftVelY{jj});
end
%}


specFac = 1 / pi * dTime^2;
dom = om(2)-om(1);

% % yes, this one..  sum over positive and negative frequencies;
%  for either positive or negative, take real part of dot product
% SpecFun = @(eft,dipft) ...
%   -1/2 * real(conj(eft{1}) .* dipft{1}) .* specFac ...
%   -1/2 * real(conj(eft{2}) .* dipft{2}) .* specFac;

% % no, not this one
% SpecFun = @(eft,dipft) -1/2 * specFac * real( ...

% iDer:  1st,2nd,3rd ddt derivatives of dipole moment in slots 1,2,3
% check iDer setting in TimeFT too
%
%$$ iDer = 3;    % for iDer=1 in TimeFT
% iDer = 2;      % for iDer=2 in TimeFT
iDer = 1;        % for iDer=3 in TimeFT

% differentiating causes problem at zero in FT for low frequency.
% now setting iDer = 1 here and iDer = 0 in TimeFT1.m

% need iDer = 2 for smooth result here for AnE at high frequency
iDer = 2; 

% but iDer = 1 gives smoother AnE at low frequency
iDer = 1;

  function x = ZeroZero(om, x)
    if iDer > 1
      x(om==0,:) = 0;
    end
  end

SpecFun = @(eft,dipft) -1/2 * specFac * ZeroZero(om, om.^(1-iDer) .* mean(...
  (+1i).^(1-iDer) .* eft{2}(:,:) .* squeeze(dipft{1}(:,iDer,:)) + ...
  (-1i).^(1-iDer) .* eft{1}(:,:) .* squeeze(dipft{2}(:,iDer,:)), 2)) ...
  ;

% SpecFun = @(eft,dipft) ...
%   -1 * real(eft{2} .* dipft{1}) .* specFac ;

% SpecFun = @(eft,dipft) ...
%  -1 * mean(conj(eft{1}) .* dipft{1},2) .* specFac;

specVelX = SpecFun(eFTX, ftVelX);
specLenX = SpecFun(eFTX, ftLenX);
specVelY = SpecFun(eFTY, ftVelY);
specLenY = SpecFun(eFTY, ftLenY);
specVelP = SpecFun(eFTP, ftVelP);
specLenP = SpecFun(eFTP, ftLenP);
specVelM = SpecFun(eFTM, ftVelM);
specLenM = SpecFun(eFTM, ftLenM);

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
      yrangeP = 1.5*max(abs(specVelP)/yunit);
      yrangeM = 1.5*max(abs(specVelM)/yunit);
      yrangeX = 1.5*max(abs(specVelX)/yunit);
      yrangeY = 1.5*max(abs(specVelY)/yunit);
      yrange = max(1e-16,min([yrangeP,yrangeM,yrangeX,yrangeY]));
      ylim([-yrange yrange]);
    end
  end
  
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
    yticks([0,totEnd,10*totEnd,100*totEnd,1000*totEnd,10000*totEnd,100000*totEnd,1000000*totEnd]);
    yticklabels({'0',...
      sprintf('%5.3e',totEnd),sprintf('%5.3e',10*totEnd),...
      sprintf('%5.3e',100*totEnd),sprintf('%5.3e',1000*totEnd),sprintf('%5.3e',10000*totEnd),...
      sprintf('%5.3e',100000*totEnd),sprintf('%5.3e',1000000*totEnd)...
      })
    % ytickformat('%4.2e')
    %ylim(totEnd*[-0.1,1.2])
    ylim(totEnd*[-0.1,Inf])
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


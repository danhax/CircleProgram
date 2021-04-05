%#ok<*AGROW>

function [cumVelTout,cumVelXout,cumVelYout,cumVelPout,cumVelMout, ...
  cumLenTout,cumLenXout,cumLenYout,cumLenPout,cumLenMout ] = ...
  FourierHhgCore(...
  eFieldX, ~, ~, ~, ...   %   eFieldX, eFieldY, aFieldX, aFieldY, ...
  ftVelX, ftVelY, ftLenX, ftLenY, ...
  dTime,omega,doPlot,plotNum)

[~,~,om] = TimeFT(eFieldX,dTime,0) ;

% [eFTX{1},eFTX{2},om] = TimeFT(eFieldX,dTime,0) ;
% [eFTY{1},eFTY{2}]    = TimeFT(eFieldY,dTime,0) ;

for ii = 1:2
  %$$
  ff = (-1)^(ii-1);
  % ff = 1;
  % eFTP{ii}   = 1/sqrt(2) * (eFTX{ii}   + ff * 1i * eFTY{ii});  
  ftLenP{ii} = 1/sqrt(2) * (ftLenX{ii} + ff * 1i * ftLenY{ii});
  ftVelP{ii} = 1/sqrt(2) * (ftVelX{ii} + ff * 1i * ftVelY{ii});
  
  % eFTM{ii}   = 1/sqrt(2) * (eFTX{ii}   - ff * 1i * eFTY{ii});
  ftLenM{ii} = 1/sqrt(2) * (ftLenX{ii} - ff * 1i * ftLenY{ii});
  ftVelM{ii} = 1/sqrt(2) * (ftVelX{ii} - ff * 1i * ftVelY{ii});
end

alpha = 1/137;

specFac = 4 * pi / 3 * dTime^2 * alpha^3 ;
dom = om(2)-om(1);

%$$  EmitFun = @(x) 1/2 * specFac * mean( abs(x{1}.^2) + abs(x{2}.^2), 2);

% iDer:  1st,2nd,3rd ddt derivatives of dipole moment in slots 1,2,3

% iDer = 1;  
% EmitFun = @(x) specFac * mean( squeeze(x{1}(:,iDer,:).*x{2}(:,iDer,:)), 2);

% iDer = 1;   jDer = 1;  
% iDer = 2;   jDer = 2; 
% iDer = 3;   jDer = 3; 
%
% iDer = 1;   jDer = 2; 
iDer = 1;   jDer = 3; 

% zz      = ones(size(om));
% omx     = om;
% omx(1)  = 1;
% zz(1)   = 0;

  function x = ZeroZero(iDer,jDer,om, x)
    if (iDer+jDer > 4)
      x(om==0,:) = 0;
    end
  end

EmitFun = @(x) specFac * 1/2 * ZeroZero(iDer,jDer,om, om.^(-(iDer+jDer-4)) .* ( ...
  mean( (-1i).^(iDer-jDer) .* squeeze(x{1}(:,iDer,:).*x{2}(:,jDer,:)), 2) + ...
  mean( ( 1i).^(iDer-jDer) .* squeeze(x{1}(:,jDer,:).*x{2}(:,iDer,:)), 2) ));


%EmitFun = @(x) specFac * fq * 1/2 * zz .* omx.^(-(iDer+jDer-2)) .* ( ...
%  mean( squeeze(x{1}(:,iDer,:).*x{2}(:,jDer,:)), 2) + ...
%  mean( squeeze(x{1}(:,jDer,:).*x{2}(:,iDer,:)), 2) );


emitVelX = EmitFun(ftVelX);
emitLenX = EmitFun(ftLenX);
emitVelY = EmitFun(ftVelY);
emitLenY = EmitFun(ftLenY);
emitVelP = EmitFun(ftVelP);
emitLenP = EmitFun(ftLenP);
emitVelM = EmitFun(ftVelM);
emitLenM = EmitFun(ftLenM);

CumEmitFun = @(x) cumsum(x) * dom ;

cumVelX = CumEmitFun(emitVelX);
cumVelY = CumEmitFun(emitVelY);
cumVelP = CumEmitFun(emitVelP);
cumVelM = CumEmitFun(emitVelM);
cumLenX = CumEmitFun(emitLenX);
cumLenY = CumEmitFun(emitLenY);
cumLenP = CumEmitFun(emitLenP);
cumLenM = CumEmitFun(emitLenM);

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

Posreal = @(x) max(1e-60,real(x));

if doPlot
  
  figure(499+plotNum)
  if 1==0
    loglog(...
      ...%       om, real(emitVelX), 'o', om, real(emitLenX),'-', ...
      ...%       om, real(emitVelY), 'o', om, real(emitLenY),'-', ...
      ...%       om, real(emitVelP), 'x', om, real(emitLenP),'--', ...
      ...%       om, real(emitVelM), 'x', om, real(emitLenM),'--', ...
      om, Posreal(emitLenX),'-', ...
      om, Posreal(emitLenY),'-', ...
      om, Posreal(emitLenP),'--', ...
      om, Posreal(emitLenM),'--', ...
      om, Posreal(emitVelX), 'o', ...
      om, Posreal(emitVelY), 'o', ...
      om, Posreal(emitVelP), 'x', ...
      om, Posreal(emitVelM), 'x', ...
      plotArgs{:} )
    set(gca,'FontSize',FS)
    xlim([0.7 10]*omega)
    xlabel('\omega (Hartree)')
  else
    % lom = log(om+omega)/log(omega) - 1;
    % lom = (log(om/omega + 1/2) - log(1/2)) / log(2);
    lom = om/omega;
    lom(lom>1) = log(lom(lom>1))/log(2)+1;
    semilogy(...
      lom, Posreal(emitLenX),'-', ...
      lom, Posreal(emitLenY),'-', ...
      ...%lom, Posreal(emitLenP),'--', ...
      ...%lom, Posreal(emitLenM),'--', ...
      lom, Posreal(emitVelX), 'o', ...
      lom, Posreal(emitVelY), 'o', ...
      ...%lom, Posreal(emitVelP), 'x', ...
      ...%lom, Posreal(emitVelM), 'x', ...
      plotArgs{:} )
    set(gca,'FontSize',FS)
    xlim([0 Inf])
    xlabel('x=log_2(\omega/\omega_0)+1 (x>1) or x=\omega/\omega_0 (x<1)')
    xticks([0,1,2,3,4])
    xticklabels({'H0','H1','H2','H4','H8'})
  end
  % ylim([10^-30 1])
  ylim([10^-20 10^-2])
  
  %   if ~doPlot  % it is meant that we are plotting emission.. temporary hack for logic
  %     title('|d^2/dt^2 D|^2')
  %   else
  %     title('|d/dt D|^2')
  %   end
  title('Thomson, FID, and HHG')
  
  % legend('Vx','Lx','Vy','Ly','V+','L+','V-','L-')
  % legend('Lx','Ly','L+','L-','Vx','Vy','V+','V-')
  legend('Lx','Ly','Vx','Vy')
  %$$ % ylim([10^-10 Inf])
  %$$ ylim([10^-14 10])
  %$$ ylim([10^-19 0.01])
  drawnow
end

if doPlot
  
  figure(498+plotNum)
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
  semilogx(mystuff{:},  plotArgs{:});
  set(gca,'FontSize',FS)
  if totEnd > 0
    yticks([0,totEnd]);
    yticklabels({'0',sprintf('%5.3e',totEnd)})
    % ytickformat('%4.2e')
    ylim(totEnd*[-0.1,1.2])
  end
  
  legend(mylegend{:},'Location','SouthEast');
  title('Total HHG integral HHG(\omega) d\omega')
  xlim([0.7 10]*omega)
  ylabel('Energy (Hartree)')
  xlabel('Frequency \omega (inverse Hartree)')
  
  vec_pos = get(get(gca, 'YLabel'), 'Position');
  %set(get(gca, 'YLabel'), 'Position', vec_pos + [0.7 0 0]);
  set(get(gca, 'YLabel'), 'Position', vec_pos + [0.1 0 0]);
  
  drawnow
end


end


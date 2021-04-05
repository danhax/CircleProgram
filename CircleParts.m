

function [fdmat,sdmat,ham0,potential,cosThetaDiag,sinThetaDiag,thetaVals,...
  gradpot,hesspot] = ...
  CircleParts(problemOption,nPts,fdOrder,vMax,FLIPPARITY)

dTheta = 2*pi/nPts;
thetaVals  = ( 0.5:nPts-0.5 ).' * dTheta;

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

[potential,gradpot,hesspot] = GetPotential(problemOption,thetaVals,vMax,FLIPPARITY);

kemat = -1/2 * sdmat ;

ham0 = kemat + diag(potential);

% cosThetaOp = diag(cos(thetaVals));
% sinThetaOp = diag(sin(thetaVals));

cosThetaDiag = cos(thetaVals);
sinThetaDiag = sin(thetaVals);

end




function [v,gv,hv] = GetPotential(problemOption,x,~,FLIPPARITY)

%$$      v = (-sin(2*x)/4 - cos(x)/5) * 90 / pi ;

%$$ if any(problemOption == [0,2])
if any(problemOption == 0)
  if FLIPPARITY
    v =  (-sin(2*x)/4 + cos(x)/5) * 90 / pi ;
    gv = (-cos(2*x)/2 - sin(x)/5) * 90 / pi ;
    hv = ( sin(2*x)   - cos(x)/5) * 90 / pi ;
  else
    v =  (-sin(2*x)/4 - cos(x)/5) * 90 / pi ;
    gv = (-cos(2*x)/2 + sin(x)/5) * 90 / pi ;
    hv = ( sin(2*x)   + cos(x)/5) * 90 / pi ;
  end
elseif problemOption == -44 %1
  if FLIPPARITY
    error('what? FLIPPARITY not necessary')
  end
  v  = 0 * x;
  gv = 0 * x;
  hv = 0 * x;
elseif any(problemOption == [1,2])
  %PVS = 0.0123456;
  PVS = 0;
  x   = x + PVS;
  if FLIPPARITY
    v  = cos( x)/5 * 90 / pi ;
    gv = -sin(x)/5 * 90 / pi ;
    hv = -cos(x)/5 * 90 / pi ;
  else
    v  = -cos(x)/5 * 90 / pi ;
    gv = sin( x)/5 * 90 / pi ;
    hv = cos( x)/5 * 90 / pi ;
  end
else
  error('not supported')
end

end



% function x = SmoothMax(x,vMax,ss)
% x = FlatFun((x-vMax)/ss) * ss + vMax;
% end
% 
% 
% function x = FlatFun(x)
% x = 0.5 * ( x - sqrt(x.^2+1) );
% end




function [xftpos,xftneg,omega] = TimeFT1(xin,dTime,iDer)  % dTime just for omega

fdOrder = 7;

if ~exist('iDer','var')
  % 
  iDer = 0;
  %$$ iDer = 1;
  % iDer = 2;    % small error in HHG integral, but can reduce with iDer=1
  % iDer = 3;    % unacceptable error in HHG integral
  %              %     that's with iDer=1, jDer=3 in HHG integral
end

% nx = size(xin,2);

% PadStart = @(x) [ones(fdOrder,1)*x(1,:);x];

PadStart = @(x) [repmat(x(1, :,:,:,:, :,:,:,:),fdOrder,1);x];

ChopPad = @(x) x(fdOrder+1:end-fdOrder, :,:,:,:, :,:,:,:);

[fdvec,sdvec,tdvec] = FDVec(fdOrder);
dervecs = {fdvec,sdvec,tdvec};

if iDer == 0
  DoDiff = @(x) x(1:end-fdOrder, :,:,:,:, :,:,:,:);
else
  DoDiff = @(x) ChopPad(VecMult(PadStart(x),dervecs{iDer},fdOrder)) / dTime^iDer ;
end

xdiff = DoDiff(xin);

[xftpos,xftneg,omega] = TimeFT0(xdiff,dTime);

if iDer~=0
  xftpos = xftpos ./ ( 1i * omega).^iDer;
  xftneg = xftneg ./ (-1i * omega).^iDer;
  xftpos(omega==0,:) = 0;
  xftneg(omega==0,:) = 0;
end

end


function [xftpos,xftneg,omega] = TimeFT0(xin,dTime)  % dTime just for omega

xsize = size(xin);
% nin = xsize(1);

% % filter = (nin-1 : -1 : 0)' / (nin-1);
% filter = (nin : -1 : 1)' / nin;
% % filter(:) = 1;
% xin = xin .* filter;

%$$ 
xft   = FT1(xin);

%$$ xft   = FTlsq3(xin);

%xft   = FT9(xin);
%xft   = FT9_bak2_quite_gaugeInvariant(xin);
%xft   = FT9_bak1_quite_positive(xin);
%xft   = FTsimp2(xin);

nn    = size(xft,1);
nout  = floor((nn+1)/2);

omega = (0:nout-1).' / nn  * 2*pi  / dTime ;

outsize = xsize;
outsize(1) = nout;
xftpos = zeros(outsize);
xftneg = zeros(outsize);

xftpos(:,:) = xft(1:nout,:);
xftneg(:,:) = [xft(1,:);flipud(xft(nn-nout+2:nn,:))];

end



function [xft] = FT1(xin)
halfEndMode  = 3;
isize        = size(xin);
nnn          = isize(1);
nmod         = mod(nnn,2);
osize        = isize;
if 1==0
  xin(:,:) = xin(:,:) - xin(end,:);
elseif 1==1
  xin(:,:) = xin(:,:) - xin(1,:);
end

switch halfEndMode
  case 1
    xft = fft(xin,[],1);
  case 2
    HalfEnds = @(x) [x(1,:)/2;x(2:end-1,:);x(end,:)/2];
    xft = reshape(fft(HalfEnds(xin),[],1),osize);
  case 3
    QQ       = 1;
    nPad     = (2*QQ-1)*nnn-1;
    zsize    = isize;
    zsize(1) = nPad;
    nTot      = nnn+nPad;
    osize(1) = nTot;
    xWt      = xin;
    xWt(1,:) = xWt(1,:) * 0.5;
    xWt(nnn,:) = xWt(nnn,:) * 0.5;
    zo       = zeros(zsize);
    xWt = [xWt;zo];
    xft = reshape(fft(xWt(:,:),[],1),osize);
  otherwise
    error('not supported')
end
end  % function FT1




%       function y = solveFun(x,tflag,fF,fT)
%         if isequal(tflag,'transp')
%           y = fT(x);
%         else
%           y = fF(x);
%         end
%       end
%       tol = 1e-8;
%       maxit = 100;
%       lsqrFun = @(fun,rhs,guess) lsqr(fun,rhs,tol,maxit,[],[],guess);
%       passFun = @(x,tflag) reshape(solveFun(reshape(x,osize),tflag,ifftFun,ifftTra),isize);
%       y = lsqrFun(passFun,yguess);




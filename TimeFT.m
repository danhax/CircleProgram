
function [xftpos,xftneg,omega] = TimeFT(xin,dTime)  % dTime just for omega


xsize = size(xin);
% nin = xsize(1);

% % filter = (nin-1 : -1 : 0)' / (nin-1);
% filter = (nin : -1 : 1)' / nin;
% % filter(:) = 1;
% xin = xin .* filter;

%$$ 
xft   = FT1(xin);
%$$ xft   = FTlsq1(xin);

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





function y = PrimeAbove(xin) %,divisible)
x = xin+1;
while true
  y = primes(x);
  y = y(end);
  if y<xin
    x = x+1;
  else
    break
  end
end
% if mod(y+1,divisible)~=0
%   y = PrimeAbove(y,divisible);
% end
end


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




function hv = Magnus5(starttime,dTime,Ham,nTime)

nFun = 3;

Ham = @(x) -1i*dTime*Ham(x);

[quadPts,quadWts] = GaussJacobi(nTime,0,0) ;
quadPts = (quadPts(:) + 1) / 2 ;  % now 0,1
quadWts = quadWts(:) / 2 ;        % now for 0,1 (interval 1)

times = starttime + quadPts * dTime;

quadPts = quadPts - 1/2 ;  % now -1/2,1/2

endtime = starttime + dTime;

for itime          = 1:nTime
  thisham          = Ham(times(itime));
  if itime == 1
    hamdim         = size(thisham,1);
    hams           = zeros(hamdim,hamdim,nTime);
  end
  hams(:,:,itime)  = thisham;
end

% funs = quadPts.^(0:nFun-1);

  function y = myLegendre(n,x)
    leg = legendre(n,x*2);
    leg = leg(1,:);
    leg = reshape(leg,size(x));    
    y   = leg;
  end

allfuns = zeros(nTime,nFun);
for ifun = 1:nFun
  allfuns(:,ifun) = myLegendre(ifun-1,quadPts);
end

fixEnds = 1==0;

if fixEnds
  ham0 = Ham(starttime)  ;
  ham1 = Ham(endtime)   ;
  hamAvg = (ham1+ham0)/2 ;
  hamSlope = (ham1-ham0) ;
  hams = hams - hamAvg - hamSlope.*reshape(quadPts,1,1,nTime);
  Ps = @(x) mod(x,2);
    
  ii = Ps(2:nFun-1)+1;
  jj = 1:nFun-2;
  kk = -1 * ones(1,nFun-2);
  funtrans = full([sparse(ii,jj,kk,2,nFun-2); eye(nFun-2)]);
else
  funtrans = eye(nFun);
end

funs = allfuns * funtrans;

% b * funs' = hams

% bh = reshape(hams,hamdim^2,nTime) / funs';

swt = sqrt(quadWts(:)');
bh = (reshape(hams,hamdim^2,nTime).*swt) / (funs'.*swt);

if fixEnds
  bh = bh * funtrans';
end
bh = reshape(bh, hamdim,hamdim,nFun);
if fixEnds
  bh(:,:,1) = bh(:,:,1) + hamAvg;
  bh(:,:,2) = bh(:,:,2) + hamSlope /2 ;  %  /2 because 2nd jacobi poly is 2x
end

j0 = bh(:,:,1);
j1 = bh(:,:,2);
j2 = bh(:,:,3);

j01          = j0 * j1 - j1 * j0;
j02          = j0 * j2 - j2 * j0;
j12          = j1 * j2 - j2 * j1;

j0_01        = j0 * j01 - j01 * j0;
j0_02        = j0 * j02 - j02 * j0;
j1_01        = j1 * j01 - j01 * j1;
j0__0_01     = j0 * j0_01 - j0_01 * j0;

hv = j0 ...
  - 1/6 * j01 ...
  - 1/30 * j12 ...
  + 1/60 * j0_02 ...
  - 1/60 * j1_01 ...
  + 1/360 * j0__0_01 ;

hv = hv / (-1i) / dTime;

end

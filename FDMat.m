
function [fdmat,sdmat] = FDMat(nPts,ord)

% n = 2*ord+1;

pts = -ord:ord;

[fdvec,sdvec] = FDVec(ord);

fdmat = zeros(nPts,nPts);
sdmat = zeros(nPts,nPts);
for ipt = 1:nPts
  iind = mod(pts + ipt - 1, nPts) + 1;
  fdmat(ipt,iind) = fdvec;
  sdmat(ipt,iind) = sdvec;
end

end



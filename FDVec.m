
function [fdvec,sdvec,tdvec] = FDVec(ord)

n = 2*ord+1;
pts = -ord:ord;
fac = ord^((ord-1)/ord);
pts = pts / fac;

eqns = zeros(n,n);

for ieqn = 1:n
  eqns(ieqn,:) = pts.^(ieqn-1);
end

fdvec = eqns \ [0;1;zeros(n-2,1)];

sdvec = eqns \ [0;0;2;zeros(n-3,1)];

tdvec = eqns \ [0;0;0;6;zeros(n-4,1)];

fdvec = fdvec / fac;
sdvec = sdvec / fac^2;
tdvec = tdvec / fac^3;

end



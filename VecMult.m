
function y = VecMult(x,fdvec,ord)
m = 2*ord+1;
if numel(fdvec)~=m
  error('ack')
end
% n = numel(x);
n = size(x,1);
if n < m
  error('oops toosmall')
end
y         = zeros(size(x));
for i     = 1:m
  iind    = mod((0:n-1) + (i - ord - 1), n) + 1 ;
  y(:,:)  = y(:,:) + fdvec(i) * x(iind,:);
end
end



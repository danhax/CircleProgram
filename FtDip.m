function varargout = FtDip(dTime,varargin)
if nargin-1 ~= nargout
  error('oops')
end
varargout = cell(nargout,1);

%if numel(size(varargin{1})) > 3
%  error('ack checkme')
%end

% [vLen,nDers,nVec]     = size(varargin{1});

for iarg  = 1:nargout
  [varargout{iarg}{1},varargout{iarg}{2}] = TimeFT(varargin{iarg},dTime);
  
  %   [boogie{1,1,1,1},boogie{1,1,1,2}] = TimeFT(varargin{iarg},dTime);
  %   boogie2 = cell2mat(boogie);
  %   boogie3 = permute(boogie2,[1 4 2 3]);  % now vLen,2,nDers,nVec
  %   boogie4 = mat2cell(boogie3,vLen,[1;1],ones(nDers,1),nVec);
  %   varargout{iarg} = reshape(boogie4,2,nDers);
end
end





function [xftpos,xftneg,omega] = TimeFT(varargin)
  [xftpos,xftneg,omega] = TimeFT1(varargin{:})  ;
end

% function [xftpos,xftneg,omega] = TimeFT(xin,dTime,iDer)  % dTime just for omega
%   [xftpos,xftneg,omega] = TimeFT1(xin,dTime,iDer)  ;
%   % [xftpos,xftneg,omega] = TimeFT2(xin,dTime)  ;
% end

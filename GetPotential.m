function [v] = GetPotential(problemOption,x,~,FLIPPARITY)

%$$      v = (-sin(2*x)/4 - cos(x)/5) * 90 / pi ;

if problemOption == 0
  if FLIPPARITY
    v = (-sin(2*x)/4 + cos(x)/5) * 90 / pi ;
  else
    v = (-sin(2*x)/4 - cos(x)/5) * 90 / pi ;
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

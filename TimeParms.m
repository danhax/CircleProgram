
function [dTime,nTime] = TimeParms(omega,timeRes,duration,fdOrder,longerFac)
% there are nTime-1 steps

dTime = 1/omega * timeRes;

nTime = ceil(duration/dTime) + 1;

dTime = duration/(nTime-1);

% ensure differentiability at end!
nTime = nTime + fdOrder ;

% go for longer
nTime = ceil(nTime * longerFac);

end


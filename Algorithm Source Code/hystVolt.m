function [HV]= hystVolt(curSeq,Mt,Q)
%%% INPUT %%%
% curSeq: current input
% Lbest: MODEL PARAMETERS: Plett model
% TimeIncrement: if sampling rate is 1s then value is 1/3600 (hour)
% Q: the capacity of the cell

%%% OUTPUT %%%
% HV: hysteresis voltage

Mo = 8e-9;%Lbest(1);
% Mt = Lbest(2);
gamma = 10;%Lbest(3);
TimeIncrement=1/3600;


eta = 1; % coulombic efficiency

% calculate instantaneous hysteresis using the Plett equation
s = [];
for k = 1:length(curSeq) 
    if curSeq(k) ~= 0
        s(k) = sign(curSeq(k));
    elseif curSeq(k) ==0 && k~=1
        s(k) = s(k-1);
    elseif curSeq(k) ==0 && k==1
        s(k) = 0;
    end
end

% With the following Plett equation calculate the "H" parameter.
H(1) = 0;
for k = 1:length(curSeq)-1 
    H(k+1) = exp(-abs((eta*curSeq(k)*gamma*TimeIncrement)/Q))*H(k)-( ...
        1-exp(-abs((eta*curSeq(k)*gamma*TimeIncrement)/Q)))*sign(curSeq(k));
end

% calculate the hysteresis voltage
HV = zeros(size(curSeq));
for k = 1:length(curSeq)
    HV(k) = -H(k)*Mt(k)+Mo*s(k);
end

end
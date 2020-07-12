% RMT Predictions for ND vs INR

cosq =gencos(vm,vi);
sinq = sqrt(1-cosq);
tanq = sinq/sqrt(cosq);

NDrm = NDens*NDres;

% NDens --> Ensemble Notchdepth 

NDres = abs((INR*(tanq)*sqrt(c))-(1+c)).^2;



if N*INR > (sinq)^2 & c<=1 
%Breakpoint 1
if INR == 1/(N*(sinq)^2)
    % Slope becomes -2
    
    
end
%Breakpoint 2
if INR == (1+c)^2/(c*(tanq)^2)
    %Slope becomes -1
    
end

end

if c > 1
    
    if  INR == sqrt(c)/N
        
        % phase transition threshold
        % Rolls off slope -1
        
    end
end






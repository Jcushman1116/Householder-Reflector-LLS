function [T_value] = Chebyshev(x,d)
    n = length(x); 
    T_value = zeros(n,d+1); 

    T_value(:,1) = 1; % T0 = 1 
    if d >= 1
        T_value(:,2) = x; % T1 = x 
    end 
    
    % Define the recurrance relationsip 
    for k =2:d
        T_value(:,k+1) = 2*x .* T_value(:,k) - T_value(:,k-1);   
    end 

return

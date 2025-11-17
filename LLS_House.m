function [x,R] = LLS_House(A,b)

%Input Variable 
% - A: Nonsingular Matrix with n x k dimensions 
% - b: vector with n x 1 dimensions 

%Output Variable
% - x: Solution Vector with k x 1 dimensions 

% The solution is either exact OR a projection that numerically
% approximates the true solution onto the Column space of A 

% ~ represents an empty variable slot for n 
% efficient to only store columns for computations 
[~,k] = size(A); 

% Householder Reflector Loop
 for j = 1:k
     % Compute algebraic components of H Reflector 
     % Grab active column 
     a = A(j:end,j);
     % Initialize and set e1 basis vector 
     e1 = zeros(length(a),1); 
     e1(1) = 1; 
     % Formula Decomposition for vector in Householder Reflector 
     v = a + sign(a(1))*norm(a)*e1; 
     v = v/ norm(v);

     % Apply Reflectors and Update A and b 
     % Note v' represents the transpose of the vector v 
     A(j:end,j:end) = A(j:end,j:end) - 2*v*(v'*A(j:end,j:end));
     b(j:end) = b(j:end) - 2*v*(v'*b(j:end));
 end 

 % Extraction of the Elements for computations
 % Note that R is Upper triangular 
 R = A(1:k,1:k);
  
 % Check for linear independence 
 if min(abs(diag(R))) < eps * max(abs(R(:))) * k %eps is machine precision
    error('Matrix A does not have full column rank');
 end
 c = b(1:k);
% d = b(k+1:end); % Contains residuals for undetermined systems

% Backsolve 
x = zeros(k, 1); % Initialize solution vector 
for i = k:-1:1
    if i == k
        x(i) = c(i)/R(i,i); 
    else 
        x(i) = (c(i) - R(i,i+1:k) * x(i+1:k)) / R(i,i);
    end 
end
return 


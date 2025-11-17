% Test Driver for Program 4 

% Notes: S.E. - Solution Error 
%        G.F. - Growth Factor 
%        R.E. - Residual Error 

rng('shuffle');
n_max = 30;      % maximum size of matrix
n_tests = 5;     % number of trials per size

for task = 1:7
    switch task
        case 1, label = 'General Trend Case 1: n = k';
        case 2, label = 'General Trend Case 2: n > k s.t. b in R(A)';
        case 3, label = 'General Trend Case 3: n > k s.t. b not in R(A)';
        case 4, label = 'Required Problem 1: Growth Factor Matrix';
        case 5, label = 'Required Problem 2: Arithmetic Mean';
        case 6, label = 'Required Problem 3: Monomial Basis for d = 1,2,3';   
        case 7, label = 'Required Problem 4: Chebyshev Basis for d = 1,2,3';     

    end

    fprintf('\n%s\n', label); % Display Table Setup 
    fprintf('  n |  k  | cond(A)  | Mean S.E.| Max S.E. | Mean G.F.| Max G.F  | Mean R.E.| Max R.E. \n');
    fprintf('----|-----|----------|----------|----------|----------|----------|----------|----------\n');

    for n = 5:n_max 

        % Initialize storage for each metric
        conditioningNumber = zeros(n_tests,1);
        solError = zeros(n_tests,1);
        residError = zeros(n_tests,1);
        growthFactor = zeros(n_tests,1);

        if task == 6 || task == 7 %initialze storage for polynomial cases 
            degrees = [1,2,3];
            D = length(degrees);
            cond_by_d = zeros(n_tests, D);
            sol_by_d = zeros(n_tests, D);
            resid_by_d = zeros(n_tests, D);
            growth_by_d = zeros(n_tests, D);
        end

        for trial = 1:n_tests
            switch task

                % Case 1: Square nonsingular
                case 1
                    A = rand(n,n) + eye(n);  % Adding Identity ensures diagonal dominance 
                    b = rand(n,1);           % Genrate random vector
                    x_true = A\b;            % Calculate Matlab true solution 

                % Case 2: Overdetermined, b in R(A)
                case 2
                    k = max(1, n-1);  % ensures k < n 
                    A = rand(n,k);    % Generate random Matrix 
                    A = A + eye(n,k); % Ensure Diag. Dominance
                    z = rand(k,1);    % Genrate random vector
                    b = A*z;          % Use random vector to generate b in R(A)
                    x_true = A\b;     % Calculate Matlab true solution 

                % Case 3: Overdetermined, inconsistent
                case 3 %  
                    k = max(1, n-1);   % ensures k < n 
                    A = rand(n,k);      % Generate random Matrix
                    A = A + eye(n,k);   % Ensure Diag. Dominance 
                    x0 = rand(k,1); % start with random vector
                    b1 = A*x0; % Solve for b1 vector
                    r = rand(n,1); % add pertubation r to ensure b not in R(A)
                    b = b1 + r; % b is now not in R(A)

                    x_true = A\b;
                % Case 4: Growth factor test
                case 4
                    % Initialize Matrix Storage
                    A = zeros(n,n);
                    for i = 1:n % Loop to Generate Special Growth Matrix 
                        A(i,i) = 1;
                        A(i,n) = 1;
                        for j = 1:i-1
                            A(i,j) = -1;
                        end
                    end
                    b = rand(n,1); % Random vector 
                    x_true = A\b; % Matlab exact/projection sol

                
                % Case 5: Arithmetic mean computation
                case 5
                    A = ones(n,1);      % Vector of 1's
                    b = rand(n,1);      % random b vector 
                    x_true = mean(b);   % arithmetic mean true solution 
                
                % Case 6: Polynomial Regression
                case 6
                    b = rand(n,1); 
                    x_i = linspace(-1,1,n).';  % n random points between 0 and 1

                    % Loop over polynomial degrees to build design matrices A_d = [1 x x^2 ... x^d]
                    for m = 1:length(degrees)
                        d = degrees(m);
                        k = d+1; 
                        if k > n % Skip degrees 
                            continue;
                        end

                        A_d = zeros(n,k); % Build monomial matrix 
                        for q = 0:d
                            A_d(:,q+1) = x_i.^q;
                        end 
                        
                        % Gram Matrix Structure 
                        if trial == 1 && n == 10
                            G = A_d' * A_d;
                            fprintf('\n  Gram matrix for Monomial basis (n=%d, d=%d):\n', n, d);
                            disp(G);
                        end

                        %Algorithm Call 
                       [x_comp_d, R_d] = LLS_House(A_d,b); 
                       % Matlab 
                        x_true_d = A_d \ b;

                        % Condition Number per d
                       cond_by_d(trial,m) = cond(A_d); 

                       % Solution Error 
                       if norm(x_true_d) < 1e-14
                            sol_by_d(trial, m) = norm(x_comp_d,2);
                        else
                            sol_by_d(trial, m) = norm(x_true_d - x_comp_d, 2) / norm(x_true_d,2);
                        end

                        % Residual error 
                        resid_by_d(trial, m) = norm(b - A_d*x_comp_d, 2) / norm(b,2);

                        % Growth factor 
                        growth_by_d(trial, m) = max(abs(R_d(:))) / max(abs(A_d(:)));
                    end
                    continue; % go to the next trial

                %Chebyshev Basis 
                case 7
                    b = rand(n,1);
                    
                    % Process uniform points
                    for m = 1:D
                        d = degrees(m);
                        k = d+1;
                        
                        % Uniform points
                        x_i = linspace(-1,1,n).';
                        B = Chebyshev(x_i, d);
                        
                        % Gram matrix analysis (only once)
                        if trial == 1 && n == 10
                            G = B' * B;
                            fprintf('\n  Gram matrix for Chebyshev basis with uniform points (n=%d, d=%d):\n', n, d);
                            disp(G);
                        end
                        
                        a_true = B \ b;
                        [a_comp, R_d] = LLS_House(B,b);
                        
                        cond_by_d(trial,m) = cond(B);
                        sol_by_d(trial,m) = norm(a_true - a_comp) / norm(a_true);
                        resid_by_d(trial,m) = norm(b - B*a_comp) / norm(b);
                        growth_by_d(trial,m) = max(abs(R_d(:))) / max(abs(B(:)));
                    end
                    
                    % Print uniform results
                    if trial == n_tests  % After all trials complete
                        fprintf('\n UNIFORM POINTS \n');
                        meanCond_d = mean(cond_by_d, 1);
                        meanSol_d = mean(sol_by_d, 1);
                        maxSol_d = max(sol_by_d, [], 1);
                        meanRes_d = mean(resid_by_d, 1);
                        maxRes_d = max(resid_by_d, [], 1);
                        meanGrowth_d = mean(growth_by_d, 1);
                        maxGrowth_d = max(growth_by_d, [], 1);
                        
                        for dd = 1:length(degrees)
                            d = degrees(dd);
                            k_display = d + 1;
                            fprintf('%3d | %3d | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e   <-- d=%d (uniform)\n',n, k_display, meanCond_d(dd), meanSol_d(dd), maxSol_d(dd), meanGrowth_d(dd), maxGrowth_d(dd), meanRes_d(dd), maxRes_d(dd), d);
                        end
                    end
                    
                    % Reset arrays for Chebyshev root points
                    cond_by_d_cheb = zeros(n_tests, D);
                    sol_by_d_cheb = zeros(n_tests, D);
                    resid_by_d_cheb = zeros(n_tests, D);
                    growth_by_d_cheb = zeros(n_tests, D);
                    
                    % Process Chebyshev root points
                    for m = 1:D
                        d = degrees(m);
                        k = d+1;
                        
                        if n < k
                            continue;  % Skip if not enough points
                        end
                        
                        % Use the k Chebyshev roots 
                        x_i = zeros(k,1);
                        for i = 0:k-1
                            x_i(i+1) = cos((2*i+1)*pi/(2*k));
                        end
                        
                        % use only first k elements of b
                        b_reduced = b(1:k);
                        
                        B = Chebyshev(x_i, d);
                        
                        % Gram matrix analysis
                        if trial == 1 && n == 10
                            G = B' * B;
                            fprintf('\n  Gram matrix for Chebyshev basis with root points (n=%d, d=%d):\n', n, d);
                            disp(G);
                        end
                        
                        a_true = B \ b_reduced;
                        [a_comp, R_d] = LLS_House(B, b_reduced);
                        
                        cond_by_d_cheb(trial,m) = cond(B);
                        sol_by_d_cheb(trial,m) = norm(a_true - a_comp) / norm(a_true);
                        resid_by_d_cheb(trial,m) = norm(b_reduced - B*a_comp) / norm(b_reduced);
                        growth_by_d_cheb(trial,m) = max(abs(R_d(:))) / max(abs(B(:)));
                    end
                    
                    % Print Chebyshev root results
                    if trial == n_tests
                        fprintf('\n  CHEBYSHEV ROOT POINTS \n');
                        meanCond_d = mean(cond_by_d_cheb, 1);
                        meanSol_d = mean(sol_by_d_cheb, 1);
                        maxSol_d = max(sol_by_d_cheb, [], 1);
                        meanRes_d = mean(resid_by_d_cheb, 1);
                        maxRes_d = max(resid_by_d_cheb, [], 1);
                        meanGrowth_d = mean(growth_by_d_cheb, 1);
                        maxGrowth_d = max(growth_by_d_cheb, [], 1);
                        
                        for dd = 1:length(degrees)
                            d = degrees(dd);
                            k_display = d + 1;
                            fprintf('%3d | %3d | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e   <-- d=%d (Cheb roots)\n',k_display, k_display, meanCond_d(dd), meanSol_d(dd), maxSol_d(dd),meanGrowth_d(dd), maxGrowth_d(dd), meanRes_d(dd), maxRes_d(dd), d);
                        end
                    end
                    continue
            end 
        
            % Store original A for growth factor calc
            A_orig = A;

            % Algorithm Sol vector via  LLS_House solver
            [x_comp, R] = LLS_House(A, b);

            % Conditioning Number Calculation 
            conditioningNumber(trial) = cond(A_orig);
            if norm(x_true) < 1e-14 % Set Tolerance 
                solError(trial) = norm(x_comp,2);
            else    
                solError(trial) = norm(x_true - x_comp,2) / norm(x_true,2);
            end

            % Calculate Residual error
            residError(trial) = norm((b - A_orig*x_comp),2) / norm(b,2);

            % Calculate Growth factor 
            growthFactor(trial) = max(abs(R(:))) / max(abs(A_orig(:)));

        end 
 
        % Compute metric statistics
        meanSolErr = mean(solError);
        maxSolErr = max(solError);
        meanResErr = mean(residError);
        maxResErr = max(residError);
        meanGrowth = mean(growthFactor);
        maxGrowth = max(growthFactor);
        meanCond = mean(conditioningNumber);

        % Summary Table 
        k_display = size(A_orig,2);
        fprintf('%3d | %3d | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e\n', n, k_display, meanCond, meanSolErr, maxSolErr, meanGrowth, maxGrowth, meanResErr, maxResErr);

        % Create Statistics per d for cases 6 and 7 
    if task == 6 || task == 7
        meanCond_d = mean(cond_by_d, 1);
        meanSol_d = mean(sol_by_d, 1);
        maxSol_d = max(sol_by_d, [], 1);
        meanRes_d = mean(resid_by_d, 1);
        maxRes_d = max(resid_by_d, [], 1);
        meanGrowth_d = mean(growth_by_d, 1);
        maxGrowth_d = max(growth_by_d, [], 1);

        % Print one line per degree 
        for dd = 1:length(degrees)
            d = degrees(dd);
            k_display = d + 1;
            fprintf('%3d | %3d | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e | %1.2e   <-- d=%d\n',n, k_display, meanCond_d(dd), meanSol_d(dd), maxSol_d(dd), meanGrowth_d(dd), maxGrowth_d(dd), meanRes_d(dd), maxRes_d(dd), d);
        end

     end
    end
end

%  Polynomial comparison between monomial and Chebyshev
disp('Polynomial Comparison: Monomial vs Chebyshev');
n = 20; % 20 data points 
d = 2;  % Test on a quadratic 
x_i = linspace(-1, 1, n)'; % 20 random points between -1 and 1
b = sin(pi * x_i) + 0.1*randn(n,1); %random vector with pertubation 

% Monomial basis Matrix
A_mono = zeros(n, d+1);
for q = 0:d
    A_mono(:,q+1) = x_i.^q;
end
% Solution vector 
alpha_mono = LLS_House(A_mono, b);

% Chebyshev basis
B_cheb = Chebyshev(x_i, d);
alpha_cheb = LLS_House(B_cheb, b);

% Evaluate both at test points
x_test = linspace(-1, 1, 100)'; % 100 random points between -1 and 1 
y_mono = zeros(length(x_test), 1); % Initialize & evaluate monomial vector 
for q = 0:d
    y_mono = y_mono + alpha_mono(q+1) * x_test.^q; 
end

%Initialize test matrix 
T_test = Chebyshev(x_test, d);
% evaluate cheb vector 
y_cheb = T_test * alpha_cheb;

% Find the difference between the 2 evaluation vector 
fprintf('Max difference: %.2e\n', max(abs(y_mono - y_cheb)));
fprintf('Condition number - Monomial: %.2e\n', cond(A_mono));
fprintf('Condition number - Chebyshev: %.2e\n', cond(B_cheb));

% Note about the results: the evaluation error is small meaning that we can
% use either basis to find the correct coefficient for our regression 

% AKA both can be used to solve the LLS





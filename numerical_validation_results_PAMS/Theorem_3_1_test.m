%%%%%%%%%%%%%%%%%%%%%%%%%
%% Testing Theorem 3.1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
p1_lower = 1;
p1_upper = 3;
p2_lower = 1;
p2_upper = 3;
alpha_range_shift = 2;  % alpha will be from p to p + alpha_range_shift
nu1_range = 1:3;
nu2_range = 1:3;
trunc = 200;  % Bound on |kappa| in hypergeometric functions of matrix argument
numSamples = 100000000;  % Number of Monte Carlo samples

% Set random seed for reproducibility
rng(2); % Use a specific seed for reproducibility

% Start the global timer
global_start_time = tic;

% Pre-generate random positive definite matrices for each p1 and p2
Sigma_matrices = containers.Map; % Store matrices in a dictionary for each (p1, p2)

for p1 = p1_lower:p1_upper
    for p2 = max(p1, p2_lower):p2_upper
        p = p1 + p2;  % Calculate total dimension
        % Generate a random positive definite matrix Sigma using a Wishart distribution
        df = p; % Degrees of freedom for Wishart distribution
        scaleMatrix = eye(p); % Identity matrix as scale matrix for Wishart
        Sigma = wishrnd(scaleMatrix, df); % Random p x p positive definite matrix
        
        % Normalize the matrix so that the eigenvalues are between 0 and 1
        Sigma = Sigma / (trace(Sigma) + eps); % Normalize by trace to ensure eigenvalues between 0 and 1
        
        % Store the matrix and the corresponding eigenvalues
        key = sprintf('%d_%d', p1, p2);
        Sigma_matrices(key) = Sigma;
    end
end

% Initialize results array
results = [];

% Start a parallel pool (if not already started)
pool = gcp('nocreate');  % Check for existing parallel pool
if isempty(pool)
    pool = parpool;  % Create a new parallel pool with default settings
end

% Compare the Monte Carlo expectation on theL LHS with the analytical expression on the RHS
for p1 = p1_lower:p1_upper
    for p2 = max(p1, p2_lower):p2_upper
        p = p1 + p2;  % Calculate total dimension
        % Retrieve the pre-generated Sigma matrix
        Sigma = Sigma_matrices(sprintf('%d_%d', p1, p2));
        
        % Block matrices
        Sigma11 = Sigma(1:p1, 1:p1);
        Sigma12 = Sigma(1:p1, p1+1:p);
        Sigma21 = Sigma(p1+1:p, 1:p1);
        Sigma22 = Sigma(p1+1:p, p1+1:p);

        % A = Sigma^(-1)/2
        A = inv(Sigma) / 2;
        A11 = A(1:p1, 1:p1);
        A12 = A(1:p1, p1+1:p);
        A22 = A(p1+1:p, p1+1:p);

        % P matrix
        P = A11^(-1/2) * A12 * A22^(-1/2);
        PPt = P * P'; % P*P^T

        for alpha = p:(p + alpha_range_shift)
            for nu1 = nu1_range
                for nu2 = nu2_range
                    % Parallel Monte Carlo simulation for the left-hand side
                    mcSamples = zeros(numSamples, 1);
                    parfor i = 1:numSamples
                        X = wishrnd(Sigma, alpha);
                        mcSamples(i) = det(X(1:p1, 1:p1))^nu1 * det(X(p1+1:p, p1+1:p))^nu2;
                    end
                    mcLHS = mean(mcSamples);

                    % Analytical expressions for the right-hand side
                    exactExpectation1 = det(2 * Sigma11)^nu1 * (gamma_mv(alpha/2 + nu1, p1) / gamma_mv(alpha/2, p1));
                    exactExpectation2 = det(2 * Sigma22)^nu2 * (gamma_mv(alpha/2 + nu2, p2) / gamma_mv(alpha/2, p2));
                    hypergeomTerm = mhg(trunc, 2, [-nu1, -nu2], alpha/2, eig(PPt));
                    RHS = exactExpectation1 * exactExpectation2 * hypergeomTerm;

                    % Compute the absolute error and relative error
                    absoluteError = abs(mcLHS - RHS);
                    relativeError = absoluteError / abs(mcLHS);

                    % Store the result in the table
                    results = [results; p1, p2, alpha, nu1, nu2, {sprintf('%.6f ', eig(PPt))}, mcLHS, RHS, absoluteError, relativeError];

                    % Display progress
                    fprintf(['Computed for p1=%d, p2=%d, alpha=%.2f, nu1=%.2f, nu2=%.2f:\n', ...
                             'Monte Carlo LHS = %.6f, Analytical RHS = %.6f, Abs Error = %.6f, Rel Error = %.6f\n\n'], ...
                             p1, p2, alpha, nu1, nu2, mcLHS, RHS, absoluteError, relativeError);
                end
            end
        end
    end
end

% Stop the global timer and display the total time taken
totalElapsedTime = toc(global_start_time);
fprintf('Total computation time: %.2f seconds.\n', totalElapsedTime);

% Display the table of results
disp('Results Table: [p1, p2, alpha, nu1, nu2, Eigenvalues, MC_LHS, Analytical_RHS, Absolute Error, Relative Error]');
disp(results);

% Define column names for the results
variableNames = {'p1', 'p2', 'alpha', 'nu1', 'nu2', 'Eigenvalues', 'MC_LHS', 'Analytical_RHS', 'AbsError', 'RelError'};

% Convert the results array to a table
resultsTable = array2table(results, 'VariableNames', variableNames);

% Write the table to an Excel file
writetable(resultsTable, 'Theorem_3_1_test.xlsx');

% Shut down the parallel pool (optional)
delete(pool);

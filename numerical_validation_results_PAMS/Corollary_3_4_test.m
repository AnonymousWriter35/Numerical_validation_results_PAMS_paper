%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Corollary 3.4 Test (Eq. 11) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
a_range = 0.3:0.3:3;   % Range for 'a' values
b_range = 0.3:0.3:3;   % Range for 'b' values
m_range = 1:3;         % Range for matrix size 'm'
num_eigen_points = 5;  % Number of points per eigenvalue range
trunc = 100;           % Truncation level for hypergeometric series

% Set random seed for reproducibility
rng(2);

% Function to define the range of c as a function of m
c_range = @(m) linspace((m-1)/2 + 0.3, 10, (m-1)/2 + 3.3);  % Range of c for each m

% Initialize the results array to store (a, b, m, c, validity)
results = [];

% Start the global timer
global_start_time = tic;

% Loop over a, b, m, and c ranges
for a = a_range
    for b = b_range
        for m = m_range
            % Get the c_range for the current value of m
            current_c_range = c_range(m);

            % Loop over c for the current value of m
            for c = current_c_range
                % Define the eigenvalue ranges for the current m
                eigenvalue_ranges = arrayfun(@(x) linspace(0.1, 0.9, num_eigen_points), 1:m, 'UniformOutput', false);

                % Create the grid of eigenvalue combinations
                [eigenvalue_grids{1:m}] = ndgrid(eigenvalue_ranges{:});
                eigenvalue_combinations = [];
                for j = 1:m
                    eigenvalue_combinations = [eigenvalue_combinations, eigenvalue_grids{j}(:)];
                end

                all_valid_for_case = true;  % Assume all inequalities will be valid for this case

                % Loop over all combinations of eigenvalues
                for i = 1:size(eigenvalue_combinations, 1)
                    % Extract the current eigenvalues
                    eigenvalues = eigenvalue_combinations(i, :);

                    % Compute the hypergeometric function term for the left-hand side
                    lhs_hypergeom = mhg(trunc, 2, [-a, -b], c, eigenvalues);

                    % Compute the right-hand side components
                    rhs_term_1 = 1 + det(eye(m) - diag(eigenvalues))^(c + a + b) * ...
                                 (mhg(trunc, 2, [], c, eig(a * b * diag(eigenvalues))) - 1);

                    % Check the inequalities
                    inequality_1 = lhs_hypergeom >= rhs_term_1;
                    inequality_2 = rhs_term_1 >= 1;

                    if ~(inequality_1 && inequality_2)
                        all_valid_for_case = false;
                        break;  % No need to check further if one combination fails
                    end
                end

                % Display progression in the console
                if all_valid_for_case
                    fprintf('Both inequalities are valid for a=%.2f, b=%.2f, m=%d, c=%.2f\n', a, b, m, c);
                else
                    fprintf('Inequalities failed for a=%.2f, b=%.2f, m=%d, c=%.2f\n', a, b, m, c);
                end

                % Save the result for this case (a, b, m, c) and validity
                results = [results; a, b, m, c, all_valid_for_case];
            end
        end
    end
end

% Convert results to a table
resultsTable = array2table(results, 'VariableNames', {'a', 'b', 'm', 'c', 'Valid'});

% Write the results to an Excel file
writetable(resultsTable, 'Corollary_3_4_test.xlsx');

% Stop the global timer and display the total time taken
totalElapsedTime = toc(global_start_time);
fprintf('Total computation time: %.2f seconds.\n', totalElapsedTime);

% Message when completed
fprintf('Results saved to Corollary_3_4_test.xlsx.\n');

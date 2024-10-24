function [estimate, se, ci_margin, weighted_mean] = weighthedRegres( x, y, s )
%
% x : Independent variable (predictor) data as a column vector
% y : Dependent variable (predictand) data as a column vector
% s : Standard deviations of the errors in y as a column vector

x = x(:);
y = y(:);
s = s(:);

% Compute weights
weights = 1 ./ (s.^2); % Inverse of the variance
weighted_mean = sum(weights .* y) / sum(weights);

% Fit the weighted linear regression model
model = fitlm(x(:), y(:), 'linear', 'Weights', weights);

% Display model summary
disp(model)

% Access coefficients and confidence intervals
coefficients = model.Coefficients;
% Calculate effective degrees of freedom for the residuals
residual_dof = model.DFE;           % Degrees of Freedom for Error
effective_dof = residual_dof;       % Using the residual degrees of freedom
% Display effective degrees of freedom
fprintf('Effective Degrees of Freedom (Residuals): %d\n\n', effective_dof);

% Calculate 95% confidence intervals using the model's standard errors
confidence_intervals = coefCI(model, 0.05); % 95% confidence intervals
estimate = [0; 0];
se = [0; 0];
ci_margin = [0; 0];
% Display coefficients with standard error and confidence intervals
fprintf('Coefficient Statistics (with ± error):\n');
for i = 1:height(coefficients)
    estimate(i) = coefficients.Estimate(i);
    se(i) = coefficients.SE(i);
    %
    % % Standard Error format
    fprintf('Coefficient for %s: %.4f ± %.4f (SE)\n', coefficients.Properties.RowNames{i}, estimate(i), se(i));
    % % Confidence Interval format
    lower_ci = confidence_intervals(i, 1);
    upper_ci = confidence_intervals(i, 2);
    ci_margin(i) = (upper_ci - lower_ci) / 2; % Margin for CI
    fprintf('Coefficient for %s: %.4f ± %.4f (CI: [%.4f, %.4f])\n', ...
        coefficients.Properties.RowNames{i}, estimate(i), ci_margin(i), lower_ci, upper_ci);
end


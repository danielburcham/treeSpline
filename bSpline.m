function [C] = bSpline(Q)
%{
DESCRIPTION
Fit closed cubic b-spline curve to thinned, ordered sectional 
points with acceptable error by varying the number of knots.

NOTES
Adapted from Khameneifar, F., Feng, H.Y. 2014. Airfoil profile
reconstruction under the uncertainty of inspection data points.
International Journal of Advanced Manufacturing Technology  71:675-683.

REQUIRED INPUTS
Q - mx2 matrix: thinned, ordered sectional points

OUTPUTS
C - mx2 matrix: final sectional coordinates computed from closed, cubic 
b-spline curve

Copyright 2025 Daniel C. Burcham

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
%}

% Remove duplicate points
Q = unique(Q,'rows','stable');

% Select a start/end point with minimum curvature
% Is this too local?
k = diff(Q(:,2),2);
[~,id] = min(abs(k));

q = Q([id:end 1:id-1],:);

E = 10; % RMSE
u = 0.0035; % Measurement uncertainty
p = 4; % B-spline order 4 (degree 3)
n = 6; % Initial control points

% Parameterize using chord length
t = zeros(size(q,1),1);
t(2:end) = unique(cumsum(sqrt(diff(q(:,1)).^2+diff(q(:,2)).^2))./...
    sum(sqrt(diff(q(:,1)).^2+diff(q(:,2)).^2)));

while 6*E > u
    % Generate desired final knot vector
    U = zeros(n+p,1);
    U(end-p+1:end) = 1;
    for j = 1:n-p
        v = size(q,1)/(n-p+1);
        k = fix(j*v);
        alpha = (j*v)-k;
        U(p+j) = (1-alpha)*t(k) + alpha*t(k+1);
    end
    
    N = spcol(U,4,t);
    N(end,n) = 1;
    
    % Obtain P through constrained optimization
    P = optimvar('P',n,2);
    initialPoint.P = zeros(size(P));
    problem = optimproblem;
    problem.Objective = fcn2optimexpr(@objectiveFcn,P,N,q);
    problem.Constraints.constraint1 = P(1,:) == P(end,:);
    problem.Constraints.constraint2 = P(2,:)-P(1,:) == P(end,:)-P(end-1,:);
    % problem.Constraints.constraint3 = 4*(P(1,:)-P(end-1,:)) == P(2,:)-P(end-2,:);
    options = optimoptions(problem,"Display","off");
    [solution,~,~] = solve(problem,initialPoint,'Options',options);
    
    % Solution
    C = N*solution.P;
    
    % Calculate RMSE
    E = rmse(C,q,"all");

    % Increment
    n = n+1;

end

% Create B-form spline
curve = spmak(U',solution.P');

% Compute perimeter to create a curve with even 1 mm spacing between points
t = linspace(0,1,sum(sqrt(sum(diff(C([2:end 1], :), 1, 1).^2, 2)))/0.001)';
C = fnval(curve,t)';
end

function objective = objectiveFcn(P,N,q)

objective = sum(vecnorm((N*P)-q,2,2).^2);

end
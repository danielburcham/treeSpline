function p = quad3DFitTaubin(x,y,z)
%{
DESCRIPTION
quad3DFitTaubin: General quadric surface fit

NOTES
Adapted from Taubin, G. 1991. Estimation of planar curves, surfaces and 
nonplanar space curves defined by implicit equations, with applications 
to edge and range image segmentation. IEEE Transactions on Pattern 
Analysis and Machine Intelligence  13:1115-1138.

REQUIRED INPUTS
x - mx1 vector: x coordinates of 3D points
y - mx1 vector: y coordinates of 3D points
z - mx1 vector: z coordinates of 3D points

OUTPUTS
p - 10x1 vector: parameter vector containing coefficients of implicit
quadric surface:
p(1)*x^2 + p(2)y^2 + p(3)*z^2 + p(4)*x*y + p(5)*x*z + p(6)*y*z + ...
p(7)*x + p(8)*y + p(9)*z + p(10) = 0

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

% Auxiliary variables
l = ones(numel(x),1);
o = zeros(numel(x),1);

% Data matrix
X = [ x.^2, y.^2, z.^2, x.*y, x.*z, y.*z, x, y, z, l ];

% Gradients
dx = [2*x, o, o, y, z, o, l, o, o, o];
dy = [o, 2*y, o, x, o, z, o, l, o, o];
dz = [o, o, 2*z, o, x, y, o, o, l, o];
dX = [dx; dy; dz];

% Scatter matrices
A = X'*X;
B = dX'*dX;

[V,D] = eig(A,B);
[~,ix] = min(abs(diag(D)));
p = V(:,ix);
end
function splines = treeSpline(pc,planes)
%{
DESCRIPTION
treeSpline: Tree stem profile reconstruction from unorganized noisy point
cloud data

NOTES
Adapted from Ghorbani, H., Khameneifar, F. 2021. Airfoil profile
reconstruction from unorganized noisy point cloud data. Journal of
Computational Design and Engineering  8(2):740-755.

REQUIRED INPUTS
pc - mx3 matrix: set of 3D points in point cloud
plane - 1x6 vector: [X0 Y0 Z0  DX1 DY1 DZ1  DX2 DY2 DZ2], with
   - (X0, Y0, Z0) is a point belonging to the plane
   - (DX1, DY1, DZ1) is a first direction vector
   - (DX2, DY2, DZ2) is a second direction vector
   The 2 direction vectors are normalized and orthogonal.

OUTPUTS
curves - m x 1 cell array: 2D sectional coordinates computed from closed,
cubic b-spline curve

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

splines = cell(size(planes,1),1);

for i = 1:size(splines,1)
    [s,Ns,NBs] = planeProj(pc,planes(i,:));
    Ns = normalFilter(Ns,NBs,planes(i,:));
    w = weights(s,Ns,NBs);
    [p,Np] = wlls(s,Ns,w);
    Q = reconstructProfile(p,Np);
    splines{i,1} = bSpline(Q);
    clearvars -except i pc splines planes
end

end
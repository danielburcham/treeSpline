function Nb = growRing(pc,p,m)
%{
DESCRIPTION
Expand directionally-balanced neighborhood of points - 
Algorithm 3 (Khameneifar and Feng 2017)
Progressively expand the directionally-balanced neighborhood of points 
closest to point, p, in a point cloud, pc, to obtain the
minimum number of neighboring points, m, for quadric surface fitting. 

NOTES
Adapted from Khameneifar, F., Feng, H.Y. 2017. Establishing a balanced
neighborhood of discrete points for local quadric surface fitting.
Computer-Aided Design  84:25-38.

REQUIRED INPUTS
pc - mx3 matrix: point cloud
p - 1x3 vector: point of interest
m - scalar: minimum number of points

OUTPUTS
Nb - mx3 matrix: m directionally-balanced set of neighboring points

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

Nb = neighborRing(pc,p);
while size(Nb,1) < m
    Nn = neighborRing(pc(~ismember(pc,Nb,'rows'),:),p);
    Nb = union(Nb,Nn,'rows');
end
end
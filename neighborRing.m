function N1 = neighborRing(pc,p)
%{
DESCRIPTION
Closest set of directionally-balanced neighboring points - 
Algorithm 2 (Khameneifar and Feng 2017)
Select the directionally-balanced neighborhood of points closest to point
p in a point cloud, pc. 

NOTES
Adapted from Khameneifar, F., Feng, H.Y. 2017. Establishing a balanced
neighborhood of discrete points for local quadric surface fitting.
Computer-Aided Design  84:25-38. 

REQUIRED INPUTS
pc - mx3 matrix: point cloud
p - 1x3 vector: point of interest

OUTPUTS
N1 - mx3 matrix: closest directionally-balanced set of neighboring points

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

davg = zeros(11,1);
Nb = cell(11,1); % Variable-sized neighborhood
beta = (2:-0.1:1)'; % Region of influence occupied by each point
for i=1:size(beta,1)
    Nb{i,1} = territoryClaim(pc,p,beta(i));
    dist = vecnorm(Nb{i,1}-p,2,2);
    davg(i)=mean(dist);
end
N1 = Nb{davg == min(davg),1};
end
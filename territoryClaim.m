function Nb = territoryClaim(pc,p,beta)
%{
DESCRIPTION
Territory claiming algorithm - Algorithm 1 (Khameneifar and Feng 2017)
Select a directionally-balanced neighborhood of points, Nb, surrounding 
point, p, in a point cloud slice, pcSlice, where each point in Nb 
occupies a distinct 3D region of influence defined by beta.

NOTES
Adapted from Khameneifar, F., Feng, H.Y. 2017. Establishing a balanced
neighborhood of discrete points for local quadric surface fitting.
Computer-Aided Design  84:25-38. 

REQUIRED INPUTS
pc - mx3 matrix: point cloud
p - 1x3 vector: point of interest
beta - scalar: parameter defining 3D region of influence

OUTPUTS
Nb - mx3 matrix: directionally-balanced set of neighboring points
NOTE: knnsearch largely governs the speed of code

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

[idx,dist] = knnsearch(pc,p,'K',31,'SortIndices',false);
q = sortrows([pc(idx(dist~=0),:),dist(dist~=0)'],4);
Nb = q(1,1:3);
for i=2:30
    a = vecnorm(Nb - (p*(1-(beta/2))+q(i,1:3)*(beta/2)),2,2);
    b = vecnorm(Nb - (q(i,1:3)*(1-(beta/2))+p*(beta/2)),2,2);
    R = (beta/2)*norm(p-q(i,1:3));
    if all(a > R) && all(b > R)
        Nb = union(Nb,q(i,1:3),'rows');
    end
end
end
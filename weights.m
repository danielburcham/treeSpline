function w = weights(s,Ns,NBs)
%{
DESCRIPTION
Assign weight factors, w, to each point based on the distance 
from its sectional neighborhood of points

NOTES
Adapted from Ghorbani, H., Khameneifar, F. 2021. Airfoil profile
reconstruction from unorganized noisy point cloud data. Journal of 
Computational Design and Engineering  8(2):740-755. 

REQUIRED INPUTS
s - mx2 matrix: points projected onto 2D section
Ns - mx2 matrix: normal vector set for sectional points
NBs - mx1 cell array: indices of sectional neighborhood for each point

OUTPUTS
w - mx1 vector: weighting factors for each sectional point

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

OF = zeros(size(s,1),1);
for i = 1:size(s,1)
    % Project neighboring points onto normal for each point
    dx = s(NBs{i,1},1) - s(i,1);
    dy = s(NBs{i,1},2) - s(i,2);
    tp = (dx .* Ns(i,1) + dy .* Ns(i,2)) ./ (Ns(i,1) .* Ns(i,1) + Ns(i,2) .* Ns(i,2));
    points = [s(i,1) + tp .* Ns(i,1), s(i,2) + tp .* Ns(i,2)];
    
    % Determine average distance between point and projected neighbors
    ds = mean(vecnorm(s(i,1:2)-points));

    % Determine average inner distance among projected neighbors
    b = nchoosek(1:size(NBs{i,1},1),2);
    Ds = mean(vecnorm(points(b(:,1),:)-points(b(:,2),:)));
    % Compute outlier factor
    OF(i,1) = ds/Ds;
end
sdout = std(OF);
w = exp(-1*OF(:,1)/sdout);
end
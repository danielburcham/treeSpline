function [p,Np] = wlls(s,Ns,w)
%{
DESCRIPTION
Weighted local least squares (Ghorbani and Khameneifar 2021)
Adaptively fit a local quadratic curve to sectional points using 
weighted local least squares

NOTES
Adapted from Ghorbani, H., Khameneifar, F. 2021. Airfoil profile
reconstruction from unorganized noisy point cloud data. Journal of 
Computational Design and Engineering  8(2):740-755.

REQUIRED INPUTS
s - mx2 matrix: points projected onto 2D section
Ns - mx2 matrix: normal vector set for sectional points
w - mx1 vector: weighting factors for each sectional point

OUTPUTS
p - mx2 matrix: points projected onto thin sectional curve estimated
from wlls regression
Np - mx2 matrix: improved normal vectors for each sectional point
estimated from wlls regression

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

output = [];

for i = 1:size(s,1)
    n = 3;
    rmse = 0;
    u = 0.00325;
    pos = zeros(size(s,1),1);
    for ia = 1:size(s,1)
        pos(ia,1) = sign(det([Ns(i,1) Ns(i,2);...
            s(ia,1)-s(i,1) s(ia,2)-s(i,2)]));
    end
    s1 = [s(pos < 0,1:2) w(pos < 0)];
    s2 = [s(pos > 0,1:2) w(pos > 0)];
    % Use normal vector to determine fitting orientation
    ang = atan2(Ns(i,1),-1.*Ns(i,2)) * (180/pi);
    % Compute rotation matrix for ang degree rotation about z-axis
    R1 = rotz(-ang);
    while rmse < u 
        [idx1,dist1] = knnsearch(s1(:,1:2),s(i,1:2),'K',round(n/2),...
            'NSMethod','kdtree','SortIndices',false);
        [idx2,dist2] = knnsearch(s2(:,1:2),s(i,1:2),'K',round(n/2),...
            'SortIndices',false);
        pts = [s1(idx1,1:3) dist1'; s2(idx2,1:3) dist2'];
        pts = sortrows(pts,4);
        % pts = [x_coord y_coord weight]
        pts = [s(i,1:2) w(i); pts(1:n,1:3)];
        % Transform subset of points for curve fitting
        rotPts = (R1*[pts(:,1:2) ones(size(pts,1),1)]')';
        % Fit local quadratic polynomial with weights
        A = [ones(size(rotPts,1),1) rotPts(:,1) rotPts(:,1).^2];
        b = rotPts(:,2);
        [p,~,mse] = lscov(A,b,rotPts(:,3));
        rmse = sqrt(mse);
        n = n+1;
    end
    % Calculate new points using final curve parameters
    nRotPts = zeros(size(rotPts,1),4);
    for ia = 1:size(rotPts,1)
        syms x
        poly = p(1) + p(2)*x + p(3)*x^2;
        distpoly = (x - rotPts(ia,1))^2 + (p(1) + p(2)*x + p(3)*x^2 - rotPts(ia,2))^2;
        r = solve(diff(distpoly));
        r(imag(double(r))~=0) = [];
        [~,idx] = min(double(subs(distpoly,r)));
        nRotPts(ia,1:2) = [double(r(idx)) double(subs(poly,r(idx)))];
        dY = matlabFunction(diff(poly,x));
        [N,D] = rat(feval(dY,nRotPts(ia,1)));
        if [0 -1]*[N; -1*D] >= 0
            nRotPts(ia,3:4) = [N -1*D]./norm([N -1*D]);
        else
            nRotPts(ia,3:4) = [-1*N D]./norm([-1*N D]);
        end
    end
    % Back-transform coordinates
    R2 = rotz(ang);
    % Keep original coordinates, new coordinates, and new normal vectors
    output = [output; pts(:,1:2) (R2*[nRotPts(:,1:2) ones(size(nRotPts,1),1)]')'...
        (R2*[nRotPts(:,3:4) ones(size(nRotPts,1),1)]')'];
end
output = output(:,[1:4,6:7]);
C = unique(output(:,1:2),'stable','rows');
p = zeros(size(C,1),2);
Np = zeros(size(C,1),2);
for i = 1:size(C,1)
    p(i,1:2) = mean(output(all(output(:,1:2) == C(i,:),2),3:4)); 
    NS = [output(all(output(:,1:2) == C(i,:),2),5:6)...
        zeros(size(output(all(output(:,1:2)==C(i,:),2)),1),1)];
    for ia = 1:size(NS,1)-1
        ids = nchoosek(1:size(NS,1),2);
        angs = atan2(vecnorm(cross(NS(ids(:,1),:),NS(ids(:,2),:),2),2,2),...
            dot(NS(ids(:,1),:),NS(ids(:,2),:),2));
        [~,I] = min(abs(angs));
        NS(ids(I,1),:) = NS(ids(I,1),:) + NS(ids(I,2),:);
        NS = NS(setdiff(1:size(NS,1),ids(I,2)),:);
    end
    Np(i,1:2) = NS(1:2)./norm(NS(1:2));    
end

end
function [s,Ns,NBs] = planeProj(pc,plane)
%{
DESCRIPTION
planeProj: Project points onto a plane using adaptive surface projection 
with a quadric surface fit to a balanced neighborhood of points in the 
vicinity of each point near the sectional plane

NOTES
Adapted from Khameneifar, F., Feng, H.Y. 2017. Extracting sectional
contours from scanned point clouds via adaptive surface projection.
International Journal of Production Research  55(15):4466-4480.

REQUIRED INPUTS
pc - mx3 matrix: set of 3D points in point cloud
plane - 1x6 vector: [X0 Y0 Z0  DX1 DY1 DZ1  DX2 DY2 DZ2], with
   - (X0, Y0, Z0) is a point belonging to the plane
   - (DX1, DY1, DZ1) is a first direction vector
   - (DX2, DY2, DZ2) is a second direction vector
   The 2 direction vectors are normalized and orthogonal.

OUTPUTS
s - mx2 matrix: points projected onto 2D section
Ns - mx2 matrix: normal vector set for sectional points
NBs - mx1 cell array: indices of sectional neighborhood for each point

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

% Unit normal vector for section plane
n = cross(plane(:,4:6),plane(:,7:9),2)'; 
n = n./norm(n);

% Determine signed distance between all points and plane
dist = -sum(bsxfun(@times, n', bsxfun(@minus, plane(:,1:3), pc)), 2);

% Select points within +/- k m of section plane
k = 0.06;
idx = find(abs(dist) <= k); 

% Pre-allocate arrays
s = nan([numel(idx), 3]); 
Ns = nan([numel(idx), 3]); 
NBs = cell(numel(idx), 1); 
keys = nan([numel(idx), 1]);

for i = 1:numel(idx)

    % Find balanced neighborhood for each point 
    Nb = growRing(pc(idx,:),pc(idx(i),:),75); 
    
    % Does the balanced neighborhood contain points above and below the 
    % section plane? 
    dist = -sum(bsxfun(@times, n', bsxfun(@minus, plane(:,1:3), Nb)), 2);
    if all(sign(dist) == 1|0) || all(sign(dist)== -1|0)
        continue
    else
        % Fit local quadric surface
        in = [Nb; pc(idx(i),:)];
        params = quad3DFitTaubin(in(:,1),in(:,2),in(:,3));

        % Normal vector for quadric surface evaluated at pc(idx(i),:)
        syms X Y Z
        F = params(1).*X.^2 + params(2).*Y.^2 + params(3).*Z.^2 + ...
            params(4).*X.*Y + params(5).*X.*Z + params(6).*Y.*Z + ...
            params(7).*X + params(8).*Y + params(9).*Z + params(10); 
        dX = matlabFunction(diff(F,X));
        dY = matlabFunction(diff(F,Y));
        dZ = matlabFunction(diff(F,Z));
        
        try
            [~,foot] = pqdist(pc(idx(i),:)',params);
            np = [dX(foot(1),foot(2),foot(3));...
                dY(foot(1),foot(2),foot(3)); ...
                dZ(foot(1),foot(2),foot(3))]; 
        catch
            np = [dX(pc(idx(i),1),pc(idx(i),2),pc(idx(i),3));...
                dY(pc(idx(i),1),pc(idx(i),2),pc(idx(i),3));...
                dZ(pc(idx(i),1),pc(idx(i),2),pc(idx(i),3))];
        end
        
        % Unit orthogonal and normal vectors for trajectory plane
        nz = np./norm(np); 
        ny = cross(n,np); 
        ny = ny./norm(ny); 
        nx = cross(ny,nz); 
        nx = nx./norm(nx); 
        
        % Coordinate transformation matrix
        T = [nx, ny, nz, pc(idx(i),:)'; 0, 0, 0, 1]; 
        
        % Transform quadric surface coefficients
        A = [params(1), 0.5*params(4), 0.5*params(5), 0.5*params(7);...
            0.5*params(4), params(2), 0.5*params(6), 0.5*params(8);...
            0.5*params(5), 0.5*params(6), params(3), 0.5*params(9);...
            0.5*params(7), 0.5*params(8), 0.5*params(9), params(10)]; 
        Ap = T'*A*T; 
        
        % Transform section plane coefficients
        p = [n(1,1); n(2,1); n(3,1); -1*(n(1,1)*plane(1,1) + ...
            n(2,1)*plane(1,2) + n(3,1)*plane(1,3))]; 
        pil = T'*p; 
        
        % Solve intersection between plane and curve
        eqns = [Ap(1,1)*X^2 + Ap(3,3)*Z^2 + 2*Ap(3,1)*X*Z + ...
            2*Ap(4,1)*X + 2*Ap(4,3)*Z + Ap(4,4) == 0, pil(1)*X +...
            pil(3)*Z + pil(4) == 0]; 
        vars = [X,Z]; 
        [solX, solZ] = solve(eqns,vars); 
        
        % If both are equal to zero, should you just take the closer point?
        if all(~isreal(double([solX(1) solX(2) solZ(1) solZ(2)])))
            continue % No real solutions
        elseif atan2d(dot(cross(nz,n),cross(plane(4:6),plane(7:9))),dot(nz,n)) < 15
            continue 
        elseif pil(1)*solX(1) + pil(3)*solZ(1) + pil(4) == 0 && ...
            pil(1)*solX(2) + pil(3)*solZ(2) + pil(4) == 0
            % Two solutions - choose the point closer to p
            [~,id] = min([vecnorm([double(solX(1)); double(solZ(1))])...
                vecnorm([double(solX(2)); double(solZ(2))])]);
            newPoints = T*[double(solX(id)); 0; double(solZ(id)); 1];
        elseif pil(1)*solX(1) + pil(3)*solZ(1) + pil(4) == 0
            newPoints = T*[double(solX(1)); 0; double(solZ(1)); 1];
        elseif pil(1)*solX(2) + pil(3)*solZ(2) + pil(4) == 0
            newPoints = T*[double(solX(2)); 0; double(solZ(2)); 1];
        end
        s(i,:) = newPoints(1:3)';
        [~,keys(i)] = ismember(pc(idx(i),:),pc,'rows'); 
        [~,NBs{i,1}] = ismember(Nb,pc,'rows');
        
        % Project normal vector onto sectional plane
        A1 = [plane(4:6)' plane(7:9)'];
        Ns(i,:) = A1*((A1'*A1)\A1'*nz);
    end
end
s = s(~isnan(s(:,1)),:); 
Ns = Ns(~isnan(Ns(:,1)),1:2); 
keys = keys(~isnan(keys)); 
NBs = NBs(~cellfun('isempty',NBs),1); 
for i = 1:numel(NBs)
    NBs{i,1} = NBs{i,1}(ismember(NBs{i,1},keys),:); 
    [~,NBs{i,1}] = ismember(NBs{i,1},keys);
end
% Orient all normal vectors outwards
centroid = mean(s(:,1:2));
for i = 1:size(s,1)
    v1 = s(i,1:2) - centroid; 
    [v,d] = eig(cov(s(NBs{i,1},1:2))); 
    v2 = v(:,[d(1,1) d(2,2)] == max([d(1,1) d(2,2)]))'; 
    ang1 = atan2d(dot(cross([v1 0],[v2 0]),cross(plane(4:6),plane(7:9))),...
        dot([v1 0],[v2 0])); 
    ang3 = atan2d(dot(cross([v1 0],[Ns(i,:) 0]),cross(plane(4:6),...
        plane(7:9))),dot([v1 0],[Ns(i,:) 0])); 
    if ang1 > 0
        ang2 = ang1 - 180; 
        if ang3 < ang2 || ang3 > ang1
            Ns(i,:) = -1.*Ns(i,:);
        end
    else 
        ang2 = ang1 + 180;
        if ang3 < ang1 || ang3 > ang2
            Ns(i,:) = -1.*Ns(i,:);
        end
    end
end
end
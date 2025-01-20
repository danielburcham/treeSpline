function Q = reconstructProfile(p,Np)
%{
DESCRIPTION
Order sectional profile points estimated from wlls 
regression and remove points with undesired connectivity deviating from 
underlying curve.

NOTES
Adapted from Ghorbani, H., Khameneifar, F. 2021. Airfoil profile
reconstruction from unorganized noisy point cloud data. Journal of 
Computational Design and Engineering  8(2):740-755.

REQUIRED INPUTS
p - mx2 matrix: points projected onto thin sectional curve estimated
from wlls regression
Np - mx2 matrix: improved normal vectors for each sectional point
estimated from wlls regression

OUTPUTS
Q - mx2 matrix: thinned, ordered sectional points

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

% Initial reconstruction (nodes ordered CCW)
E = zeros(size(p,1),2);
for i = 1:size(p,1)
    pos = zeros(size(p,1),1);
    for ia = 1:size(p,1)
        pos(ia,1) = sign(det([Np(i,1) Np(i,2);...
            p(ia,1)-p(i,1) p(ia,2)-p(i,2)]));
    end
    s1 = p(pos < 0,:);
    s2 = p(pos > 0,:);
    [idx1,~] = knnsearch(s1,p(i,:),'K',1);
    [idx2,~] = knnsearch(s2,p(i,:),'K',1);
    [~,Loc1] = ismember(s1(idx1,:), p(:,1:2), 'rows');
    [~,Loc2] = ismember(s2(idx2,:), p(:,1:2), 'rows');
    E(i,:) = [Loc1 Loc2];
end

% Imperfect node modification
% Find indices of rows with connections to the same preceding and trailing nodes
[C1,~,ic1] = unique(E,'rows','stable');
cts1 = [C1, accumarray(ic1,1)];
cts1 = cts1(cts1(:,3)>1,:);
ids = ismember(E,cts1(:,1:2),'rows');
ids1 = zeros(sum(cts1(:,3)),3);
ids1(:,1) = find(ids);
ids1(:,2:3) = E(ids1(:,1),:);
Tp = [-1.*Np(ids1(:,2),2) Np(ids1(:,2),1)]; % Tangent vector
% Edge vector
e = [p(ids1(:,1),1)-p(ids1(:,2),1) p(ids1(:,1),2)-p(ids1(:,2),2)]./...
    norm([p(ids1(:,1),1)-p(ids1(:,2),1) p(ids1(:,1),2)-p(ids1(:,2),2)]);
theta = atan2d(norm(cross([Tp zeros(size(Tp,1),1)],[e zeros(size(e,1),1)],2)),...
    dot([Tp zeros(size(Tp,1),1)],[e zeros(size(e,1),1)],2));
idx = zeros(size(cts1,1),1);
for i = 1:size(cts1,1)
    mx = max(theta(ismember(ids1(:,2:3),cts1(i,1:2),'rows')));
    idx(i) = ids1(theta==mx,1);
end
p(idx,:) = [];
Np(idx,:) = [];

E = zeros(size(p,1),2);
for i = 1:size(p,1)
    pos = zeros(size(p,1),1);
    for ia = 1:size(p,1)
        pos(ia,1) = sign(det([Np(i,1) Np(i,2);...
            p(ia,1)-p(i,1) p(ia,2)-p(i,2)]));
    end
    s1 = p(pos < 0,:);
    s2 = p(pos > 0,:);
    [idx1,~] = knnsearch(s1(:,1:2),p(i,:),'K',1);
    [idx2,~] = knnsearch(s2(:,1:2),p(i,:),'K',1);
    [~,Loc1] = ismember(s1(idx1,:), p(:,1:2), 'rows');
    [~,Loc2] = ismember(s2(idx2,:), p(:,1:2), 'rows');
    E(i,:) = [Loc1 Loc2];
end

% Find indices of rows with connections to the same preceding nodes
[C2,~,ic2] = unique(E(:,1),'stable');
cts2 = [C2, accumarray(ic2,1)];
cts2 = cts2(cts2(:,2)>1,:);
[~,ids] = ismember(E(:,1),cts2(:,1),'rows');
ids2 = zeros(sum(cts2(:,2)),3);
ids2(:,1) = find(ids);
ids2(:,2:3) = E(ids2(:,1),:);
Tp = [-1.*Np(ids2(:,2),2) Np(ids2(:,2),1)]; % Tangent vector
e = [p(ids2(:,1),1)-p(ids2(:,2),1) p(ids2(:,1),2)-p(ids2(:,2),2)]./...
    norm([p(ids2(:,1),1)-p(ids2(:,2),1) p(ids2(:,1),2)-p(ids2(:,2),2)]);
theta = atan2d(norm(cross([Tp zeros(size(Tp,1),1)],[e zeros(size(e,1),1)],2)),...
    dot([Tp zeros(size(Tp,1),1)],[e zeros(size(e,1),1)],2));
idx = zeros(size(cts2,1),1);
for i=1:size(cts2,1)
    [~,ix] = max(theta(ids2(:,2)==cts2(i,1)));
    tmp = ids2(ids2(:,2)==cts2(i,1));
    idx(i) = tmp(ix,1);
end
p(idx,:) = [];
Np(idx,:) = [];

E = zeros(size(p,1),2);
for i = 1:size(p,1)
    pos = zeros(size(p,1),1);
    for ia = 1:size(p,1)
        pos(ia,1) = sign(det([Np(i,1) Np(i,2);...
            p(ia,1)-p(i,1) p(ia,2)-p(i,2)]));
    end
    s1 = p(pos < 0,:);
    s2 = p(pos > 0,:);
    [idx1,~] = knnsearch(s1(:,1:2),p(i,:),'K',1);
    [idx2,~] = knnsearch(s2(:,1:2),p(i,:),'K',1);
    [~,Loc1] = ismember(s1(idx1,:), p(:,1:2), 'rows');
    [~,Loc2] = ismember(s2(idx2,:), p(:,1:2), 'rows');
    E(i,:) = [Loc1 Loc2];
end

% Find indices of rows with connections to the same trailing nodes
[C3,~,ic3] = unique(E(:,2),'stable');
cts3 = [C3, accumarray(ic3,1)];
cts3 = cts3(cts3(:,2)>1,:);
[~,ids] = ismember(E(:,2),cts3(:,1),'rows');
ids3 = zeros(sum(cts3(:,2)),3);
ids3(:,1) = find(ids);
ids3(:,2:3) = E(ids3(:,1),:);
Tp = [-1.*Np(ids3(:,2),2) Np(ids3(:,2),1)]; % Tangent vector
e = [p(ids3(:,1),1)-p(ids3(:,2),1) p(ids3(:,1),2)-p(ids3(:,2),2)]./...
    norm([p(ids3(:,1),1)-p(ids3(:,2),1) p(ids3(:,1),2)-p(ids3(:,2),2)]);
theta = atan2d(norm(cross([Tp zeros(size(Tp,1),1)],[e zeros(size(e,1),1)],2)),...
    dot([Tp zeros(size(Tp,1),1)],[e zeros(size(e,1),1)],2));
idx = zeros(size(cts3,1),1);
for i=1:size(cts3,1)
    [~,ix] = max(theta(ids3(:,3)==cts3(i,1)));
    tmp = ids3(ids3(:,3)==cts3(i,1));
    idx(i) = tmp(ix,1);
end
p(tmp(ix,1),:) = [];
Np(tmp(ix,1),:) = [];

E = zeros(size(p,1),2);
for i = 1:size(p,1)
    pos = zeros(size(p,1),1);
    for ia = 1:size(p,1)
        pos(ia,1) = sign(det([Np(i,1) Np(i,2);...
            p(ia,1)-p(i,1) p(ia,2)-p(i,2)]));
    end
    s1 = p(pos < 0,:);
    s2 = p(pos > 0,:);
    [idx1,~] = knnsearch(s1(:,1:2),p(i,:),'K',1);
    [idx2,~] = knnsearch(s2(:,1:2),p(i,:),'K',1);
    [~,Loc1] = ismember(s1(idx1,:), p(:,1:2), 'rows');
    [~,Loc2] = ismember(s2(idx2,:), p(:,1:2), 'rows');
    E(i,:) = [Loc1 Loc2];
end

orderP = zeros(size(p,1),1);
orderP(1) = 1;
for i = 2:size(orderP,1)
    orderP(i)=E(orderP(i-1),2);
end
Q = p(orderP,:);

end
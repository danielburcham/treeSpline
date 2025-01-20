function Ns = normalFilter(Ns,NBs,plane)
n = zeros(size(Ns,1),2);
for i = 1:size(Ns,1)
    NB = [Ns(NBs{i,1},:) zeros(size(NBs{i,1},1),1); Ns(i,:) 0];
    ids = nchoosek(1:size(NB,1),2);
    angs = atan2d(dot(cross(...
        NB(ids(:,1),:),NB(ids(:,2),:)),...
        repmat(cross(plane(4:6),plane(7:9)),size(ids,1),1),2),...
        dot(NB(ids(:,1),:),NB(ids(:,2),:),2));
    [~,I] = min(abs(angs));
    ni = NB(ids(I,1),:) + NB(ids(I,2),:);
    NB([ids(I,:)],:) = [];

    while ~isempty(NB)
        angs = atan2d(dot(cross(repmat(ni,size(NB,1),1),NB),...
            repmat(cross(plane(4:6),plane(7:9)),size(NB,1),1),2),...
            dot(repmat(ni,size(NB,1),1),NB,2));
        [~,I] = min(abs(angs));
        ni = ni + NB(I,:);
        NB(I,:) = [];
    end
    n(i,:) = ni(1:2)./norm(ni(1:2));
end
Ns = n;
end
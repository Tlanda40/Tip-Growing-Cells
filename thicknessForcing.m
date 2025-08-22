function [V, T, volume] = thicknessForcing(Vold, Told, poleColumn)
    % Step 1: Get free boundary
    V = Vold;
    T = Told;
    TR = triangulation(T, V);
    [F, V_bound] = freeBoundary(TR);
    
    % Step 2: Clean up vertex indexing
    [uniqueVerts, ~, ic] = unique(F(:));
    V_clean = V_bound(uniqueVerts, :);
    F_clean = reshape(ic, size(F));
    
    % Step 3: Create triangulation for boundary faces
    TR_surf = triangulation(F_clean, V_clean);
    faceAdj = TR_surf.neighbors;
    
    % Step 4: Efficient BFS to find connected components
    numFaces = size(F_clean, 1);
    componentId = zeros(numFaces, 1);
    currentComponent = 0;
    
    for i = 1:numFaces
        if componentId(i) > 0
            continue;
        end
        % Start new component
        currentComponent = currentComponent + 1;
        queue = i;
        componentId(i) = currentComponent;
    
        while ~isempty(queue)
            f = queue(end);
            queue(end) = [];
    
            neighbors = faceAdj(f, :);
            neighbors = neighbors(~isnan(neighbors));  % remove NaNs
            for n = neighbors
                if componentId(n) == 0
                    componentId(n) = currentComponent;
                    queue(end+1) = n; %#ok<AGROW>
                end
            end
        end
    end
    
    % Step 5: Select the largest and smallest components
    counts = histcounts(componentId, 1:max(componentId)+1);
    [~, largestComponent] = max(counts);
    [~, smallestComponent] = min(counts);
    keepFaces1 = find(componentId == largestComponent);
    F_final1 = F_clean(keepFaces1, :);
    keepFaces2 = find(componentId == smallestComponent);
    F_final2 = F_clean(keepFaces2, :);
    
    % Step 6: Reindex vertices
    [uniqueVerts21, ~, ic21] = unique(F_final1(:));
    V_final1 = V_clean(uniqueVerts21, :);
    F_final1 = reshape(ic21, size(F_final1));
    [uniqueVerts22, ~, ic22] = unique(F_final2(:));
    V_final2 = V_clean(uniqueVerts22, :);
    F_final2 = reshape(ic22, size(F_final2));

    %[OuterDist,~,in2outNorms] = point2trimesh('Faces', F_final2, 'Vertices', V_final2, 'QueryPoints', V_final1);
    
    %{
    % Step 7: Visualize manifold 1 (outer)
    figure;
    trisurf(F_final1, V_final1(:,1), V_final1(:,2), V_final1(:,3), ...
        'FaceAlpha', 0.6);
    title('manifold mesh 1');
    axis equal;
    
    % Step 8: Visualize manifold 2 (inner)
    figure;
    trisurf(F_final2, V_final2(:,1), V_final2(:,2), V_final2(:,3), ...
        'FaceAlpha', 0.6);
    title('manifold mesh 2');
    axis equal;
    %}

    V_final2nnIdx = knnsearch(V_final2, V_final1);
    V_final2nn = V_final2(V_final2nnIdx,:); % Nearest nighbors in 2 for each 1
    
    faceCtroid1 = (V_final1(F_final1(:,1),:) + V_final1(F_final1(:,2),:) + V_final1(F_final1(:,3),:)) / 3;
    nearestInside1 = knnsearch(poleColumn, faceCtroid1);
    insideVec1 = faceCtroid1 - poleColumn(nearestInside1,:);
    v1cross1 = V_final1(F_final1(:,2),:) - V_final1(F_final1(:,1),:);
    v2cross1 = V_final1(F_final1(:,3),:) - V_final1(F_final1(:,1),:);
    Nvec1 = cross(v1cross1, v2cross1, 2);
    Nvecnorm1 = Nvec1 ./ vecnorm(Nvec1, 2, 2);
    dotProd1 = sum(insideVec1 .* Nvecnorm1, 2);
    negNorm1 = dotProd1 < 0;
    Nvecnorm1(negNorm1,:) = -1*Nvecnorm1(negNorm1,:);
    nodeNormals1 = zeros(size(V_final1,1),3);
    for i = 1:3
        idx = F_final1(:, i);
        nodeNormals1(:,1) = nodeNormals1(:,1) + accumarray(idx, Nvecnorm1(:,1), [size(V_final1,1),1], @sum, 0);
        nodeNormals1(:,2) = nodeNormals1(:,2) + accumarray(idx, Nvecnorm1(:,2), [size(V_final1,1),1], @sum, 0);
        nodeNormals1(:,3) = nodeNormals1(:,3) + accumarray(idx, Nvecnorm1(:,3), [size(V_final1,1),1], @sum, 0);
    end
    nodeNormals1 = nodeNormals1 ./ vecnorm(nodeNormals1, 2, 2);
    
    thicknessCorr1 = 0.0001 - abs(sum(nodeNormals1 .* (V_final1 - V_final2nn), 2));
    thicknessCorr1(thicknessCorr1 < 0.000005) = 0;
    % new part
    V_final1Upper = V_final1(:,3) > 0.003;
    V_final1new = V_final1;
    V_final1new(V_final1Upper,:) = V_final1(V_final1Upper,:) + (0.5 * thicknessCorr1(V_final1Upper)) .* nodeNormals1(V_final1Upper,:);
    %V_final1new = V_final1 + (0.5 * thicknessCorr1) .* nodeNormals1;  % Outward extrusion
    [~, idx_map] = ismember(V_final1, V, 'rows');
    V(idx_map,:) = V_final1new;
    
    
    V_final1nnIdx = knnsearch(V_final1new, V_final2);
    V_final1nn = V_final1(V_final1nnIdx,:); % Nearest nighbors in 1 for each 2
    
    faceCtroid2 = (V_final2(F_final2(:,1),:) + V_final2(F_final2(:,2),:) + V_final2(F_final2(:,3),:)) / 3;
    nearestInside2 = knnsearch(poleColumn, faceCtroid2);
    insideVec2 = faceCtroid2 - poleColumn(nearestInside2,:);
    v1cross2 = V_final2(F_final2(:,2),:) - V_final2(F_final2(:,1),:);
    v2cross2 = V_final2(F_final2(:,3),:) - V_final2(F_final2(:,1),:);
    Nvec2 = cross(v1cross2, v2cross2, 2);
    Nvecnorm2 = Nvec2 ./ vecnorm(Nvec2, 2, 2);
    dotProd2 = sum(insideVec2 .* Nvecnorm2, 2);
    negNorm2 = dotProd2 > 0;
    F_final2(negNorm2,:) = F_final2(negNorm2, [1 3 2]);
    Nvecnorm2(negNorm2,:) = -1*Nvecnorm2(negNorm2,:);
    nodeNormals2 = zeros(size(V_final2,1),3);
    for i = 1:3
        idx = F_final2(:, i);
        nodeNormals2(:,1) = nodeNormals2(:,1) + accumarray(idx, Nvecnorm2(:,1), [size(V_final2,1),1], @sum, 0);
        nodeNormals2(:,2) = nodeNormals2(:,2) + accumarray(idx, Nvecnorm2(:,2), [size(V_final2,1),1], @sum, 0);
        nodeNormals2(:,3) = nodeNormals2(:,3) + accumarray(idx, Nvecnorm2(:,3), [size(V_final2,1),1], @sum, 0);
    end
    nodeNormals2 = nodeNormals2 ./ vecnorm(nodeNormals2, 2, 2);
    
    thicknessCorr2 = 0.0001 - abs(sum(nodeNormals2 .* (V_final2 - V_final1nn), 2));
    thicknessCorr2(thicknessCorr2 < 0.000005) = 0;
    % new part
    V_final2Upper = V_final2(:,3) > 0.003;
    V_final2new = V_final2;
    V_final2new(V_final2Upper,:) = V_final2(V_final2Upper,:) + (0.5 * thicknessCorr2(V_final2Upper)) .* nodeNormals2(V_final2Upper,:);
    %V_final1new = V_final1 + (0.5 * thicknessCorr1) .* nodeNormals1;  % Outward extrusion
    [~, idx_map] = ismember(V_final2, V, 'rows');
    V(idx_map,:) = V_final2new;
    
    volp1 = V_final2new(F_final2(:,1), :);
    volp2 = V_final2new(F_final2(:,2), :);
    volp3 = V_final2new(F_final2(:,3), :);
    signed_volumes = dot(volp1, cross(volp2, volp3, 2), 2);
    volume = abs(sum(signed_volumes) / 6.0);
    %{
    % Initialize femodel with single domain geometry
    model2 = femodel("AnalysisType","structuralStatic");
    model2.Geometry = fegeometry(V, T);
    
    % Plot geometry
    figure;
    pdegplot(model2.Geometry, FaceLabels='on', FaceAlpha=0.5);
    axis equal;
    title('initial single-domain geometry')
    %}
end
%{
% Initialize femodel with single domain geometry
model = femodel("AnalysisType","structuralStatic","Geometry", "HollowHemisphere.step");
g = model.Geometry;

% Plot geometry
figure;
pdegplot(g, FaceLabels='on', FaceAlpha=0.5);
axis equal;
title('initial single-domain geometry')

% Generate tetrahedral mesh, can specify resolution
model = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.00025, 'GeometricOrder', 'linear');

% Plot mesh
figure;
pdemesh(model,'FaceAlpha', 0.3);
title('Tetrahedral Mesh of Direct hollow sphere');

% Define pole column for knnsearch
N = 8;  % Choose however many points you want
p1 = [0, 0, 0.002];
p2 = [0, 0, 0.008];
t = linspace(0, 1, N)';   % N values from 0 to 1 (column vector)
poleColumn = (1 - t) * p1 + t * p2; 
%}
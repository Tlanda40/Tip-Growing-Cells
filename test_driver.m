% Helper function to define extensibility profile
function y = ext(beta, x, y)
    C1=abs(beta(1));
    C2=abs(beta(2));
    a1=beta(3);
    a2=beta(4);
    y=C1*exp(-((sqrt(x.^2+y.^2))/a1).^2)+C2*exp(-((sqrt(x.^2+y.^2))/a2).^2);
end

% Helper function to define stiffness
function y = youngMod(Ext)
    y = 500 ./ Ext;
end

% Initialize femodel with single domain geometry
model = femodel("AnalysisType","structuralStatic","Geometry", "HollowHemisphere.step");
g = model.Geometry;

% Plot geometry
figure;
pdegplot(g, FaceLabels='on', FaceAlpha=0.5);
axis equal;
title('initial single-domain geometry')

% Generate tetrahedral mesh, can specify resolution
model = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.0002, 'GeometricOrder', 'linear');
size(model.Mesh.Nodes)

% Plot mesh
figure;
pdemesh(model,'FaceAlpha', 0.3);
title('Tetrahedral Mesh of Direct hollow sphere');

% Step 1: Get free boundary
V = model.Mesh.Nodes';
T = model.Mesh.Elements';
TR = triangulation(T, V);
[F, V] = freeBoundary(TR);

% Step 2: Clean up vertex indexing
[uniqueVerts, ~, ic] = unique(F(:));
V_clean = V(uniqueVerts, :);
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

% Step 5: Select the largest component (assumed outer surface)
counts = histcounts(componentId, 1:max(componentId)+1);
[~, largestComponent] = max(counts);
keepFaces = find(componentId == largestComponent);
F_final = F_clean(keepFaces, :);

% Step 6: Reindex vertices
[uniqueVerts2, ~, ic2] = unique(F_final(:));
V_final = V_clean(uniqueVerts2, :);
F_final = reshape(ic2, size(F_final));

figure;
trisurf(F_final, V_final(:,1), V_final(:,2), V_final(:,3), ...
    'FaceAlpha', 0.6);

% Find index of pole vertex
poleCoord = [0, 0, 0.01];           % Your pole location
poleIdx = knnsearch(V_final, poleCoord);  % Closest vertex on the surface
fprintf("Pole index is %d, pole coords are [%d, %d, %d]\n", poleIdx, V_final(poleIdx, 1), V_final(poleIdx, 2), V_final(poleIdx, 3));

targetCoord = [-0.00545175, -0.00450885, 0.007044]; % a target point
targetIdx = knnsearch(V_final, targetCoord); % the closest node to target
targetIdx

%write_off('init_mesh_test.off', V_final, F_final);

geodistances = readmatrix('distances.txt');  % Nx2 matrix
vertex_indices = geodistances(:,1) + 1;      % CGAL is 0-indexed, MATLAB is 1-indexed
D = geodistances(:,2);

% Visualize arclength by color
figure;
trisurf(F_final, V_final(:,1), V_final(:,2), V_final(:,3), D, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9);
colormap(turbo); colorbar;
title('Geodesic Distance from Pole');
hold on;
axis equal;

P = readmatrix("path.txt");
plot3(P(:,1), P(:,2), P(:,3), 'r-', 'LineWidth', 2);
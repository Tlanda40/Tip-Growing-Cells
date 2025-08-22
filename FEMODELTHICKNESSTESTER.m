function vols = tetVolume(nodes, elements)
% nodes: 3×N array of node coordinates (x; y; z)
% elements: 4×M array of tetrahedral elements (indices into nodes)

    p1 = nodes(:, elements(1, :));
    p2 = nodes(:, elements(2, :));
    p3 = nodes(:, elements(3, :));
    p4 = nodes(:, elements(4, :));

    % Compute signed volume of each tetrahedron
    v = dot(cross(p2 - p1, p3 - p1, 1), p4 - p1, 1) / 6;

    vols = v;  % 1×M array of signed volumes
end
%{
function writeSTL(filename, faces, vertices, faceNorms)
    % Simple ASCII STL exporter
    fid = fopen(filename, 'w');
    fprintf(fid, 'solid matlab_export\n');
    
    for i = 1:size(faces, 1)
        v1 = vertices(faces(i, 1), :);
        v2 = vertices(faces(i, 2), :);
        v3 = vertices(faces(i, 3), :);
        
        % Access face normal
        n = faceNorms(i,:);
        
        fprintf(fid, '  facet normal %.6f %.6f %.6f\n', n);
        fprintf(fid, '    outer loop\n');
        fprintf(fid, '      vertex %.6f %.6f %.6f\n', v1);
        fprintf(fid, '      vertex %.6f %.6f %.6f\n', v2);
        fprintf(fid, '      vertex %.6f %.6f %.6f\n', v3);
        fprintf(fid, '    endloop\n');
        fprintf(fid, '  endfacet\n');
    end
    
    fprintf(fid, 'endsolid matlab_export\n');
    fclose(fid);
    
    fprintf('STL written to %s\n', filename);
end
%}

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

% TOPOLOGICAL MANIFOLD EXTRACTION PROCESS FOR CGAL GEODESIC
% Step 1: Get free boundary
V = model.Mesh.Nodes';
T = model.Mesh.Elements';
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
V_final1new = V_final1 + (0.5 * thicknessCorr1) .* nodeNormals1;  % Outward extrusion
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
V_final2new = V_final2 + (thicknessCorr2) .* nodeNormals2;  % inward extrusion
[~, idx_map] = ismember(V_final2, V, 'rows');
V(idx_map,:) = V_final2new;

model2 = femodel("AnalysisType","structuralStatic");
model2.Geometry = fegeometry(V, T);
figure;
pdegplot(model2.Geometry, FaceLabels='on', FaceAlpha=0.5);
axis equal;
title('initial single-domain geometry')



%{
% Step 9: Calculate average vertices
V_final2nnIdx = knnsearch(V_final2, V_final1);
V_final2nn = V_final2(V_final2nnIdx,:);
avgXCoords = mean([V_final1(:,1),V_final2nn(:,1)],2);
avgYCoords = mean([V_final1(:,2),V_final2nn(:,2)],2);
avgZCoords = mean([V_final1(:,3),V_final2nn(:,3)],2);
V_avg = [avgXCoords, avgYCoords, avgZCoords];


% Correct face orientations and extract node normals
faceCtroid = (V_avg(F_final1(:,1),:) + V_avg(F_final1(:,2),:) + V_avg(F_final1(:,3),:)) / 3;
nearestInside = knnsearch(poleColumn, faceCtroid);
insideVec = faceCtroid - poleColumn(nearestInside,:);
v1cross = V_avg(F_final1(:,2),:) - V_avg(F_final1(:,1),:);
v2cross = V_avg(F_final1(:,3),:) - V_avg(F_final1(:,1),:);
Nvec = cross(v1cross, v2cross, 2);
Nvecnorm = Nvec ./ vecnorm(Nvec, 2, 2);
dotProd = sum(insideVec .* Nvecnorm, 2);
negNorm = dotProd < 0;
F_fixed = F_final1;
Nvecnorm(negNorm,:) = -1*Nvecnorm(negNorm,:);
F_fixed(negNorm,:) = F_final1(negNorm, [1 3 2]);
nodeNormals = zeros(size(V_avg,1),3);
for i = 1:3
    idx = F_fixed(:, i);
    nodeNormals = nodeNormals + accumarray(idx, Nvecnorm(:,1), [size(V_avg,1),1], @sum, 0);
    nodeNormals(:,2) = nodeNormals(:,2) + accumarray(idx, Nvecnorm(:,2), [size(V_avg,1),1], @sum, 0);
    nodeNormals(:,3) = nodeNormals(:,3) + accumarray(idx, Nvecnorm(:,3), [size(V_avg,1),1], @sum, 0);
end
nodeNormals = nodeNormals ./ vecnorm(nodeNormals, 2, 2);

% Step 8: Visualize average manifold
figure;
trisurf(F_fixed, avgXCoords, avgYCoords, avgZCoords, ...
    'FaceAlpha', 0.6);
title('manifold mesh avg');
axis equal;

thickness = 0.0001;  % Total shell thickness
V_inner = V_avg - (0.5*thickness) * nodeNormals;  % Inward extrusion
V_outer = V_avg + (0.5*thickness) * nodeNormals;  % Outward extrusion
%}

%{
F2_shifted = fliplr(F_fixed) + size(V_outer,1); % shift indices of second mesh
V_combined = [V_outer; V_inner];
F_combined = [F_fixed; F2_shifted];
TRdouble = triangulation(F_combined, V_combined);
stlwrite(TRdouble, 'doubleManifold.stl');


TRinner = triangulation(fliplr(F_fixed), V_inner);
TRouter = triangulation(F_fixed, V_outer);
stlwrite(TRinner, 'manifoldInn.stl');
stlwrite(TRinner, 'manifoldOut.stl');
%}
%{
figure;
trisurf(F_fixed, V_inner(:,1), V_inner(:,2), V_inner(:,3), ...
    'FaceAlpha', 1, 'FaceColor', 'red');
title('manifold mesh 2');
axis equal;
hold on;
trisurf(F_fixed, V_outer(:,1), V_outer(:,2), V_outer(:,3), ...
    'FaceAlpha', 0.6);
title('manifold mesh 2');
axis equal;
%}




%{
edata = load('elements.mat');
pelements = edata.elements;
ndata = load('nodes.mat');
pnodes = ndata.nodes;

model2 = femodel("AnalysisType","structuralStatic");
model2.Geometry = fegeometry(pnodes, pelements);

% Plot geometry
figure;
pdegplot(model2.Geometry, FaceLabels='on', FaceAlpha=0.5);
axis equal;
title('initial single-domain geometry')
%}






%{
figure;
tetramesh(E_new, V_new);
axis equal; view(3);
%}
%{
V_new = [V_outer; V_inner];
F_new = [F_fixed, F_fixed + size(V_outer, 1)];
E_new = zeros(3*size(F_fixed, 1), 4);
for i = 1:size(F_avg, 1)
    E_new(3*i-2, :) = [F_new(i, 2), F_new(i, 4), F_new(i, 5), F_new(i, 6)];
    E_new(3*i-1, :) = [F_new(i, 1), F_new(i, 2), F_new(i, 4), F_new(i, 6)];
    E_new(3*i, :) = [F_new(i, 1), F_new(i, 2), F_new(i, 3), F_new(i, 6)];
end

ENewVols = tetVolume(V_new', E_new')';
bad = ENewVols < 0;
E_new(bad, :) = E_new(bad, [1 2 4 3]);
%}

%model2 = femodel("AnalysisType","structuralStatic");
%model2.Geometry = fegeometry(V_new, E_new);

%{
figure;
tetramesh(E_new, V_new);
axis equal; view(3);
%}

%{
%{
figure;
trisurf(F_avg, V_inner(:,1), V_inner(:,2), V_inner(:,3), ...
    'FaceAlpha', 1, 'FaceColor', 'red');
title('manifold mesh 2');
axis equal;
hold on;
trisurf(F_avg, V_outer(:,1), V_outer(:,2), V_outer(:,3), ...
    'FaceAlpha', 0.6);
title('manifold mesh 2');
axis equal;
%}

%{
nodes_combined = [V_outer; V_inner];
faces_inner_shifted = fliplr(F_final1) + size(V_outer, 1);
faces_combined = [F_final1; faces_inner_shifted];


maxh = 7.5e-12;
[mesh_nodes, mesh_tets] = surf2mesh( ...
    nodes_combined, ...
    faces_combined, ...
    [], [], ...         % Bounding box and holes (optional here)
    0, ...
    struct('maxvol', 1e-3, 'method', 'tetgen'));

figure;
tetramesh(mesh_tets, mesh_nodes);
axis equal; view(3);
%}
%{
V_new = [V_outer; V_inner];
F_new = [F_avg, F_avg + size(V_inner, 1)];
E_new = zeros(3*size(F_avg, 1), 4);
for i = 1:size(F_avg, 1)
    E_new(3*i-2, :) = [F_new(i, 2), F_new(i, 4), F_new(i, 5), F_new(i, 6)];
    E_new(3*i-1, :) = [F_new(i, 1), F_new(i, 2), F_new(i, 4), F_new(i, 6)];
    E_new(3*i, :) = [F_new(i, 1), F_new(i, 2), F_new(i, 3), F_new(i, 6)];
end

ENewVols = tetVolume(V_new', E_new')';
bad = ENewVols < 0;
E_new(bad, :) = E_new(bad, [1 2 4 3]);
%}

V_all = [V_inner; V_outer];
num_nodes = size(V_inner, 1);
num_faces = size(F_avg, 1);
% Loop over each triangle
tets = zeros(numFaces*3, 4);  % 3 tets per face
tetCount = 1;

for i = 1:num_faces
    % Get local node indices for current triangle
    inner_nodes = F_avg(i, :);                 % indices into V_inner
    outer_nodes = F_avg(i, :) + num_nodes;     % shift indices for V_outer

    % Combine into prism: [i1 i2 i3 | o1 o2 o3]
    % So: p = [1 2 3 4 5 6] = [inner1 inner2 inner3 outer1 outer2 outer3]
    prism_nodes = [inner_nodes, outer_nodes];

    % Define tetrahedra using node indices (not coordinates!)
    local_tets = [
        prism_nodes([1 2 3 6]);
        prism_nodes([1 3 5 6]);
        prism_nodes([1 4 5 6]);
    ];

    tets(tetCount:tetCount+2, :) = local_tets;
    tetCount = tetCount + 3;
end

% 'tets' now contains your shell volume mesh

figure;
tetramesh(tets, V_all);
%TR = triangulation(E_new, V_new);


model2 = femodel("AnalysisType","structuralStatic");
model2.Geometry = fegeometry(V_all, tets);
%}
%{
function checkTetMeshConnectivity(nodes, tets)
    %CHECKTETMESHCONNECTIVITY Performs basic connectivity checks on a tetrahedral mesh
    %
    % Inputs:
    %   nodes - Nx3 array of vertex positions
    %   tets  - Mx4 array of tetrahedra (vertex indices into nodes)
    
    fprintf('--- Mesh Connectivity Report ---\n');
    fprintf('Nodes: %d\n', size(nodes,1));
    fprintf('Tetrahedra: %d\n', size(tets,1));
    
    % --- Step 1: Extract all faces
    F = [tets(:,[1 2 3]);
         tets(:,[1 2 4]);
         tets(:,[1 3 4]);
         tets(:,[2 3 4])];
     
    % --- Step 2: Sort for uniqueness (regardless of winding)
    F_sorted = sort(F, 2);
    [unique_faces, ~, ic] = unique(F_sorted, 'rows');
    face_counts = accumarray(ic, 1);
    
    % --- Step 3: Boundary and Non-manifold check
    num_boundary_faces = sum(face_counts == 1);
    num_nonmanifold_faces = sum(face_counts > 2);
    
    fprintf('Boundary faces: %d\n', num_boundary_faces);
    if num_nonmanifold_faces > 0
        fprintf('❌ Non-manifold faces: %d\n', num_nonmanifold_faces);
    else
        fprintf('✅ No non-manifold faces found.\n');
    end
    
    % --- Step 4: Check component connectivity using graph theory
    E = [tets(:,[1 2]); tets(:,[1 3]); tets(:,[1 4]);
         tets(:,[2 3]); tets(:,[2 4]); tets(:,[3 4])];
    E = unique(sort(E,2), 'rows');
    G = graph(E(:,1), E(:,2));
    component_labels = conncomp(G);
    num_components = max(component_labels);
    
    fprintf('Connected components: %d\n', num_components);
    if num_components > 1
        fprintf('❌ Mesh is disconnected.\n');
    else
        fprintf('✅ Mesh is fully connected.\n');
    end
    
    % --- Optional: Show histogram of node degree
    node_degree = degree(G);
    fprintf('Node degree range: [%d - %d]\n', min(node_degree), max(node_degree));
    
    % --- Optional: Return outputs?
    % Could add outputs here if needed
    
    fprintf('--- End of Report ---\n');
end

checkTetMeshConnectivity(V_new, E_new);
%}
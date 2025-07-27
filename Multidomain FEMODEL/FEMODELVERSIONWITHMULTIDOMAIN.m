%{
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
%}

% Helper function to define extensibility profile
function y = ext(beta, x)
    C1=abs(beta(1));
    C2=abs(beta(2));
    a1=beta(3);
    a2=beta(4);
    y=C1*exp(-(x/a1).^2)+C2*exp(-(x/a2).^2);
end

% Helper function to define stiffness
function y = youngMod(Ext)
    y = 500 ./ Ext;
end

% import single domain geometry
gmSingleDomain = fegeometry("HollowHemisphere.step");
gmSingleDomain = generateMesh(gmSingleDomain);

% mesh single domain for subdivision
%figure;
%pdemesh(gmSingleDomain);
msh = gmSingleDomain.Mesh;
elements = msh.Elements;

% Create matrix of centroids of elements, calculate xy distance from 0
centroids = zeros(3,size(elements,2));
xyDist = zeros(1,size(elements,2));
for i = 1:size(elements,2)
    centroids(1,i) = mean(msh.Nodes(1,elements(1:4,i)));
    centroids(2,i) = mean(msh.Nodes(2,elements(1:4,i)));
    centroids(3,i) = mean(msh.Nodes(3,elements(1:4,i)));
    xyDist(1,i) = sqrt(centroids(1,i)^2 + centroids(2,i)^2);
end

% create multi domain geometry
ElementIdToRegionId = 1:size(elements,2);
gmMultiDomain = fegeometry(msh.Nodes',elements',ElementIdToRegionId);
%figure;
%pdegplot(gmMultiDomain,FaceAlpha=0.5)

% Initialize femodel with multi domain geoemtry
model = femodel("AnalysisType","structuralStatic","Geometry", gmMultiDomain);
g = model.Geometry;

% Plot geometry
figure;
pdegplot(g, FaceAlpha=0.5);
axis equal;
title('initial multi-domain geometry')

% Find boundary faces (faces belonging to only one cell)
faceUseCount = zeros(g.NumFaces, 1);  % face usage counter

% Step 1: Loop over all cells and count face appearances
for cellIdx = 1:g.NumCells
    theseFaces = g.cellFaces(cellIdx);  % get faces of this cell
    faceUseCount(theseFaces) = faceUseCount(theseFaces) + 1;
end

% Step 2: Boundary faces are those used in only one cell
boundaryFaces = find(faceUseCount == 1);

fprintf('Found %d boundary faces.\n', numel(boundaryFaces));

if exist('components_cache.mat', 'file')
    load('components_cache.mat', 'components');
else
    % Fallback: compute it again if file is missing
    % Step 1: Build face adjacency based on shared edges
    % Initialize adjacency matrix for boundary faces only
    n = numel(boundaryFaces);
    adj = sparse(n, n);  % sparse saves memory for large data
    
    for i = 1:n
        fi = boundaryFaces(i);
        edges_i = g.faceEdges(fi);
        
        for j = i+1:n
            fj = boundaryFaces(j);
            edges_j = g.faceEdges(fj);
        
            % Check if the two faces share an edge
            if ~isempty(intersect(edges_i, edges_j))
                adj(i, j) = 1;
                adj(j, i) = 1;
            end
        end
        fprintf("  At boundary face %d of %d\n", i, n);
    end

    % Step 2: BFS to find connected faces starting from an initial boundary face
    visited = false(n, 1);       % Track visited boundary face indices
    components = {};             % Each entry holds one group of boundary face indices
    
    for i = 1:n
        if ~visited(i)
            queue = i;
            visited(i) = true;
            component = i;
    
            while ~isempty(queue)
                current = queue(1);
                queue(1) = [];  % Dequeue

                neighbors = find(adj(current, :));
                for nb = neighbors
                    if ~visited(nb)
                        visited(nb) = true;
                        queue(end+1) = nb;
                        component(end+1) = nb;
                    end
                end
            end

            components{end+1} = boundaryFaces(component);  % map back to global face indices
        end
    end
    save('components_cache.mat', 'components');
end

% Determine which component stores inner/outer faces by size
if size(components{1}) > size(components{2})
    outerFaces = components{1};
    innerFaces = components{2};
else
    outerFaces = components{2};
    innerFaces = components{1};
end

fprintf('Found %d outer faces.\n', numel(outerFaces));
fprintf('Found %d inner faces.\n', numel(innerFaces));

% Generate tetrahedral mesh, can specify resolution
mesh = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.0002, 'GeometricOrder', 'linear');

% Plot mesh
figure;
pdemesh(model,'FaceAlpha', 0.3);
title('Tetrahedral Mesh of Direct hollow sphere');

%insert structural properties. For now, it will be some kind of basic
%thing, with uniform structural properties (as opposed to properties that
%are a function of arc length.

%E = 1e4;            % Base Young's modulus (Pa)
nu = 0.3;           % Poisson's ratio
rho = 1000;         % Mass Density, only relevant in a transient-solid model
p = 1e2;            % Pressure
beta0 = [0.00375,0.13875,0.01,0.001]'; % Inputs for extensibility profile

%{
%Find initial volume enclosed (volume of void inside shell)
face3Nodes = findNodes(mesh, 'region', 'Face', 3);
face3NodeCoords = mesh.Nodes(:, face3Nodes)';
face4Nodes = findNodes(mesh, 'region', 'Face', 4);
face4NodeCoords = mesh.Nodes(:, face4Nodes)';
InnerFaceNodes_combined = [face3NodeCoords; face4NodeCoords];
shp = alphaShape(InnerFaceNodes_combined(:,1), InnerFaceNodes_combined(:,2), InnerFaceNodes_combined(:,3));
[faces, points] = boundaryFacets(shp);
shp.Alpha = 0.01;
figure;
h = plot(shp)
h.FaceAlpha = 0.3;
title('Inner void')
axis equal;

Vi = volume(shp); %this is the initial volume enclosed
%V = mesh.volume() (this is the volume of the shell itself)

%Here, we normalize pressure with respect to volume
PV_ratio = p / Vi;
% And later, we apply structuralBoundaryLoad(model, 'Face', [3,4], 'Pressure', PV_ratio * Vi);
%}

% This is where the FE expansion is set up.
% Here are the material properties
if exist('MaterialPropertiesExt4_cache.mat', 'file')
    load('MaterialPropertiesExt4_cache.mat', 'matProps');
    model.MaterialProperties = matProps;
else
    % Fallback: compute it again if file is missing
    for i = 1:g.NumCells
        model.MaterialProperties(i).PoissonsRatio = nu;
        model.MaterialProperties(i).YoungsModulus = ...
            youngMod(ext(beta0, xyDist(i)));
        fprintf('Applied properties to %d of %d cells.\n', i, g.NumCells);
    end
    matProps = model.MaterialProperties;
    save('MaterialPropertiesExt4_cache.mat', 'matProps');
end

% Applying Neumann boundary conditions to innerFaces
model.FaceLoad(innerFaces) = faceLoad(Pressure = p); % Pressure in Pascals

% Finding the bottom vertices, and applying Dirichlet boundary conditions
z_thresh = 0.00105;
indices = find(g.Vertices(:,3) <= z_thresh);
model.VertexBC(indices) = vertexBC(Constraint="fixed");

fprintf('BCs and Material Properties set.');

% solve
if exist('resultsExt4_cache.mat', 'file')
    load('resultsExt4_cache.mat', 'results');
else
    % Fallback: compute it again if file is missing
    results = solve(model);
    save('resultsExt4_cache.mat', 'results');
end

fprintf('Solved!');

% visualize
u = results.Displacement;
figure;
pdeplot3D(model.Mesh,'ColorMapData', u.Magnitude, 'FaceAlpha', 0.3, ...
    'Deformation', u, ...
    'DeformationScaleFactor', 1);
title('Displacement colormap')

%{
result = solve(model);
u = result.Displacement;
umag = sqrt(u.ux.^2 + u.uy.^2 + u.uz.^2);
figure;
pdeplot3D(model,'ColorMapData', umag, 'FaceAlpha', 0.3, ...
    'Deformation', u, ...
    'DeformationScaleFactor', 1);
title('Displacement colormap')

displacedNodes = mesh.Nodes + [u.ux'; u.uy'; u.uz'];
elements = mesh.Elements;
elementsLinear = elements(1:4, :);
vols = tetVolume(displacedNodes, elementsLinear);
bad = vols < 0;
% Flip nodes 3 and 4 for inverted tets (common fix)
elementsLinear(:, bad) = elementsLinear([1 2 4 3], bad);

%Now, we want to loop with normalized pressure

for i = 1:8
    newModel = createpde('structural','static-solid');
    geometryFromMesh(newModel, displacedNodes, elementsLinear);
    mesh = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.0002, 'GeometricOrder', 'linear')
    face3Nodes = findNodes(mesh, 'region', 'Face', 3);
    face3NodeCoords = mesh.Nodes(:, face3Nodes)';
    face4Nodes = findNodes(mesh, 'region', 'Face', 4);
    face4NodeCoords = mesh.Nodes(:, face4Nodes)';
    InnerFaceNodes_combined = [face3NodeCoords; face4NodeCoords];
    shp = alphaShape(InnerFaceNodes_combined(:,1), InnerFaceNodes_combined(:,2), InnerFaceNodes_combined(:,3));
    [faces, points] = boundaryFacets(shp);
    shp.Alpha = 0.01;
    Vn = volume(shp);
    structuralProperties(newModel, 'YoungsModulus', E, 'PoissonsRatio', nu);
    structuralBoundaryLoad(newModel, 'Face', [3,4], 'Pressure', PV_ratio * Vn);
    structuralBC(newModel,"Face",[2,4],"Constraint","fixed");
    result = solve(newModel);
    u = result.Displacement;
    umag = sqrt(u.ux.^2 + u.uy.^2 + u.uz.^2);
    figure;
    pdeplot3D(newModel,'ColorMapData', umag, 'FaceAlpha', 0.3, ...
        'Deformation', u, ...
        'DeformationScaleFactor', 1);
    title('Displacement colormap')
    displacedNodes = mesh.Nodes + [u.ux'; u.uy'; u.uz'];
    elements = mesh.Elements;
    elementsLinear = elements(1:4, :);
    vols = tetVolume(displacedNodes, elementsLinear);
    bad = vols < 0;
    % Flip nodes 3 and 4 for inverted tets (common fix)
    elementsLinear(:, bad) = elementsLinear([1 2 4 3], bad);
end

%}
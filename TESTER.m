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

% Extract Nx3 vertices matrix
V = model.Mesh.Nodes';  %transpose because nodes is 3xN

% Extract Mx3 faces matrix
T = model.Mesh.Elements';  %This is a Mx4 tetrahedra matrix
function F = getSurfaceFacesFromTetrahedra(T)  %This takes T to F
    % Extract faces from Mx4 tetrahedra matrix
    f = [T(:, [1 2 3]);
             T(:, [1 2 4]);
             T(:, [1 3 4]);
             T(:, [2 3 4])];
    f = sort(f, 2); % Sorts each row (face)
    f = unique(f, 'rows'); % Remove duplicates
    F = f;
end
F = getSurfaceFacesFromTetrahedra(T);

% Find index of pole vertex
poleCoord = [0, 0, 0.01];           % Your pole location
poleIdx = knnsearch(V, poleCoord);  % Closest vertex on the surface
fprintf("Pole index is %d, pole coords are [%d, %d, %d]\n", poleIdx, V(poleIdx, 1), V(poleIdx, 2), V(poleIdx, 3));


%THIS MAKES AN ALPHASHAPE
%{
% Make alphashape
shp = alphaShape(model.Mesh.Nodes', 1);  % choose alpha carefully (this is too coarse)
figure;
plot(shp);
%}

%{
%THIS IS THE WORKING METHOD USING DIJKSTRA'S ALGORITHM

% Create V (vertices) and T (tetrahedra) for graph creation
V = model.Mesh.Nodes';
T = model.Mesh.Elements';

% Mesh to Graph
function G = meshToSparseGraph(V, T)
    % Extract edges from Mx4 tetrahedra matrix (each row has 4 indices of V)
    E = [T(:, [1, 2]); T(:, [2, 3]); T(:, [3, 4]); T(:, [4, 1]); T(:, [1, 3]); T(:, [2, 4])];  % 6Mx2, each row is pair of vertex indices
    E = sort(E, 2);  % Sorts each row (edge)
    E = unique(E, 'rows');  % Remove duplicates

    % Compute edge lengths
    edgeLengths = vecnorm(V(E(:,1),:) - V(E(:,2),:), 2, 2);

    % Create sparse adjacency graph
    G = graph(E(:,1), E(:,2), edgeLengths, size(V,1));
end

% Find index of pole vertex
poleCoord = [0, 0, 0.01];           % Your pole location
poleIdx = knnsearch(V, poleCoord);  % Closest vertex on the surface
fprintf("Pole index is %d, pole coords are [%d, %d, %d]\n", poleIdx, V(poleIdx, 1), V(poleIdx, 2), V(poleIdx, 3));

% Compute lengths of geodesics. s stores arclength.
G = meshToSparseGraph(V, T);
s = distances(G, poleIdx);  % s(i) is geodesic distance to vertex i

% Extract faces from tetrahedral elements to visualize s with trisurf
function F = getSurfaceFacesFromTetrahedra(T)
    % Extract faces from Mx4 tetrahedra matrix
    f = [T(:, [1 2 3]);
             T(:, [1 2 4]);
             T(:, [1 3 4]);
             T(:, [2 3 4])];
    f = sort(f, 2); % Sorts each row (face)
    f = unique(f, 'rows'); % Remove duplicates
    F = f;
end

F = getSurfaceFacesFromTetrahedra(T);

% Visualize arclength by color
figure;
trisurf(F, V(:,1), V(:,2), V(:,3), s, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.4);
colormap(turbo); colorbar;
title('Geodesic Distance from Pole');
hold on;
axis equal;

% visualize one geodesic, to check
targetCoord = [-0.00545175, -0.00450885, 0.007044]; % a target point
targetIdx = knnsearch(V, targetCoord); % the closest node to target
pathIdx = shortestpath(G, poleIdx, targetIdx); % a vector of points along the geodesic
pathCoords = V(pathIdx, :);  % Nx3, their coordinates
plot3(pathCoords(:,1), pathCoords(:,2), pathCoords(:,3), ...
      'r-', 'LineWidth', 3);

%}


%insert structural properties. For now, it will be some kind of basic
%thing, with uniform structural properties (as opposed to properties that
%are a function of arc length.

E = 1e6;            % Base Young's modulus (Pa)
nu = 0.3;           % Poisson's ratio
rho = 1000;         % Mass Density, only relevant in a transient-solid model
p = 1e2;            % Pressure
beta0 = [0.00375,0.13875,0.01,0.001]'; % Inputs for extensibility profile

% This is where the FE expansion is set up.

% Here are the material properties
model.MaterialProperties.PoissonsRatio = nu;

val = @(location,state) youngMod(ext(beta0, location.x, location.y));
model.MaterialProperties.YoungsModulus = val;

% Applying Neumann boundary conditions to innerFaces
model.FaceLoad(3) = faceLoad(Pressure = p); % Pressure in Pascals

% Finding the bottom vertices, and applying Dirichlet boundary conditions
model.FaceBC([2,4]) = faceBC(Constraint="fixed");
% Finding the bottom vertices, and applying Dirichlet boundary conditions
%z_thresh = 0.00105;
%indices = find(g.Vertices(:,3) <= z_thresh);
%model.VertexBC(indices) = vertexBC(Constraint="fixed");

fprintf('BCs and Material Properties set.\n');
results = solve(model);
fprintf('Solved!\n');

% visualize
u = results.Displacement;
figure;
pdeplot3D(model.Mesh,'ColorMapData', u.Magnitude, 'FaceAlpha', 0.3, ...
    'Deformation', u, ...
    'DeformationScaleFactor', 1);
title('Displacement colormap')
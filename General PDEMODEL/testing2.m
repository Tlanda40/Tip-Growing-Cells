% Use importGeometry for PDE
model = createpde(3);  %3 dimensions
importGeometry(model, "HollowHemisphere.step")

% Plot geometry
pdegplot(model.Geometry,'FaceLabels','on');
axis equal;
title('initial geo')

model.Geometry.boundingBox

% Generate tetrahedral mesh, can specify resolution
%mesh = generateMesh(model, 'Hmax', 0.02, 'Hmin', 0.0002);

mesh = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.0002, 'GeometricOrder', 'linear')

%Plot mesh
figure;
pdemesh(model,'FaceAlpha', 0.3);
title('Tetrahedral Mesh of Direct hollow sphere');


function idx = index9(i,j)
    % Map (i,j) to 1..9 index for flattened gradient (same as before)
    idx = 3*(j-1) + i;
end

function Cmat = elasticityTensor3D(location, state)
    E = 1e9;
    nu = 0.3;

    lambda = E*nu / ((1+nu)*(1-2*nu));
    mu = E/(2*(1+nu));

    n = numel(location.x);
    Cmat = zeros(81, n); % 9x9 flattened matrix per point

    % Build 4th order tensor C (3x3x3x3)
    C4 = zeros(3,3,3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    C4(i,j,k,l) = lambda * (i==j)*(k==l) + mu * ((i==k)*(j==l) + (i==l)*(j==k));
                end
            end
        end
    end

    % Reshape C4 to 9x9 matrix according to PDE Toolbox ordering
    % Index mapping from (i,k) to row, (j,l) to column:
    % row = 3*(i-1)+k
    % col = 3*(j-1)+l
    C9 = zeros(9,9);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    row = 3*(i-1) + k;
                    col = 3*(j-1) + l;
                    C9(row, col) = C4(i,j,k,l);
                end
            end
        end
    end

    for idx = 1:n
        Cmat(:, idx) = C9(:); % same isotropic tensor everywhere
    end
end

p = 1e2;           % Pressure

%Helpers for computing triangle normals for faces 3 and 4
face3Nodes = findNodes(mesh, 'region', 'Face', 3);
face3NodeCoords = mesh.Nodes(:, face3Nodes)';
face4Nodes = findNodes(mesh, 'region', 'Face', 4);
face4NodeCoords = mesh.Nodes(:, face4Nodes)';
InnerFaceNodes_combined = [face3NodeCoords; face4NodeCoords];
shp = alphaShape(InnerFaceNodes_combined(:,1), InnerFaceNodes_combined(:,2), InnerFaceNodes_combined(:,3));
[faces, points] = boundaryFacets(shp);
shp.Alpha = 0.01;

% Compute triangle normals
triangleNormals = zeros(size(faces,1),3);
for i = 1:size(faces,1)
    v1 = points(faces(i,1), :);
    v2 = points(faces(i,2), :);
    v3 = points(faces(i,3), :);
    edge1 = v2 - v1;
    edge2 = v3 - v1;
    n = cross(edge1, edge2);
    triangleNormals(i, :) = n / norm(n);
end

% Compute node normals by averaging adjacent triangle normals
nodeNormals = zeros(size(points));
countAdj = zeros(size(points,1),1);
for i = 1:size(faces,1)
    for j = 1:3
        idx = faces(i,j);
        nodeNormals(idx, :) = nodeNormals(idx, :) + triangleNormals(i, :);
        countAdj(idx) = countAdj(idx) + 1;
    end
end
nodeNormals = nodeNormals ./ countAdj; % average
for i=1:size(nodeNormals,1)
    nodeNormals(i,:) = nodeNormals(i,:) / norm(nodeNormals(i,:));
end
% Now nodeNormals is a Nx3 array of unit normals at each surface node.

% These interpolants map positions → normal components
interpNX = scatteredInterpolant(points(:,1), points(:,2), points(:,3), nodeNormals(:,1), 'linear', 'none');
interpNY = scatteredInterpolant(points(:,1), points(:,2), points(:,3), nodeNormals(:,2), 'linear', 'none');
interpNZ = scatteredInterpolant(points(:,1), points(:,2), points(:,3), nodeNormals(:,3), 'linear', 'none');

gFunction = @(location, state) p * [ ...
    interpNX(location.x, location.y, location.z); 
    interpNY(location.x, location.y, location.z); 
    interpNZ(location.x, location.y, location.z)];


%This is where the FE expansion is set up.
specifyCoefficients(model, ...
    'm', 0, ...
    'd', 0, ...
    'c', @elasticityTensor3D, ...
    'a', 0, ...
    'f', @(location, state) repmat([0; 0; 0], 1, numel(location.x)));  %no body force (gravity is irrelevant?)

%structuralBoundaryLoad(model, 'Face', [3,4], 'Pressure', PV_ratio * Vi);
applyBoundaryCondition(model, ...
    'neumann', ...
    'Face', [3,4], ...
    'g', gFunction, ...
    'q', zeros(3,3));
%face1 is outer curve, face3 is inner curve, face2 is outer flat, face4 is
%inner flat

%structuralBC(model,"Face",[2,4],"Constraint","fixed");
applyBoundaryCondition(model, ...
    'dirichlet', ...
    'Face', [2,4], ...
    'u', [0;0;0]);

model.SolverOptions.ReportStatistics = 'on';
model.SolverOptions.MaxIterations = 50;
model.SolverOptions.ResidualTolerance = 1e-6;

numNodes = size(model.Mesh.Nodes, 2);
numDOFs = numNodes * 3;  % 3 equations
u0 = zeros(numDOFs, 1);
setInitialConditions(model, 0)

result = solvepde(model);

u = result.nodalSolution;
umag = sqrt(u.ux.^2 + u.uy.^2 + u.uz.^2);
figure;
pdeplot3D(model,'ColorMapData', umag, 'FaceAlpha', 0.3, ...
    'Deformation', u, ...
    'DeformationScaleFactor', 1);
title('Displacement colormap')
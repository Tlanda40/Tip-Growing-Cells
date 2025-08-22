% Create 3 PDEs for 3D elasticity (u_x, u_y, u_z)
model = createpde(3);

% Geometry: unit cube
importGeometry(model, "HollowHemisphere.step")
pdegplot(model, 'FaceLabels', 'on'); view(3);
%mesh = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.0002, 'GeometricOrder', 'linear')
mesh = generateMesh(model, 'GeometricOrder', 'linear')

% Material properties
E = 1e5;      % Young's modulus in Pascals
nu = 0.3;       % Poisson's ratio
rho = 1000;     % Density in kg/m^3
g = 9.81;       % Gravity in m/s^2

% Lame parameters
lambda = E*nu / ((1 + nu)*(1 - 2*nu));
mu = E / (2*(1 + nu));

% Generate c coefficient matrix (isotropic elasticity tensor)
c = elasticityC3D(mu, lambda);  % see function below

p = -1e5;           % Pressure

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


% Specify PDE coefficients
specifyCoefficients(model, ...
    'm', 0, ...
    'd', 0, ...
    'c', c, ...
    'a', 0, ...
    'f', [0; 0; 0]);

% Boundary conditions: fix bottom face (z = 0)
%applyBoundaryCondition(model, 'dirichlet', 'Face', [2,4], 'u', [0; 0; 0]);


applyBoundaryCondition(model, ...
    'neumann', ...
    'Face', [3,4], ...
    'g', gFunction, ...
    'q', zeros(3,3));

% Assemble global stiffness K and load F immediately
u0 = zeros(3, size(mesh.Nodes, 2));  % initial guess
state.u = u0;  % create the required state structure
state.time = 0;
FEM = assembleFEMatrices(model, 'nullspace', state);
K   = FEM.Kc;
F   = FEM.Fc;

% Pin bottom-half nodes (z<=0)
meshNodes = model.Mesh.Nodes';    % N×3 array of [x,y,z]
bottomIDs   = find(meshNodes(:,3) <= 0.002);  % indices of nodes on bottom half
nN   = size(meshNodes,1);
dofs = [ bottomIDs; bottomIDs + nN; bottomIDs + 2*nN ];
for d = dofs.'
    K(d,:) = 0;  K(:,d) = 0;  K(d,d) = 1;
    F(d)   = 0;
end

% Solve and update the node positions
uvec     = K\F;
ux       = uvec(         1:nN);
uy       = uvec(nN+1:2*nN);
uz       = uvec(2*nN+1:3*nN);
newPts   = meshNodes + [ux,uy,uz];  % N×3

% Make Displacement for visualization
Displacement.ux = ux;
Displacement.uy = uy;
Displacement.uz = uz;

% Plot the final deformed sphere
figure;
umag = sqrt(ux.^2 + uy.^2 + uz.^2);
umag(isnan(umag)) = 0;
pdeplot3D(model,'ColorMapData',umag, 'FaceAlpha', 0.3, 'Deformation', Displacement, 'DeformationScaleFactor', 1);
title('Final displacement magnitude');
axis equal; view(30,20);

function Cmat = elasticityC3D(mu, lambda)
% Returns 3D elasticity tensor in PDE Toolbox format (3x3x3x3)
% for an isotropic linear elastic material
Cmat = zeros(81, 1); % 9x9 flattened matrix per point
c = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                c(i,j,k,l) = lambda * (i==j)*(k==l) + ...
                             mu * ((i==k)*(j==l) + (i==l)*(j==k));
            end
        end
    end
end
C9 = zeros(9,9);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    row = 3*(i-1) + k;
                    col = 3*(j-1) + l;
                    C9(row, col) = c(i,j,k,l);
                end
            end
        end
    end

    for idx = 1:1
        Cmat(:, idx) = C9(:); % same isotropic tensor everywhere
    end

end
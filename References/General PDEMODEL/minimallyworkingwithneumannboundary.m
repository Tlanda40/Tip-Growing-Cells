% Clear everything
clear; clc;

% Create 3D vector-valued PDE model
model = createpde(3);

% Our geometry
gm = importGeometry(model, "HollowHemisphere.step");
model.Geometry = gm;
pdegplot(model, 'FaceLabels', 'on', 'FaceAlpha', 0.4); view(3);
title('Inspect Face IDs');
%mesh = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.0002, 'GeometricOrder', 'linear');
mesh = generateMesh(model, 'GeometricOrder', 'linear');
meshNodes = mesh.Nodes';    % N×3 array of [x,y,z]

% Material properties
E = 1e4;     % Young's modulus
nu = 0.3;    % Poisson's ratio

% Lame parameters
lambda = E * nu / ((1 + nu)*(1 - 2*nu));
mu = E / (2 * (1 + nu));

% Define elasticity tensor
c = elasticityC3D(mu, lambda);

% Set PDE coefficients (no body force)
specifyCoefficients(model, ...
    'm', 0, ...
    'd', 0, ...
    'c', c, ...
    'a', 0, ...
    'f', [0; 0; 0]);

% Use pdegplot to find face labels
figure;
pdegplot(model, 'FaceLabels', 'on'); view(3); title('Face Labels');

% Fix the bottom face (Dirichlet BC)
%applyBoundaryCondition(model, 'dirichlet', 'Face', [2,4], 'u', [0; 0; 0]);

% Assemble FEM Matrices
FEM = assembleFEMatrices(model,'nullspace');
K   = FEM.Kc;
F   = FEM.Fc;

bottomIDs   = find(meshNodes(:,3) <= 0.005);  % indices of nodes on bottom half
nN   = size(meshNodes,1);
dofs = [ bottomIDs; bottomIDs + nN; bottomIDs + 2*nN ];
for d = dofs.'
    K(d,:) = 0;  K(:,d) = 0;  K(d,d) = 1;
    F(d)   = 0;
end

% Apply downward pressure on the top face (Neumann BC)
pressure = 1e6;  % Pascals
applyBoundaryCondition(model, 'neumann', ...
    'Face', [3,4], ...
    'g', [0; 0; -pressure], ...
    'q', zeros(3,3));

% Get displacement result
uvec     = K\F;
ux       = uvec(1:nN);
uy       = uvec(nN+1:2*nN);
uz       = uvec(2*nN+1:3*nN);

% Make Displacement for visualization
Displacement.ux = ux;
Displacement.uy = uy;
Displacement.uz = uz;

% Plot the final deformed sphere
figure;
umag = sqrt(ux.^2 + uy.^2 + uz.^2);
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
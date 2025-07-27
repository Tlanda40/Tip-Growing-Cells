%  Parameters
R0       = 1.0;        % Initial sphere radius
E        = 50e3;      % Young's modulus (Pa)
nu       = 0.3;        % Poisson's ratio
pressure = 1e3;        % Internal pressure (Pa)
nIter    = 5;          % Number of iterations
h0       = R0/6;       % Initial grid spacing

model = createpde('structural','static-solid');
gm = multisphere([0.099 0.1],Void=[true,false]);
model.Geometry = gm;

% Plot geometry
figure; % Plot with face labels
pdegplot(model.Geometry,'FaceLabels','on', 'FaceAlpha', 0.2);
axis equal;
title('initial geo')
figure; % Plot with cell labels
pdegplot(model, 'CellLabels', 'on', 'FaceAlpha', 0.2);
view(3); axis equal; title('Cell IDs');
model.Geometry.NumCells

%mesh = generateMesh(model,'Hmax', 0.003, 'Hmin', 0.002,'GeometricOrder', 'linear');
mesh = generateMesh(model,'Hmax', 0.009, 'Hmin', 0.007,'GeometricOrder', 'linear');

% Plot mesh
figure;
pdemesh(model,'FaceAlpha', 0.3);
title('Tetrahedral Mesh of Direct hollow sphere');

structuralProperties(model, ...
    'YoungsModulus',E, ...
    'PoissonsRatio',nu);

structuralBoundaryLoad(model, ...
  'Face', 1, ...
  'Pressure', pressure);

% Assemble global stiffness K and load F immediately
FEM = assembleFEMatrices(model,'nullspace');
K   = FEM.Kc;
F   = FEM.Fc;

% Pin bottom-half nodes (z<=0)
meshNodes = model.Mesh.Nodes';    % N×3 array of [x,y,z]
bottomIDs   = find(meshNodes(:,3) <= 0);  % indices of nodes on bottom half
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
pdeplot3D(model,'ColorMapData',umag, 'FaceAlpha', 0.3, 'Deformation', Displacement, 'DeformationScaleFactor', 1);
title('Final displacement magnitude');
axis equal; view(30,20);
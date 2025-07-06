% sphere_iterative_pressure_V2.m
% Simulate a sphere under internal pressure,
% remeshing only by re-tetrahedralizing your displaced nodes.
%% Parameters
R0       = 1.0;        % Initial sphere radius
E        = 200e3;      % Young's modulus (Pa)
nu       = 0.3;        % Poisson's ratio
pressure = 1e1;        % Internal pressure (Pa)
nIter    = 5;          % Number of iterations
h0       = R0/6;       % Initial grid spacing
for k = 1:2
    fprintf('Iteration %d/%d: h = %.3f\n', k, nIter, h0);
    %% 1) Build Delaunay mesh of the current sphere
    [X,Y,Z] = ndgrid(-R0:h0:R0);
    P = [X(:),Y(:),Z(:)];
    P = P(sum(P.^2,2)<=R0^2, :);
    DT = delaunayTriangulation(P);
    nodes = DT.Points';           % 3×N
    elems = DT.ConnectivityList'; % 4×M
    %% 2) Create PDE model + import that mesh as BOTH geometry & mesh
    model = createpde('structural','static-solid');
    geometryFromMesh(model, nodes, elems);
    structuralProperties(model, ...
        'YoungsModulus',E, ...
        'PoissonsRatio',nu);
    %% 3) Apply internal pressure on *all* boundary faces
    structuralBoundaryLoad(model, ...
      'Face',1:model.Geometry.NumFaces, ...
      'Pressure',-pressure);
    %% 4) Assemble global stiffness K and load F immediately
    FEM = assembleFEMatrices(model,'nullspace');
    K   = FEM.Kc;
    F   = FEM.Fc;
    %% 4) Pin left-half nodes (x<=0)
    meshNodes = model.Mesh.Nodes';    % N×3 array of [x,y,z]
    leftIDs   = find(meshNodes(:,1) <= 0);  % indices of nodes on left half
    nN   = size(meshNodes,1);
    dofs = [ leftIDs; leftIDs + nN; leftIDs + 2*nN ];
    for d = dofs.'
        K(d,:) = 0;  K(:,d) = 0;  K(d,d) = 1;
        F(d)   = 0;
    end
    %% 6) Solve and update the node positions
    uvec     = K\F;
    ux       = uvec(         1:nN);
    uy       = uvec(nN+1:2*nN);
    uz       = uvec(2*nN+1:3*nN);
    newPts   = meshNodes + [ux,uy,uz];  % N×3
    %% 7) Prepare for next iter: new radius + finer h
    R0 = max(sqrt(sum(newPts.^2,2)));
    h0 = h0 * 0.9;
end
%% Plot the final deformed sphere
figure;
mag = sqrt(ux.^2 + uy.^2 + uz.^2);
pdeplot3D(model,'ColorMapData',mag);
title('Final displacement magnitude');
axis equal; view(30,20);
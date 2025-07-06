% Use importGeometry for PDE
model = createpde('structural','static-solid');
gm = multisphere([0.99 1],Void=[true,false])

model.Geometry = gm

% Plot geometry
pdegplot(model.Geometry,'FaceLabels','on');
axis equal;
title('initial geo')

model.Geometry.boundingBox

% Generate tetrahedral mesh, can specify resolution
%mesh = generateMesh(model, 'Hmax', 0.002, 'Hmin', 0.0002);

mesh = generateMesh(model, 'Hmin', 0.1)

figure;
pdemesh(model,'FaceAlpha', 0.3);
title('Tetrahedral Mesh of Direct hollow sphere');


%insert structural properties. For now, it will be some kind of basic
%thing, with uniform structural properties (as opposed to properties that
%are a function of arc length.


E = 200e3;            % Base Young's modulus (Pa), applies to whole cylinder
nu = 0.3;           % Poisson's ratio

structuralProperties(model, 'YoungsModulus', E, 'PoissonsRatio', nu);


%{
%here, we want to get all points at z=0 and less.

bottomHalf = findNodes(mesh,"box",[-1 1],[-1 1],[-1 0])

model.NumVertices

for i = 1:size(bottomHalf,2)
    %bottomHalf(1,i)
    structuralBC(model,'Vertex',bottomHalf(1,i),'Constraint','fixed');
    %model.VertexBC(bottomHalf(i)) = vertexBC(Constraint="fixed")
end

%model.VertexBC(VertexID) = vertexBC(Name=Value)
%}

structuralBoundaryLoad(model, 'Face', 1, 'Pressure', -1e1);
%Assemble global stiffness K and load F immediately
    FEM = assembleFEMatrices(model,'nullspace');
    K   = FEM.Kc;
    F   = FEM.Fc;
% 4) Pin left-half nodes (x<=0)
    meshNodes = model.Mesh.Nodes';    % N×3 array of [x,y,z]
    leftIDs   = find(meshNodes(:,1) <= 0);  % indices of nodes on left half
    nN   = size(meshNodes,1);
    dofs = [ leftIDs; leftIDs + nN; leftIDs + 2*nN ];
    for d = dofs.'
        K(d,:) = 0;  K(:,d) = 0;  K(d,d) = 1;
        F(d,:)   = 0;
    end
result = solve(model);
pdeplot3D(model,'ColorMapData',result.Displacement.Magnitude)
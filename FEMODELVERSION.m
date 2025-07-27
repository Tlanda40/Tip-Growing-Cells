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
%size(model.Mesh.Nodes)

% Plot mesh
figure;
pdemesh(model,'FaceAlpha', 0.3);
title('Tetrahedral Mesh of Direct hollow sphere');

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

fprintf('BCs and Material Properties set.');
results = solve(model);
fprintf('Solved!');

% visualize
u = results.Displacement;
figure;
pdeplot3D(model.Mesh,'ColorMapData', u.Magnitude, 'FaceAlpha', 0.3, ...
    'Deformation', u, ...
    'DeformationScaleFactor', 1);
title('Displacement colormap')
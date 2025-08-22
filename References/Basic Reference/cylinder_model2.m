%cylinder_model2.m

%%
% Parameters
radius = 1;      % Radius of the cylinder (m)
height = 0.5;     % Height of the cylinder (m)
E = 1e1;        % Young's modulus (Pa)% 0.6MPa
nu = 0.3;       % Poisson's ratio

% Spring stiffnesses (N/m³)
kN_top = 1e5;   kS_top = 1e5;    % Top face
kN_bottom = 1e5; kS_bottom = 1e5;% Bottom face
kN_side = 1e25;  kS_side = 1e25;   % Side face

pressureValue = -1; % Pressure (Pa) %10

% Create PDE model
model = createpde('structural', 'static-solid');

% Create cylinder geometry
gm = multicylinder(radius, height);
model.Geometry = gm;

% Generate mesh
generateMesh(model, 'Hmax', 0.3);

% Assign material properties
structuralProperties(model, 'YoungsModulus', E, 'PoissonsRatio', nu);

% Apply spring boundary conditions
structuralBoundaryLoad(model, 'Face', 2, 'TranslationalStiffness', [kS_top, kN_top, kN_top]);%top
structuralBoundaryLoad(model, 'Face', 1, 'TranslationalStiffness', [kS_bottom, kN_bottom, kN_bottom]);%bottom
structuralBoundaryLoad(model, 'Face', 3, 'TranslationalStiffness', [kS_side, kN_side, kN_side]);%side

%structuralBC(model,"Vertex",1,"YDisplacement",0);
%structuralBC(model,"Vertex",2,"YDisplacement",0);
%structuralBC(model,"Vertex",2,"XDisplacement",0);

% Apply pressure on side wall
structuralBoundaryLoad(model, 'Face', [1:2], 'Pressure', pressureValue);

% Solve the model
result = solve(model);

% Plot the deformed shape
figure;
pdeplot3D(model, ...
    'ColorMapData', result.VonMisesStress, ...
    'Deformation', result.Displacement, ...
    'DeformationScaleFactor', 10);
%title('Deformed Cylinder under Pressure');


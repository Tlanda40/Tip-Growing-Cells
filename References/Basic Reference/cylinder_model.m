
% Parameters
E = 1e9;            % Base Young's modulus (Pa), applies to whole cylinder
nu = 0.3;           % Poisson's ratio
p_value = 1e5;      % Internal pressure (Pa)
n_steps = 100;
pressures = linspace(p_value,p_value*1000,100);

% Create the model
model = createpde('structural','static-solid');

% Geometry: full cylinder
R = 1;              % radius
H = 2;              % height
gm = multicylinder(R, H);
model.Geometry = gm;

% Assign material properties to the entire domain
structuralProperties(model, 'YoungsModulus', E, 'PoissonsRatio', nu);

% Plot geometry with face labels
figure;
pdegplot(model, 'FaceLabels', 'on');
title('Cylinder with Face Labels');

% Apply boundary conditions
% Top face (F2) → fixed
structuralBC(model, 'Face', 2, 'Constraint', 'fixed');

% Side wall (F3) → apply pressure
% You want this to be stiff but still deformable, which is already handled via the high E
%structuralBoundaryLoad(model, 'Face', 3, 'Pressure', p_value);
structuralBoundaryLoad(model, 'Face', 1, 'Pressure', p_value);

% Bottom face (F1) → more compliant
% We simulate this by letting it move freely, but we can later restrict it or apply a soft support if needed
% (Optionally, you could simulate extra compliance using a spring foundation:
%k_soft = (1e9)/10;
structuralBoundaryLoad(model, ...
                          "Face",1, ...
                          "TranslationalStiffness",[1 1 1])

% Mesh and solve
generateMesh(model, 'Hmax', 0.5);
result = solve(model);

% Plot displacement magnitude
figure;
pdeplot3D(model, 'ColorMapData', result.Displacement.Magnitude);
title('Displacement Magnitude');

%%
figure;
for i = 1:n_steps
    % Remove existing boundary load and apply new pressure
    %structuralBoundaryLoad(model, 'Face', 3, 'Pressure', pressures(i));
    
    structuralBoundaryLoad(model, ...
                          "Face",1, ...
                          "Pressure",pressures(i),...
                          "TranslationalStiffness",[1 1 1])

    %structuralBoundaryLoad(model, 'Face', 1, 'Pressure', pressures(i));


    % Solve the model
    result = solve(model);
    
    % Clear the axis and replot the updated solution:
    clf;
    pdeplot3D(model, 'ColorMapData', result.Displacement.Magnitude, 'Deformation', result.Displacement);
    title(sprintf('Step %d/%d - Pressure = %.0f Pa', i, n_steps, pressures(i)));
    drawnow;
    
    % Optionally, pause a bit to slow down the animation:
    pause(0.1);
end

gm = multisphere([0.99 1],Void=[true,false])
pdegplot(gm,CellLabels="on",FaceAlpha=0.5)


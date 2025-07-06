function vols = tetVolume(nodes, elements)
% nodes: 3×N array of node coordinates (x; y; z)
% elements: 4×M array of tetrahedral elements (indices into nodes)

    p1 = nodes(:, elements(1, :));
    p2 = nodes(:, elements(2, :));
    p3 = nodes(:, elements(3, :));
    p4 = nodes(:, elements(4, :));

    % Compute signed volume of each tetrahedron
    v = dot(cross(p2 - p1, p3 - p1, 1), p4 - p1, 1) / 6;

    vols = v;  % 1×M array of signed volumes
end

function Extensibility = ext(location)
    x = location.x
    y = location.y
    A1=0.27;
    A2=10;
    w1=0.1;
    w2=0.1;
    g1x = A1*exp(-(x/w1).^2);
    g2x = A2*exp(-(x/w2).^2);
    Gmixx = g1x + g2x;

    g1y = A1*exp(-(y/w1).^2);
    g2y = A2*exp(-(y/w2).^2);
    Gmixy = g1y + g2y;

    G=Gmixx.*Gmixy;

    Extensibility = G;
end

% Use importGeometry for PDE
model = createpde('structural','static-solid');
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

%insert structural properties. For now, it will be some kind of basic
%thing, with uniform structural properties (as opposed to properties that
%are a function of arc length.

E = 1e5;         % Base Young's modulus (Pa)
nu = 0.3;           % Poisson's ratio
rho = 1000;         % Mass Density, only relevant in a transient-solid model
p = 1e2;            % Pressure

%Find initial volume enclosed (volume of void inside shell)
face3Nodes = findNodes(mesh, 'region', 'Face', 3);
face3NodeCoords = mesh.Nodes(:, face3Nodes)';
face4Nodes = findNodes(mesh, 'region', 'Face', 4);
face4NodeCoords = mesh.Nodes(:, face4Nodes)';
InnerFaceNodes_combined = [face3NodeCoords; face4NodeCoords];
shp = alphaShape(InnerFaceNodes_combined(:,1), InnerFaceNodes_combined(:,2), InnerFaceNodes_combined(:,3));
[faces, points] = boundaryFacets(shp);
shp.Alpha = 0.01;
figure;
h = plot(shp)
h.FaceAlpha = 0.3;
title('Inner void')
axis equal;

Vi = volume(shp); %this is the initial volume enclosed
%V = mesh.volume() (this is the volume of the shell itself)

%Here, we normalize pressure with respect to volume
PV_ratio = p / Vi;

%This is where the FE expansion is set up.
structuralProperties(model, 'YoungsModulus', E, 'PoissonsRatio', nu);
structuralBoundaryLoad(model, 'Face', [3,4], 'Pressure', PV_ratio * Vi);
%face1 is outer curve, face3 is inner curve, face2 is outer flat, face4 is
%inner flat
structuralBC(model,"Face",[2,4],"Constraint","fixed");
%structuralIC(model, 'Velocity', [0; 0; 0], 'Displacement', [0; 0; 0]);
%tlist = 0:0.1:0.1;  % simulate from t = 0 to 1 s
result = solve(model);
u = result.Displacement;
umag = sqrt(u.ux.^2 + u.uy.^2 + u.uz.^2);
figure;
pdeplot3D(model,'ColorMapData', umag, 'FaceAlpha', 0.3, ...
    'Deformation', u, ...
    'DeformationScaleFactor', 1);
title('Displacement colormap')

displacedNodes = mesh.Nodes + [u.ux'; u.uy'; u.uz'];
elements = mesh.Elements;
elementsLinear = elements(1:4, :);
vols = tetVolume(displacedNodes, elementsLinear);
bad = vols < 0;
% Flip nodes 3 and 4 for inverted tets (common fix)
elementsLinear(:, bad) = elementsLinear([1 2 4 3], bad);

%Now, we want to 

for i = 1:5
    newModel = createpde('structural','static-solid');
    geometryFromMesh(newModel, displacedNodes, elementsLinear);
    mesh = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.0002, 'GeometricOrder', 'linear')
    face3Nodes = findNodes(mesh, 'region', 'Face', 3);
    face3NodeCoords = mesh.Nodes(:, face3Nodes)';
    face4Nodes = findNodes(mesh, 'region', 'Face', 4);
    face4NodeCoords = mesh.Nodes(:, face4Nodes)';
    InnerFaceNodes_combined = [face3NodeCoords; face4NodeCoords];
    shp = alphaShape(InnerFaceNodes_combined(:,1), InnerFaceNodes_combined(:,2), InnerFaceNodes_combined(:,3));
    [faces, points] = boundaryFacets(shp);
    shp.Alpha = 0.01;
    Vn = volume(shp);
    structuralProperties(newModel, 'YoungsModulus', E, 'PoissonsRatio', nu);
    structuralBoundaryLoad(newModel, 'Face', [3,4], 'Pressure', PV_ratio * Vn);
    structuralBC(newModel,"Face",[2,4],"Constraint","fixed");
    result = solve(newModel);
    u = result.Displacement;
    umag = sqrt(u.ux.^2 + u.uy.^2 + u.uz.^2);
    figure;
    pdeplot3D(newModel,'ColorMapData', umag, 'FaceAlpha', 0.3, ...
        'Deformation', u, ...
        'DeformationScaleFactor', 1);
    title('Displacement colormap')
    displacedNodes = mesh.Nodes + [u.ux'; u.uy'; u.uz'];
    elements = mesh.Elements;
    elementsLinear = elements(1:4, :);
    vols = tetVolume(displacedNodes, elementsLinear);
    bad = vols < 0;
    % Flip nodes 3 and 4 for inverted tets (common fix)
    elementsLinear(:, bad) = elementsLinear([1 2 4 3], bad);
end

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
    x = location.x;
    y = location.y;
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

%insert structural properties. For now, it will be some kind of basic
%thing, with uniform structural properties (as opposed to properties that
%are a function of arc length.

%{
function Cmat = elasticityTensor3D(location, state)

    % Spatially varying modulus
    E = ext(location);  %refer to ext function; this is the Young's Modulus
    nu = 0.3; %poisson's ratio

    % Lame parameters
    lambda = (E .* nu) ./ ((1 + nu) .* (1 - 2 * nu));
    mu     = E ./ (2 * (1 + nu));

    n = numel(location.x);
    Cmat = zeros(36,n);  % 6x6 tensor at each point

    for i = 1:n
        l = lambda(i);
        m = mu(i);
        C = [l+2*m, l,     l,     0, 0, 0;
             l,     l+2*m, l,     0, 0, 0;
             l,     l,     l+2*m, 0, 0, 0;
             0,     0,     0,     m, 0, 0;
             0,     0,     0,     0, m, 0;
             0,     0,     0,     0, 0, m];
        Cmat(:,i) = C(:);  % flatten 6x6 matrix to 36x1 vector column-wise
    end
end
%}

function idx = index9(i,j)
    % Map (i,j) to 1..9 index for flattened gradient (same as before)
    idx = 3*(j-1) + i;
end

function Cmat = elasticityTensor3D(location, state)

    % Spatially varying modulus
    E = ext(location);  %refer to ext function; this is the Young's Modulus
    nu = 0.3; %poisson's ratio

    n = numel(location.x);
    Cmat = zeros(81, n); % 9x9 flattened matrix per point

    for i = 1:n
        Ei = E(i);
        lambda = (Ei * nu) / ((1 + nu)*(1 - 2*nu));
        mu = Ei / (2 * (1 + nu));
        
        % Construct the 9x9 tensor C at point i
        C = zeros(9,9);
        
        % Populate C according to PDE Toolbox indexing:
        % Ordering of indices (row/col): 
        % (1,1) (2,1) (3,1) (1,2) (2,2) (3,2) (1,3) (2,3) (3,3)
        % Using symmetry and isotropy:
        
        % Fill diagonal blocks (normal stress-strain)
        for a = 1:3
            for b = 1:3
                idx1 = index9(a,a); % e.g. (1,1), (2,2), (3,3)
                idx2 = index9(b,b);
                C(idx1, idx2) = lambda;
            end
            idx = index9(a,a);
            C(idx, idx) = lambda + 2*mu;
        end
        
        % Fill shear terms:
        shearIdx = [2,3,4,6,7,8]; % indices corresponding to shear components in the 9x9 matrix
        
        for s = shearIdx
            C(s,s) = mu;
        end
        
        % Flatten column-wise and store:
        Cmat(:, i) = C(:);
    end
end


rho = 1000;        % Mass Density, only relevant in a transient-solid model
p = 1e2;           % Initial Pressure


%Find initial volume enclosed (volume of void inside shell)
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

%{
figure;
h = plot(shp)
h.FaceAlpha = 0.3;
title('Inner void')
axis equal;

Vi = volume(shp); %this is the initial volume enclosed
%V = mesh.volume() (this is the volume of the shell itself)

%Here, we normalize pressure with respect to volume
PV_ratio = p / Vi;
%}

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

result = solvepde(model);
u = result.nodalSolution;
umag = sqrt(u.ux.^2 + u.uy.^2 + u.uz.^2);
figure;
pdeplot3D(model,'ColorMapData', umag, 'FaceAlpha', 0.3, ...
    'Deformation', u, ...
    'DeformationScaleFactor', 1);
title('Displacement colormap')

%{

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
%}
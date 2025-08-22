% Helper function to calculate element volume
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
    y = 5000 ./ Ext;
end

% Helper function tet 4xN to face 3x4N
function Y = tetToFace(tets)
    Y = zeros(3,4*size(tets,2));
    for i = 1:size(tets,2)
        Y(:,(4*i)-3) = [tets(1,i); tets(2,i); tets(3,i)];
        Y(:,(4*i)-2) = [tets(1,i); tets(2,i); tets(4,i)];
        Y(:,(4*i)-1) = [tets(1,i); tets(3,i); tets(4,i)];
        Y(:,(4*i)) = [tets(2,i); tets(3,i); tets(4,i)];
    end
end

% Helper Function takes mxN elts to only columns with val
function Y = selectEltsWithVal(elts, val)
    logicalElts = elts == val;
    for j = 1:size(elts,2)
        for i = 1:size(elts,1)
            if logicalElts(i,j) == 1
                logicalElts(:,j) = ones(size(elts,1),1);
                break;
            end
        end
    end
    Y = elts .* logicalElts;
    Y(:,all(Y == 0))=[];
end

% Helper Function to remove all columns that have a duplicate, and their
% duplicate (assumes no nonzero matrix elements)
function Y = rmDupeCols(M)
    Y = sort(M);
    for j = 1:(size(Y,2)-1)
        if Y(1,j) ~= 0
            for i = (j+1):size(Y,2)
                if isequal(Y(:,j), Y(:,i))
                    Y(:,j) = zeros(size(Y,1),1);
                    Y(:,i) = zeros(size(Y,1),1);
                    break;
                end
            end
        end             
    end
    Y(:,all(Y == 0))=[];
end

% Helper function to compute triangle normals (expects 3xN and 3xM),
% returns 3xN
function Y = triangleNormals(faces, points)
    Y = zeros(size(faces));
    for j = 1:size(faces,2)
        v1 = points(:, faces(1,j));
        v2 = points(:, faces(2,j));
        v3 = points(:, faces(3,j));
        edge1 = v2 - v1;
        edge2 = v3 - v1;
        n = cross(edge1, edge2);
        Y(:,j) = n / norm(n);
    end
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
model = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.00028, 'GeometricOrder', 'linear');
mesh = model.Mesh;
%size(model.Mesh.Nodes)

% Plot mesh
figure;
pdemesh(model,'FaceAlpha', 0.3);
title('Tetrahedral Mesh of Direct hollow sphere');

%insert structural properties. For now, it will be some kind of basic
%thing, with uniform structural properties (as opposed to properties that
%are a function of arc length.

nu = 0.3;           % Poisson's ratio
p = 1e2;            % Pressure
p1 = -4e2;   
p2 = -3.948e2;
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

face3Nodes = findNodes(mesh, 'region', 'Face', 1);
face3NodeCoords = mesh.Nodes(:, face3Nodes)';
face4Nodes = findNodes(mesh, 'region', 'Face', 2);
face4NodeCoords = mesh.Nodes(:, face4Nodes)';
P = [face3NodeCoords; face4NodeCoords];

shp = alphaShape(P(:,1), P(:,2), P(:,3));
shp.Alpha = criticalAlpha(shp, 'all-points')+0.0001;
figure;
h = plot(shp);
h.FaceAlpha = 0.3;
Vi = volume(shp);

shpinner = alphaShape(mesh.Nodes');
shpinner.Alpha = criticalAlpha(shpinner, 'all-points');
Vinner = volume(shpinner);

PV_ratio = p / Vi;

PV_ratio1init = p1 / Vinner;
PV_ratio2init = p2 / Vinner;

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

displacedNodes = mesh.Nodes + [u.ux'; u.uy'; u.uz'];
elements = mesh.Elements;
vols = tetVolume(displacedNodes, elements);
bad = vols < 0;
elements(:, bad) = elements([1 2 4 3], bad);
warning('off', 'MATLAB:alphaShape:DupPointsBasicWarnId');
timeSteps = 100;
%frames(timeSteps) = struct('cdata', [], 'colormap', []);
%v = VideoWriter('my_animation.mp4', 'MPEG-4');
%v.FrameRate = 10;
%open(v);
figure;
i = 1;
thicknessTick = 0;
while i < timeSteps
    i
    model = femodel("AnalysisType","structuralStatic");
    model.Geometry = fegeometry(displacedNodes', elements');
    %{
    figure;
    pdegplot(model.Geometry, FaceLabels='on', FaceAlpha=0.5);
    axis equal;
    title(['loop geometry ', num2str(i)])
    %}
    %model = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.00014, 'GeometricOrder', 'linear');
    model = generateMesh(model, 'Hmax', 0.00015, 'GeometricOrder', 'linear');
    mesh = model.Mesh;
    %{
    figure;
    pdeplot3D(mesh, 'FaceAlpha', 0.3);
    %}
    faceNodeIdx = cell(1, 4);
    maxvec = [0, 0, 0, 0];
    faces = [1, 2, 3, 4];
    for j = 1:4
        faceNodeIdx{j} = findNodes(mesh, 'region', 'Face', j);
        maxvec(j) = max(mesh.Nodes(3, faceNodeIdx{j}));
    end
    [~, idx] = sort(maxvec, 'Descend');
    facesSorted = faces(idx);
    P = mesh.Nodes(:,[faceNodeIdx{facesSorted(2)}, faceNodeIdx{facesSorted(3)}]);
    shp = alphaShape(P');
    shp.Alpha = criticalAlpha(shp, 'all-points')+0.0001;
    %{
    figure;
    h = plot(shp);
    h.FaceAlpha = 0.3;
    %}
    V = volume(shp);
    
    shpinner = alphaShape(mesh.Nodes');
    shpinner.Alpha = criticalAlpha(shpinner, 'all-points');
    Vinner = volume(shpinner);
    
    % max z pole finder (worse)
    %[~, poleIdx] = max(mesh.Nodes(3,:));

    % projected pole finder (better)
    tempOuterCoords = mesh.Nodes([1, 2],faceNodeIdx{facesSorted(1)});
    poleIdxOuter = knnsearch(tempOuterCoords',[0,0]);
    poleIdx = faceNodeIdx{facesSorted(1)}(1,poleIdxOuter);

    poleCoords = mesh.Nodes(:,poleIdx);

    %{
    % Extracting the pole boundary faces for normal calculation
    PoleElements = selectEltsWithVal(mesh.Elements, poleIdx);
    PoleElementsFaces = tetToFace(PoleElements);
    PoleFaces = selectEltsWithVal(PoleElementsFaces, poleIdx);
    PoleBoundaryFaces = rmDupeCols(PoleFaces);
    PoleBFNormals = triangleNormals(PoleBoundaryFaces, mesh.Nodes);
    avgNormal = mean(PoleBFNormals, 2);
    avgNormal = avgNormal / norm(avgNormal);
    %}

    tempInnerCoords = mesh.Nodes(:,faceNodeIdx{facesSorted(2)});
    innerIdx = knnsearch(tempInnerCoords',poleCoords');
    innerCoords = tempInnerCoords(:,innerIdx);
    %thickness = abs(dot((poleCoords - innerCoords), avgNormal))
    thickness = abs(dot((poleCoords - innerCoords), [0; 0; 1]))

    model.MaterialProperties.PoissonsRatio = nu;
    model.MaterialProperties.YoungsModulus = @(location,state) youngMod(ext(beta0, location.x, location.y));
    
    if thickness < 2e-05
        if thicknessTick ~= 0
            PV_ratio1 = PV_ratio1 * 2;
            PV_ratio2 = PV_ratio1 * 2;
        else
            PV_ratio1 = PV_ratio1init;
            PV_ratio2 = PV_ratio2init;
        end
        model.FaceLoad(facesSorted(1)) = faceLoad(Pressure = PV_ratio2 * Vinner); 
        model.FaceBC(facesSorted(2)) = faceBC(Constraint="fixed");
        mvForward = 0;
        thicknessTick = thicknessTick + 1;
    else
        model.FaceLoad(facesSorted(2)) = faceLoad(Pressure = PV_ratio * V);
        mvForward = 1;
        thicknessTick = 0;
    end
    
    %model.FaceLoad(facesSorted(2)) = faceLoad(Pressure = PV_ratio * V);
    %mvForward = 1;
    model.FaceBC([facesSorted(3),facesSorted(4)]) = faceBC(Constraint="fixed");
    results = solve(model);
    u = results.Displacement;
    %figure;
    clf;
    pdeplot3D(model.Mesh,'ColorMapData', u.Magnitude, 'FaceAlpha', 0.3, ...
        'Deformation', u, ...
        'DeformationScaleFactor', 1);
    title(['Loop displacement colormap ', num2str(i)]);
    drawnow;
    %frames(i) = getframe(gcf);
    %frame = getframe(gcf);
    %writeVideo(v, frame);
    displacedNodes = model.Mesh.Nodes + [u.ux'; u.uy'; u.uz'];
    elements = model.Mesh.Elements;
    if mvForward == 1
        i = i + 1;
    end
end
warning('on', 'MATLAB:alphaShape:DupPointsBasicWarnId');
%close(v);


%{
model2 = femodel("AnalysisType","structuralStatic");
model2.Geometry = fegeometry(displacedNodes', elements');
figure;
pdegplot(model2.Geometry, FaceLabels='on', FaceAlpha=0.5);
axis equal;
title('initial single-domain geometry 2')
model2 = generateMesh(model2, 'Hmax', 0.0003, 'Hmin', 0.00028, 'GeometricOrder', 'linear');
mesh = model2.Mesh;

faceNodeCoords = cell(1, 4);
for i = 1:4
    tempNodes = findNodes(mesh, 'region', 'Face', i);
    faceNodeCoords{i} = mesh.Nodes(:, tempNodes)';
end
maxvec = [0, 0, 0, 0];
for i = 1:4
    maxvec(i) = max(faceNodeCoords{i}(:,3));
end
faces = [1, 2, 3, 4];
[maxvecSorted, idx] = sort(maxvec, 'Descend');
facesSorted = faces(idx)

P = [faceNodeCoords{facesSorted(2)}; faceNodeCoords{facesSorted(3)}];



shp2 = alphaShape(P(:,1), P(:,2), P(:,3));
shp2.Alpha = criticalAlpha(shp2, 'all-points')+0.0001
figure;
h = plot(shp2);
h.FaceAlpha = 0.3;
V2 = volume(shp2)


model2.MaterialProperties.PoissonsRatio = nu;
model2.MaterialProperties.YoungsModulus = @(location2,state) youngMod(ext(beta0, location2.x, location2.y));
model2.FaceLoad(facesSorted(2)) = faceLoad(Pressure = PV_ratio * V2);
model2.FaceBC([facesSorted(3),facesSorted(4)]) = faceBC(Constraint="fixed");
results2 = solve(model2);
u2 = results2.Displacement;
figure;
pdeplot3D(model2.Mesh,'ColorMapData', u2.Magnitude, 'FaceAlpha', 0.3, ...
    'Deformation', u2, ...
    'DeformationScaleFactor', 1);
title('Displacement colormap')

displacedNodes2 = model2.Mesh.Nodes + [u2.ux'; u2.uy'; u2.uz'];
elements2 = model2.Mesh.Elements;
vols = tetVolume(displacedNodes2, elements2);
bad = vols < 0;
elements2(:, bad) = elements2([1 2 4 3], bad);
model3 = femodel("AnalysisType","structuralStatic");
model3.Geometry = fegeometry(displacedNodes2', elements2');
figure;
pdegplot(model3.Geometry, FaceLabels='on', FaceAlpha=0.5);
axis equal;
title('initial single-domain geometry 2')

model3 = generateMesh(model3, 'Hmax', 0.0003, 'Hmin', 0.00028, 'GeometricOrder', 'linear');
mesh = model3.Mesh;
figure;
pdeplot3D(mesh, 'FaceAlpha', 0.3);

face1Nodes = findNodes(mesh, 'region', 'Face', 1);
face1NodeCoords = mesh.Nodes(:, face1Nodes)';
face4Nodes = findNodes(mesh, 'region', 'Face', 4);
face4NodeCoords = mesh.Nodes(:, face4Nodes)';
P = [face1NodeCoords; face4NodeCoords];
%{
[k1, av1] = convhull(P);
V3 = av1
figure;
trisurf(k1, P(:,1), P(:,2), P(:,3), ...
    'FaceColor', 'cyan', 'FaceAlpha', 0.7);
axis equal;
%}

shp3 = alphaShape(P(:,1), P(:,2), P(:,3));
shp3.Alpha = criticalAlpha(shp3, 'all-points')+0.0001
figure;
h = plot(shp3);
h.FaceAlpha = 0.3;
V3 = volume(shp3)


model3.MaterialProperties.PoissonsRatio = nu;
model3.MaterialProperties.YoungsModulus = @(location3,state) youngMod(ext(beta0, location3.x, location3.y));
model3.FaceLoad(1) = faceLoad(Pressure = PV_ratio * V3);
model3.FaceBC([3,4]) = faceBC(Constraint="fixed");
results3 = solve(model3);
u3 = results3.Displacement;
figure;
pdeplot3D(model3.Mesh,'ColorMapData', u3.Magnitude, 'FaceAlpha', 0.3, ...
    'Deformation', u3, ...
    'DeformationScaleFactor', 1);
title('Displacement colormap')

displacedNodes3 = model3.Mesh.Nodes + [u3.ux'; u3.uy'; u3.uz'];
elements3 = model3.Mesh.Elements;
vols = tetVolume(displacedNodes3, elements3);
bad = vols < 0;
elements3(:, bad) = elements3([1 2 4 3], bad);
model4 = femodel("AnalysisType","structuralStatic");
model4.Geometry = fegeometry(displacedNodes3', elements3');
figure;
pdegplot(model4.Geometry, FaceLabels='on', FaceAlpha=0.5);
axis equal;
title('initial single-domain geometry 2')
model4 = generateMesh(model4, 'Hmax', 0.0003, 'Hmin', 0.00028, 'GeometricOrder', 'linear');
mesh = model4.Mesh;
figure;
pdeplot3D(mesh, 'FaceAlpha', 0.3);

face2Nodes = findNodes(mesh, 'region', 'Face', 2);
face2NodeCoords = mesh.Nodes(:, face2Nodes)';
face3Nodes = findNodes(mesh, 'region', 'Face', 3);
face3NodeCoords = mesh.Nodes(:, face3Nodes)';
P = [face2NodeCoords; face3NodeCoords];
%{
[k1, av1] = convhull(P);
V4 = av1
figure;
trisurf(k1, P(:,1), P(:,2), P(:,3), ...
    'FaceColor', 'cyan', 'FaceAlpha', 0.7);
axis equal;
%}

shp4 = alphaShape(P(:,1), P(:,2), P(:,3));
shp4.Alpha = criticalAlpha(shp4, 'all-points')+0.0001
figure;
h = plot(shp4);
h.FaceAlpha = 0.3;
V4 = volume(shp4)


model4.MaterialProperties.PoissonsRatio = nu;
model4.MaterialProperties.YoungsModulus = @(location4,state) youngMod(ext(beta0, location4.x, location4.y));
model4.FaceLoad(2) = faceLoad(Pressure = PV_ratio * V4);
model4.FaceBC([3,4]) = faceBC(Constraint="fixed");
results4 = solve(model4);
u4 = results4.Displacement;
figure;
pdeplot3D(model4.Mesh,'ColorMapData', u4.Magnitude, 'FaceAlpha', 0.3, ...
    'Deformation', u4, ...
    'DeformationScaleFactor', 1);
title('Displacement colormap')
%}
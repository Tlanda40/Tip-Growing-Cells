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
function y = ext(beta, x)
    C1=abs(beta(1));
    C2=abs(beta(2));
    C3=abs(beta(3));
    a1=beta(4);
    a2=beta(5);
    a3=beta(6);
    y=C1*exp(-(x/a1).^2)+C2*exp(-(x/a2).^2)+C3*exp(-(x/a3).^2);
    %y=C1*exp(-(x/a1).^2)+C2*exp(-(x/a2).^2);
    %y=C1*exp(-(x/a1).^2);
end

% Helper function to define stiffness
function y = youngMod(Ext)
    y = 15000 ./ Ext;
    %y = (8.003 - Ext)* 2e5;
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

% Helper function takes Mesh to Graph
function G = meshToSparseGraph(V, T)
    % Extract edges from Mx4 tetrahedra matrix (each row has 4 indices of V)
    E = [T(:, [1, 2]); T(:, [2, 3]); T(:, [3, 4]); T(:, [4, 1]); T(:, [1, 3]); T(:, [2, 4])];  % 6Mx2, each row is pair of vertex indices
    E = sort(E, 2);  % Sorts each row (edge)
    E = unique(E, 'rows');  % Remove duplicates

    % Compute edge lengths
    edgeLengths = vecnorm(V(E(:,1),:) - V(E(:,2),:), 2, 2);

    % Create sparse adjacency graph
    G = graph(E(:,1), E(:,2), edgeLengths, size(V,1));
end

%insert structural properties.

nu = 0.3;           % Poisson's ratio
p = 100;            % Pressure
p1 = -4e3;
%beta0 = [0.002,0.13875,0.0025,0.0025]'; % Inputs for extensibility profile
% beta0 = [0.00375,0.65,0.005,0.002]'; % Inputs for extensibility profile
%beta0 = [0.002,8,0.0075,0.0005]';
%beta0 = [0.65,0.002]';
beta0 = [0.002,0.648,2,0.0075,0.002,0.0007]';

warning('off', 'MATLAB:alphaShape:DupPointsBasicWarnId');
timeSteps = 2;
%frames(timeSteps) = struct('cdata', [], 'colormap', []);
thicknesses = zeros(1, timeSteps);
correctionTimes = [];
heights = zeros(1, timeSteps);
volumes = zeros(1, timeSteps);
v = VideoWriter('my_animation.avi');
v.FrameRate = 15;
v.Quality = 100;
open(v);
figure;
i = 1;
thicknessTick = 0;
while i <= timeSteps
    fprintf('Time Step %d\n', i-1);

    if i > 1
        % Initialize femodel with geometry from previous iter
        model = femodel("AnalysisType","structuralStatic");
        model.Geometry = fegeometry(displacedNodes', elements');
    else
        % Initialize femodel with initial geometry
        model = femodel("AnalysisType","structuralStatic","Geometry", "HollowHemisphere.step");
    end
    %{
    % Plot geometry
    figure;
    pdegplot(model.Geometry, FaceLabels='on', FaceAlpha=0.5);
    axis equal;
    title(['geometry ', num2str(i)])
    %}

    % Generate tetrahedral mesh, can specify resolution
    %model = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.00014, 'GeometricOrder', 'linear');
    %model = generateMesh(model, 'Hmax', 0.0003, 'Hmin', 0.00028, 'GeometricOrder', 'linear');
    model = generateMesh(model, 'Hmax', 0.0003, 'GeometricOrder', 'linear');
    mesh = model.Mesh;
    %{
    % Plot mesh
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
    % Plot alphashape
    figure;
    h = plot(shp);
    h.FaceAlpha = 0.3;
    %}
    V = volume(shp);
    volumes(1,i) = V;

    shpinner = alphaShape(mesh.Nodes');
    shpinner.Alpha = criticalAlpha(shpinner, 'all-points');
    Vinner = volume(shpinner);
    
    if i == 1
        PV_ratio = p / V;
        PV_ratio1init = p1 / Vinner;
    end
    
    % max z pole finder (worse)
    %[~, poleIdx] = max(mesh.Nodes(3,:));

    % projected pole finder (better)
    tempOuterCoords = mesh.Nodes([1, 2],faceNodeIdx{facesSorted(1)});
    poleIdxOuter = knnsearch(tempOuterCoords',[0,0]);
    poleIdx = faceNodeIdx{facesSorted(1)}(1,poleIdxOuter);

    poleCoords = mesh.Nodes(:,poleIdx);
    heights(1,i) = poleCoords(3,1);
    %fprintf("Pole index is %d, pole coords are [%d, %d, %d]\n", poleIdx, poleCoords(1, 1), poleCoords(2, 1), poleCoords(3, 1));

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
    thicknesses(1,i) = thickness;
    criticalThickness = 7.5e-05;
    if thickness < criticalThickness
        if thicknessTick ~= 0
            PV_ratio1 = PV_ratio1 * 2;
        else
            PV_ratio1 = PV_ratio1init;
        end
        model.FaceLoad(facesSorted(1)) = faceLoad(Pressure = PV_ratio1 * Vinner); 
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
    
    % Compute lengths of geodesics. s stores arclength.
    G = meshToSparseGraph(mesh.Nodes', mesh.Elements');
    s = distances(G, poleIdx);  % s(i) is geodesic Dijkstra distance to vertex i

    %{
    % Faces for trisurf visualization
    F = tetToFace(mesh.Elements)';
    F = sort(F, 2); % Sorts each row (face)
    F = unique(F, 'rows'); % Remove duplicates

    % Visualize arclength by color
    figure;
    trisurf(F, mesh.Nodes(1,:), mesh.Nodes(2,:), mesh.Nodes(3,:), s, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.4);
    colormap(turbo); colorbar;
    title('Dijkstra Geodesic Distance from Pole');
    axis equal;
    %}

    model.MaterialProperties.PoissonsRatio = nu;

    kdtree = createns(mesh.Nodes', 'NSMethod', 'kdtree');
    val = @(location,state) youngMod(ext(beta0, s(knnsearch(kdtree, ...
    [location.x(:), location.y(:), location.z(:)]))));

    model.MaterialProperties.YoungsModulus = val;

    results = solve(model);
    u = results.Displacement;
    %figure;
    clf;
    pdeplot3D(model.Mesh,'ColorMapData', u.Magnitude, 'FaceAlpha', 0.3, ...
        'Deformation', u, ...
        'DeformationScaleFactor', 1);
    clim([0 4e-5]);
    axis([-0.011 0.011 -0.011 0.011 -0.001 0.025]);
    view(90, 0);
    ax.CameraPositiion = [0.8, 0, 0.0125];
    %camtarget([0, 0, 0.0125]);
    if thickness < criticalThickness
        title(['Loop displacement colormap ', num2str(i-1)]);
        correctionTimes = [correctionTimes, i-1];
    else
        title(['Loop displacement colormap ', num2str(i)]);
    end
    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
    displacedNodes = model.Mesh.Nodes + [u.ux'; u.uy'; u.uz'];
    elements = model.Mesh.Elements;
    if mvForward == 1
        i = i + 1;
    end
end
warning('on', 'MATLAB:alphaShape:DupPointsBasicWarnId');
close(v);
ffmpegPath = '/opt/homebrew/bin/ffmpeg';
inFile = 'my_animation.avi';
outFile = 'my_animation.mp4';
delete('my_animation.mp4');
cmd = sprintf('"%s" -i "%s" -vcodec libx264 -crf 18 "%s"', ffmpegPath, inFile, outFile);
system(cmd);
delete('my_animation.avi');


    % Plot alphashape
    figure;
    h = plot(shp);
    h.FaceAlpha = 0.3;
    

time = 0:1:timeSteps-1;
figure;
plot(time,thicknesses);
hold on;
for i = 1:size(correctionTimes,2)
    xline(correctionTimes(i), 'r--');
end
hold off;
title('Thickness over time');
figure;
plot(time,heights);
hold on;
for i = 1:size(correctionTimes,2)
    xline(correctionTimes(i), 'r--');
end
hold off;
title('Height over time');
figure;
dheights = zeros(1,timeSteps);
for i = 2:timeSteps
    dheights(1,i) = heights(1,i) - heights(1,i-1);
end
plot(time,dheights);
hold on;
for i = 1:size(correctionTimes,2)
    xline(correctionTimes(i), 'r--');
end
hold off;
title('Change in Height over time');
figure;
plot(time,volumes);
hold on;
for i = 1:size(correctionTimes,2)
    xline(correctionTimes(i), 'r--');
end
hold off;
title('Volume over time');
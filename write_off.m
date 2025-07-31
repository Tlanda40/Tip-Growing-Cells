function write_off(filename, V, F)
    
    % Check inputs
    if size(V, 2) ~= 3
        error("V (vertices) must be an Nx3 array of 3D coords");
    end

    if size(F, 2) ~= 3
        error("F (faces) must be an Nx3 array of indices into V");
    end

    if ~endsWith(filename, '.off')
        filename = [filename, '.off'];
    end

    % Open file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error("Could not open file %s for writing.", filename);
    end

    % Write header ('OFF' and number of verts, edges, faces)
    fprintf(fid, "OFF\n");
    fprintf(fid, "%d %d 0\n", size(V,1), size(F,1));

    % Write vertices
    fprintf(fid, "%f %f %f\n", V');

    % Write faces
    fprintf(fid, "3 %d %d %d\n", (F - 1)');

    fclose(fid);
end
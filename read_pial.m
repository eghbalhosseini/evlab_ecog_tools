function mesh = read_pial(filename)
    fid = fopen(filename,'r');
    if fid<0
        error('Could not open %s', filename);
    end
    % First line: <nVertices> <nFaces> <unused>
    hdr = fscanf(fid, '%d %d %d\n', 3);
    nV = hdr(1);
    nF = hdr(2);
    % Next nV lines: x y z coordinates
    verts = fscanf(fid, '%f %f %f\n', [3, nV])';
    % Next nF lines: v1 v2 v3 indices (0-based)
    faces0 = fscanf(fid, '%d %d %d\n', [3, nF])';
    fclose(fid);
    % Convert to 1-based indexing for MATLAB
    mesh.vertices = verts;
    mesh.faces    = faces0 + 1;
end
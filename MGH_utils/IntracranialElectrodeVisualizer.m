classdef IntracranialElectrodeVisualizer < handle
    % IntracranialElectrodeVisualizer – ECoG visualization using FreeSurfer
    % Updated to work with organized electrode data structure
    
    properties
        freesurferHome char
        subjectsDir    char
        reconDir       char
        elecDataDir    char
        surfaces       struct
        annotations    struct
        electrodes     table
        config         struct
        figureHandles  = []
        axesHandles    = []
    end
    
    methods
        function obj = IntracranialElectrodeVisualizer(opts)
            arguments
                opts.freesurferHome char = getenv('FREESURFER_HOME')
                opts.subjectsDir    char = getenv('SUBJECTS_DIR')
                opts.reconDir       char = getenv('RECONDIR')
                opts.elecDataDir    char = ''
                opts.config         struct = struct()
            end
            obj.freesurferHome = opts.freesurferHome;
            obj.subjectsDir    = opts.subjectsDir;
            obj.reconDir       = opts.reconDir;
            obj.elecDataDir    = opts.elecDataDir;
            obj.config         = obj.mergeStructs(obj.getDefaultConfig(), opts.config);
            obj.verifyEnvironment();
        end
        
        function loadSubject(obj, subjectName, cfg)
            arguments
                obj
                subjectName char
                cfg.surfaceTypes cell = {'pial','inflated','white'}
                cfg.hemispheres  cell = {'lh','rh'}
            end
            opts = obj.mergeStructs(obj.config, cfg);
            obj.config.subject = subjectName;
            
            % Ensure top-level subject field exists in obj.surfaces
            obj.surfaces = obj.ensureField(obj.surfaces, subjectName);
            
            for h = opts.hemispheres
                hemi = h{1};
                % Ensure hemisphere field exists for this subject
                obj.surfaces.(subjectName) = obj.ensureField(obj.surfaces.(subjectName), hemi);
                
                % For each surface type, load surfaces from file
                for st = opts.surfaceTypes
                    surfType = st{1};
                    surfFile = fullfile(obj.subjectsDir, subjectName, 'surf', sprintf('%s.%s', hemi, surfType));
                    [v, f] = obj.read_surf(surfFile);
                    obj.surfaces.(subjectName).(hemi).(surfType) = struct('vertices', v, 'faces', f);
                end
            end
        end

        
        function loadElectrodes(obj, subjectList, cfg)
            arguments
                obj
                subjectList cell = {}
                cfg.coordType char {mustBeMember(cfg.coordType,{'coords','mni_coords'})} = 'coords'
                cfg.file  char = ''
                cfg.table table = table.empty()
            end
            opts = obj.mergeStructs(obj.config, cfg);
            if ~isempty(opts.table)
                obj.electrodes = opts.table;
            elseif ~isempty(opts.file)
                obj.electrodes = readtable(opts.file,'Delimiter','\t');
            elseif ~isempty(subjectList)
                obj.electrodes = obj.loadFromOrganizedData(subjectList, opts.coordType);
            else
                obj.electrodes = table();
            end
        end
        
        function elecTable = loadFromOrganizedData(obj, subjectList, coordType)
            arguments
                obj
                subjectList cell
                coordType char = 'coords'
            end
            allCoords = []; allNames = {}; allSubjects = {}; allLabels = {};
            for subID = subjectList(:)'
                matFile = fullfile(obj.elecDataDir, subID{1}, 'electrodes_combined.mat');
                if ~isfile(matFile)
                    warning('No combined data file for subject %s', subID{1});
                    continue;
                end
                C = load(matFile,'combined').combined;
                coords = obj.pickField(C, coordType, 'coords');
                names  = obj.pickField(C, 'names', arrayfun(@(i) sprintf('%s_%d',subID{1},i),1:size(coords,1),'UniformOutput',false)');
                labels = obj.pickField(C, 'labels', repmat({''}, size(coords,1),1));
                
                allCoords    = [allCoords; coords];
                allNames     = [allNames;  names(:)];
                allSubjects  = [allSubjects; repmat(subID, size(coords,1),1)];
                allLabels    = [allLabels;  labels(:)];
                fprintf('Loaded %d electrodes for subject %s\n', size(coords,1), subID{1});
            end
            elecTable = table(allCoords(:,1), allCoords(:,2), allCoords(:,3), ...
                allNames, allSubjects, allLabels, ...
                'VariableNames', {'x','y','z','name','subject','anatomical_label'});
        end
        
        function plotElectrodeDensity(obj, subjectList, sigma_mm, surfaceName, cfg)
            arguments
                obj
                subjectList cell
                sigma_mm double = 2
                surfaceName char = 'fsaverage'
                cfg.surfaceType char {mustBeMember(cfg.surfaceType, {'pial','inflated'})} = 'pial'
                cfg.threshold double = 0.1
                cfg.colormap = flipud(hot(256));
                cfg.showElectrodes logical = false
            end
            opts = obj.mergeStructs(obj.config, cfg);
            % Load pial and display surfaces
            surfPath = fullfile(obj.subjectsDir, surfaceName, 'surf');
            pial = obj.readBothHemispheres(surfPath, 'pial');
            displaySurf = obj.readBothHemispheres(surfPath, opts.surfaceType);
            curv = struct(...
                'lh', read_curv(fullfile(surfPath,'lh.curv')), ...
                'rh', read_curv(fullfile(surfPath,'rh.curv')));
            
            coords = obj.aggregateCoords(subjectList, 'coords');
            if isempty(coords), error('No electrode coordinates loaded.'); end
            
           D_pial.lh = obj.computeElectrodeDensity( ...
                 pial.lh.vertices, coords, sigma_mm);
            D_pial.rh = obj.computeElectrodeDensity( ...
                 pial.rh.vertices, coords, sigma_mm);

            [D.lh, D.rh] = obj.normalizeAndProject(D_pial, pial, displaySurf, opts.surfaceType);
            
            obj.plotDensityHemisphere('lh', displaySurf.lh, curv.lh, coords, D.lh, opts);
            obj.plotDensityHemisphere('rh', displaySurf.rh, curv.rh, coords, D.rh, opts);
        end
        
        function plotElectrodesByActivity(obj, subjectList, surfaceName, activityType, cfg)
            arguments
                obj
                subjectList cell
                surfaceName char = 'fsaverage'
                activityType char {mustBeMember(activityType,{'ictal','interictal'})} = 'ictal'
                cfg.surfaceType char {mustBeMember(cfg.surfaceType,{'pial','inflated'})} = 'inflated'
                cfg.colormap char = 'hot'
            end
            opts = obj.mergeStructs(obj.config, cfg);
            surfPath = fullfile(obj.subjectsDir, surfaceName, 'surf');
            hemiData = struct();
            for h = ["lh","rh"]
                [verts,faces] = obj.read_surf(fullfile(surfPath, h + "." + opts.surfaceType));
                hemiData.(h).verts = verts; hemiData.(h).faces = faces;
                hemiData.(h).curv  = read_curv(fullfile(surfPath, h + ".curv"));
            end
            
            [coords, act] = obj.aggregateCoordsActivity(subjectList, activityType);
            fig = figure('Color','w','Position',[100 100 1200 600]);
            hemiList = ["lh", "rh"];
            for idx = 1:2
                hemi = hemiList(idx);
                ax = subplot(1,2,idx); hold(ax,'on');
                obj.plotGraySurface(ax, hemiData.(hemi).verts, hemiData.(hemi).faces, hemiData.(hemi).curv);
                sel = (coords(:,1)<0 & hemi=="lh") | (coords(:,1)>0 & hemi=="rh");
                obj.scatterActivity(ax, coords(sel,:), act(sel), opts.colormap);
                title(ax, upper(hemi),"FontSize",14,"FontWeight","bold");
            end
            legend({'Inactive','Active'},'Position',[0.45 0.05 0.1 0.05]);
        end
        
        function loadAtlas(obj, subjectName, atlasName, cfg)
            arguments
                obj
                subjectName char
                atlasName   char
                cfg.hemispheres cell = {'lh','rh'}
            end
            opts = obj.mergeStructs(obj.config, cfg);
            obj.config.subject = subjectName;
            
            % Ensure the top-level field exists
            if isempty(fieldnames(obj.annotations))
                obj.annotations = struct(subjectName, struct());
            elseif ~isfield(obj.annotations, subjectName)
                obj.annotations.(subjectName) = struct();
            end
            
            for h = opts.hemispheres
                hemi = h{1};
                % Ensure the hemisphere field exists before adding atlas
                if ~isfield(obj.annotations.(subjectName), hemi) || isempty(obj.annotations.(subjectName).(hemi))
                    obj.annotations.(subjectName).(hemi) = struct();
                end
                
                % Now safe to assign deeply
                fn = fullfile(obj.subjectsDir, subjectName, 'label', sprintf('%s.%s.annot', hemi, atlasName));
                obj.annotations.(subjectName).(hemi).(atlasName) = obj.readAnnotation(fn);
            end
        end

        
        function plotSurface(obj, subjectName, surfaceType, hemisphere, cfg)
            arguments
                obj
                subjectName char
                surfaceType char {mustBeMember(surfaceType,{'pial','inflated','white'})}
                hemisphere  char {mustBeMember(hemisphere,{'lh','rh','both'})}
                cfg          struct = struct()
            end
            opts = obj.mergeStructs(obj.mergeStructs(obj.config,cfg), struct('subject',subjectName,'surfaceType',surfaceType,'hemisphere',hemisphere));
           if strcmpi(opts.hemisphere, 'both')
                hems = ["lh", "rh"];
            else
                hems = string(opts.hemisphere);
            end

            for hemi = hems
                surf = obj.surfaces.(opts.subject).(hemi).(opts.surfaceType);
                fig = obj.createSurfaceFigure(surf.vertices, surf.faces, opts, hemi);
                obj.figureHandles(end+1) = fig;
            end
        end
        
        function plotElectrodes(obj, colorVar, cutoff, hemisphere, cfg)
            arguments
                obj
                colorVar   double = ones(height(obj.electrodes),1)
                cutoff     double = 0.2
                hemisphere char {mustBeMember(hemisphere,{'lh','rh','both'})} = 'both'
                cfg        struct = struct()
            end
            opts = obj.mergeStructs(obj.config, cfg);
            cmap = obj.createLanguageColormap();
            if strcmpi(opts.hemisphere, 'both')
                hems = ["lh", "rh"];
            else
                hems = string(opts.hemisphere);
            end

            for hemi = hems
                fig = findobj(obj.figureHandles,'Tag',hemi);
                if isempty(fig), error("No surface figure tagged '%s' found.",hemi); end
                ax = findall(fig(1),'Type','axes'); hold(ax,'on');
                sel = (hemi=="lh" & obj.electrodes.x<0) | (hemi=="rh" & obj.electrodes.x>0);
                coords = [obj.electrodes.x(sel),obj.electrodes.y(sel),obj.electrodes.z(sel)];
                vals   = colorVar(sel);
                high = vals>=cutoff;
                scatter3(ax,coords(high,1),coords(high,2),coords(high,3),opts.elec_size,vals(high),'filled');
                colormap(ax,cmap); caxis(ax,[cutoff 1]);
                scatter3(ax,coords(~high,1),coords(~high,2),coords(~high,3),opts.elec_size/2,[0.7 0.7 0.7],'filled');
            end
        end
        
        function plotAnnotation(obj, hemi, atlasName, regionNames, cfg)
            arguments
                obj
                hemi       char {mustBeMember(hemi,{'lh','rh','both'})}
                atlasName  char
                regionNames  cell  % cell array of region name strings
                cfg        struct = struct()
            end
            
            opts = obj.mergeStructs(obj.config, cfg);
            subj = opts.subject;
            if strcmpi(hemi,'both')
                hemis = {'lh','rh'};
            else
                hemis = {hemi};
            end
            
            for hIdx = 1:length(hemis)
                h = hemis{hIdx};
                fh = findobj(obj.figureHandles, 'Tag', h);
                if isempty(fh)
                    error('No surface figure tagged ''%s'' found.', h);
                end
                figure(fh(1));
                ax = findall(fh(1), 'Type', 'axes');
                hold(ax, 'on');
                annotStruct = obj.annotations.(subj).(h).(atlasName);
                
                % Map region names to IDs
                labelIDs = [];
                for rNameCell = regionNames
                    rName = rNameCell{1};
                    % Find label matching the region name
                    matchIdx = find(strcmpi(annotStruct.colortable.struct_names, rName));
                    if ~isempty(matchIdx)
                        labelIDs(end+1) = matchIdx; % IDs are probably 0-based in labels
                    else
                        warning('Region name "%s" not found in atlas labels for %s', rName, h);
                    end
                end
                patches = findall(ax,'Type','Patch');
                obj.applyAnnotationPatch(patches, annotStruct, labelIDs);
            end
        end
    
        function D = computeElectrodeDensity(~, vertices, electrodes, sigma_mm)
            D = zeros(size(vertices,1),1);
            for v=1:size(vertices,1)
                d2 = sum((electrodes - vertices(v,:)).^2,2);
                D(v) = sum(exp(-d2/(2*sigma_mm^2)));
            end
        end

        function D = computeElectrodeDensityGeodesic(obj, vertices, faces, electrodes, sigma_mm)
        % Compute density via geodesic distances on the mesh
        arguments
            obj
            vertices double
            faces    double
            electrodes double
            sigma_mm double = 2
        end
        
        N = size(vertices,1);
        D = zeros(N,1);
        
        % Build adjacency graph once
        adj = obj.buildSurfaceAdjacency(vertices, faces);
        cutoff = 3 * sigma_mm;  % limit distances for speed
        
        for e = 1:size(electrodes,1)
            % Find nearest surface vertex in Euclidean space
            [~, src] = min( vecnorm(vertices - electrodes(e,:), 2, 2) );
            
            % Compute geodesic distances truncated at cutoff
            gdist = obj.dijkstraGeodesic(adj, src, cutoff);
            
            % Gaussian kernel contribution
            mask = gdist <= cutoff;
            D(mask) = D(mask) + exp(-(gdist(mask).^2)/(2*sigma_mm^2));
        end
        
        % Normalize to [0,1]
        if max(D) > 0
            D = D / max(D);
        end
    end
    
    end
    
    methods (Access=private)
        function S = loadHemisphereSurfaces(obj, subj, hemi, types)
            S = struct();
            for t = types
                surfType = t{1};
                fn = fullfile(obj.subjectsDir, subj, 'surf', sprintf('%s.%s', hemi, surfType));
                [v,f] = obj.read_surf(fn);
                S.(surfType) = struct('vertices',v,'faces',f);
            end
        end
        
        function data = readBothHemispheres(obj, surfPath, surfType)
            data = [];
            [vertices,faces] = obj.read_surf(fullfile(surfPath,strcat('lh.',surfType)));
            data.lh.vertices = vertices;
            data.lh.faces = faces;
           [vertices,faces] = obj.read_surf(fullfile(surfPath,strcat('rh.',surfType)));
            data.rh.vertices = vertices;
            data.rh.faces = faces;
        end
        
        function coords = aggregateCoords(obj, subjectList, fieldName)
            coords = [];
            for s = subjectList(:)'
                mf = fullfile(obj.elecDataDir, s{1}, 'electrodes_combined.mat');
                if isfile(mf)
                    C = load(mf,'combined').combined;
                    coords = [coords; obj.pickField(C, fieldName, [])];
                end
            end
        end
        
        function [coords, act] = aggregateCoordsActivity(obj, subjectList, activityType)
            coords = []; act = [];
            for s = subjectList(:)'
                mf = fullfile(obj.elecDataDir, s{1}, 'electrodes_combined.mat');
                if isfile(mf)
                    C = load(mf,'combined').combined;
                    coords = [coords; obj.pickField(C,'mni_coords',obj.pickField(C,'coords',[]))];
                    act = [act; obj.pickField(C, activityType+"_activity", zeros(size(coords,1),1))];
                end
            end
        end
        
        function val = pickField(~, structIn, fld, defaultVal)
            if isfield(structIn, fld)
                val = structIn.(fld);
            else
                val = defaultVal;
            end
        end
        
        

        
        function [D_lh_norm,D_rh_norm] = normalizeAndProject(obj, D_pial, pial, dispSurf, surfType)
            maxD = max([D_pial.lh; D_pial.rh]);
            D_pial.lh = D_pial.lh./maxD; D_pial.rh = D_pial.rh./maxD;
            if strcmp(surfType,'inflated')
                D_lh_norm = obj.projectDensity(D_pial.lh,pial.lh.vertices,dispSurf.lh.vertices);
                D_rh_norm = obj.projectDensity(D_pial.rh,pial.rh.vertices,dispSurf.rh.vertices);
            else
                D_lh_norm = D_pial.lh; D_rh_norm = D_pial.rh;
            end
        end
        
        function D_inf = projectDensity(~, D_pial, verts_pial, verts_inf)
            idx = knnsearch(createns(verts_pial,'NSMethod','kdtree'), verts_inf);
            D_inf = D_pial(idx);
        end
        
        function plotDensityHemisphere(obj, hemi, surf, curv, coords, D, opts)
            fig = figure('Name', upper(hemi)+" Density",'Color','w','Position',[100,100,600,600]);
            ax = axes('Parent',fig); hold(ax,'on');
            Cgray = obj.create_gray_curvature_colors(curv);
            patch('Vertices',surf.vertices,'Faces',surf.faces,'FaceVertexCData',Cgray, ...
                  'FaceColor','flat','EdgeColor','none','Parent',ax);
            Dm = D; alpha = Dm>=opts.threshold; Dm(~alpha)=NaN;
            patch('Vertices',surf.vertices,'Faces',surf.faces,'FaceVertexCData',Dm, ...
                  'FaceColor','interp','EdgeColor','none','FaceAlpha','flat','FaceVertexAlphaData',double(alpha),'Parent',ax);
            if opts.showElectrodes
                side = (coords(:,1) < 0 & strcmp(hemi, 'lh')) | (coords(:,1) > 0 & strcmp(hemi, 'rh'));


                scatter3(ax,coords(side,1),coords(side,2),coords(side,3),20,'k','filled');
            end
            obj.applyLighting(ax, hemi);
            axis(ax,'equal','off');
            xlim(ax, [min(surf.vertices(:,1)), max(surf.vertices(:,1))]);
            ylim(ax, [min(surf.vertices(:,2)), max(surf.vertices(:,2))]);
            zlim(ax, [min(surf.vertices(:,3)), max(surf.vertices(:,3))]);
            colormap(ax,opts.colormap); caxis(ax,[opts.threshold 1]);
            colorbar('peer',ax,'Location','eastoutside');
            title(ax,sprintf('%s Hemisphere %s Density',upper(hemi),opts.surfaceType),'FontSize',14);
        end
        
        function fig = createSurfaceFigure(obj, verts, faces, opts, hemi)
            if min(faces(:))==0, faces=faces+1; end
            fig = figure('Color',opts.backgroundColor,'Position',opts.figurePosition,'Tag',hemi); 
            ax = axes('Parent',fig); hold(ax,'on');
            if strcmp(opts.surfaceType,'inflated')
                sulc = read_curv(fullfile(obj.subjectsDir,opts.subject,'surf',sprintf('%s.curv',hemi)));
                C = repmat(opts.gyrusColor,numel(sulc),1);
                mask = sulc<opts.sulcThreshold;
                C(mask,:) = repmat(opts.sulcusColor,sum(mask),1);
                fc = 'flat';
            else
                C = opts.surfColor; fc = 'flat';
            end
            patch('Vertices',verts,'Faces',faces,'FaceVertexCData',C,'FaceColor',fc, ...
                  'EdgeColor','none','FaceAlpha',opts.alpha,'Parent',ax);
            axis(ax,'equal','off');
            xlim(ax, [min(verts(:,1)), max(verts(:,1))]);
            ylim(ax, [min(verts(:,2)), max(verts(:,2))]);
            zlim(ax, [min(verts(:,3)), max(verts(:,3))]);

            lighting(ax,'phong'); material(ax,'dull');
            obj.applyLighting(ax,hemi);
            title(ax,sprintf('%s %s Surface',upper(hemi),opts.surfaceType));
            rotate3d(fig,'on');
        end
        
        function annot = readAnnotation(~,fname)
            if endsWith(fname,'.mat')
                annot = load(fname,'annot').annot;
            else
                [~,labels,ct] = read_annotation(fname);
                annot.label      = labels;
                annot.colortable = ct;
            end
        end
        
        function S = ensureField(~, S, fld)
            if isempty(S)
                % Create a new struct with the desired field
                S = struct(fld, struct());
            elseif ~isfield(S, fld)
                S.(fld) = struct();
            end
        end

        
        function cfg = getDefaultConfig(~)
            cfg = struct( ...
                'subject','fsaverage','surfaceType','pial','hemisphere','both', ...
                'figurePosition',[100 100 1200 800], 'backgroundColor',[1 1 1], ...
                'surfColor',[0.8 0.8 0.8],'alpha',.6, ...
                'sulcThreshold',0,'gyrusColor',[.9 .9 .9],'sulcusColor',[.6 .6 .6], ...
                'elec_size',40 ...
            );
        end
        
        function out = mergeStructs(~, a, b)
            out = a;
            for f = string(fieldnames(b))'
                out.(f) = b.(f);
            end
        end
        
        function verifyEnvironment(obj)
            assert(isfolder(obj.freesurferHome),'Invalid FREESURFER_HOME');
            assert(isfolder(obj.subjectsDir),   'Invalid SUBJECTS_DIR');
            if ~isempty(obj.elecDataDir)
                assert(isfolder(obj.elecDataDir), 'Invalid electrode data directory');
            end
        end
        
        function applyLighting(~, ax, hemi)
            lighting(ax,'phong'); material(ax,'dull'); axis(ax,'equal','off');
            if hemi=='lh'
                view(ax,[-90 0]); camlight(ax,0,0);
            else
                view(ax,[90 0]); camlight(ax,90,0);
            end
        end
        
        function applyAnnotationPatch(~, hSurf, annot, regIDs)
            % Set colors for the specified regions
            V = numel(annot.label);
            C = repmat([0.7 0.7 0.7], V, 1);  % default gray
            
            ct = annot.colortable.table(:, 1:3)/255;
            
            for rID = regIDs
                mask = (annot.label == annot.colortable.table(rID, 5));
                if any(mask)
                    C(mask, :) = repmat(ct(rID+1, :), sum(mask), 1);
                end
            end
        
            set(hSurf, 'FaceVertexCData', C, 'FaceColor', 'flat');
        end
        function applyAnnotation(~, hSurf, annot, regions, opts)
            ct = annot.colortable.table(:,1:3) / 255;  % Color table normalized
            labels = annot.label;
            V = numel(labels);
            C = repmat([0.7 0.7 0.7], V, 1);  % Default gray color
            
            maxID = size(ct,1) - 1;
            valid = regions(regions >= 1 & regions <= maxID);
            
            for r = valid(:)'
                mask = labels == r;
                C(mask,:) = repmat(ct(r+1,:), sum(mask), 1);  % assign region colors
            end
            
            set(hSurf, 'FaceVertexCData', C, 'FaceColor', 'flat');
        end

        function C = create_gray_curvature_colors(~, curv)
            g = (curv-min(curv))/(max(curv)-min(curv));
            C = repmat(g,1,3);
        end
        
        function [v, f] = read_surf(~, fname)
        if ~exist(fname, 'file'), error('Missing %s', fname); end
        
        % Call FreeSurfer's read_surf function
        [v,f] = freesurfer_read_surf(fname);
        
        % If cell type return, unwrap (depends on freesurfer_read_surf implementation)
        if iscell(v), v = v{1}; end
        if iscell(f), f = f{1}; end
        
        % Convert to double if needed
        v = double(v);
        f = double(f);
        
        % MATLAB indexing: use 1-based indexing
        if min(f(:)) == 0, f = f + 1; end
    end

        function cmap = createLanguageColormap(~)
            n = 256; x = linspace(0,1,n)';
            cmap = zeros(n,3);
            for k=1:n
                if x(k)<.25
                    cmap(k,:) = [1 1 4*x(k)];
                elseif x(k)<.6
                    cmap(k,:) = [1 1-(x(k)-.25)/.35 0];
                else
                    t=(x(k)-.6)/.4; cmap(k,:)=[1-t 0 0];
                end
            end
        end
        function adj = buildSurfaceAdjacency(~, vertices, faces)
            % Ensure faces are 1-based
            if min(faces(:)) == 0
                faces = faces + 1;
            end
            
            % Build a sparse adjacency matrix of edge lengths
            edges = [faces(:,[1 2]); faces(:,[2 3]); faces(:,[3 1])];
            edges = sort(edges,2);  % ensure (i<j)
            [E, ~, ~] = unique(edges, 'rows');
            i = E(:,1); j = E(:,2);
            
            v1 = vertices(i,:); 
            v2 = vertices(j,:);
            w  = sqrt(sum((v1 - v2).^2, 2));
            
            N = size(vertices,1);
            adj = sparse([i; j], [j; i], [w; w], N, N);
        end


    function dist = dijkstraGeodesic(~, adj, src, maxd)
        % Compute geodesic distances from a source vertex up to maxd
        %
        % adj:  sparse N×N adjacency matrix of edge lengths
        % src:  source vertex index (scalar)
        % maxd: maximum distance to compute (truncate)
        
        N = size(adj,1);
        dist = inf(N,1);
        dist(src) = 0;
        
        visited = false(N,1);
        Q = 1:N;  % naive priority queue
        
        while ~isempty(Q)
            % Extract unvisited node with minimal dist
            [~, idx] = min(dist(Q));
            u = Q(idx);
            if dist(u) > maxd
                break;
            end
            visited(u) = true;
            Q(idx) = [];  % remove u
            
            % Relax neighbors
            nbrs = find(adj(u,:));
            for v = nbrs
                alt = dist(u) + adj(u,v);
                if alt < dist(v) && alt <= maxd
                    dist(v) = alt;
                end
            end
        end
    end
    
    end

    
end

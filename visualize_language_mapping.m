function visualize_language_mapping(results, lh, rh, cutoff, figsize)
    arguments
        results table
        lh struct
        rh struct
        cutoff double {mustBePositive, mustBeLessThan(cutoff,1)} = 0.2
        figsize double {mustBeVector, mustBePositive} = [900, 600]
    end
    
    mycmap = lanamap();
    for h = 1:2
        figure('Color','w','Position',[100 100 figsize]);
        ax = axes('Position',[0 0 1 1],'Visible','off');
        hold(ax,'on');
        
        if h == 1
            mesh = lh; view(ax, [-90 0]); idx = results.R < 0;
            title(ax,'Left Hemisphere');
        else
            mesh = rh; view(ax, [90 0]); idx = results.R > 0;
            title(ax,'Right Hemisphere');
        end
        
        plot_surface(mesh, ax);
        plot_electrodes(results, idx, cutoff, ax);
        
        colormap(ax, mycmap);
        caxis(ax, [cutoff 0.8]);
        add_colorbar(mycmap, cutoff, ax);
    end
end

function c = lanamap()
    nColors = 256;
    c = zeros(nColors,3);
    for k = 1:nColors
        x = (k-1)/(nColors-1);
        if x < 0.25
            c(k,:) = [1, 1, x*4];
        elseif x < 0.6
            c(k,:) = [1, 1-(x-0.25)/0.35, 0];
        else
            t = (x-0.6)/0.4;
            c(k,:) = [1-t, 0, 0];
        end
    end
end

function plot_surface(mesh, ax)
    faces = mesh.faces + 1;
    trisurf(faces, mesh.vertices(:,1), mesh.vertices(:,2), mesh.vertices(:,3), ...
        'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'Parent', ax);
    lighting(ax, 'gouraud');
    camlight(ax, 'headlight');
    axis(ax, 'equal', 'off', 'tight');
end

function plot_electrodes(results, idx, cutoff, ax)
    color_idx = idx & (results.pLanA >= cutoff);
    if any(color_idx)
        scatter3(results.R(color_idx), results.A(color_idx), results.S(color_idx), ...
            40, results.pLanA(color_idx), 'filled', 'Parent', ax);
    end
    gray_idx = idx & (results.pLanA < cutoff);
    if any(gray_idx)
        scatter3(results.R(gray_idx), results.A(gray_idx), results.S(gray_idx), ...
            8, [0 0 0], 'filled', 'Parent', ax);
    end
end

function add_colorbar(cmap, cutoff, ax)
    cb = colorbar(ax, 'Position', [0.92 0.3 0.02 0.4]);
    colormap(ax, cmap);
    caxis(ax, [cutoff 0.8]);
    cb.Label.String = 'LanA probability';
    cb.Ticks = linspace(cutoff, 0.8, 5);
end

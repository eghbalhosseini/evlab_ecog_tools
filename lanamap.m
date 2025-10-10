function lmap = lanamap(nColors)
% Custom colormap similar to publication figure

lmap = zeros(nColors,3);

% Define the gradient: white to yellow to red to black
for k = 1:nColors
    x = (k-1)/(nColors-1); 
    if x < 0.25
        % White to yellow
        lmap(k,:) = [1, 1, x*4]; % white to yellow
    elseif x < 0.6
        % Yellow to red
        lmap(k,:) = [1, 1-(x-0.25)/0.35, 0]; % yellow fading to red
    else
        % Red to black
        t = (x-0.6)/0.4;
        lmap(k,:) = [1-t, 0, 0]; % red to black
    end
end

% Clamp values to [0,1]
lmap = max(0, min(1, lmap));

end
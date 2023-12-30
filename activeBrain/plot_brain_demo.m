%% General parameter
param.subject_name           = '046';       
param.electrodes_pos         = 'right';     %right | left
param.kernel                 = 'linear';    %linear | gaussian
param.param                  = 10;
param.cutoff                 = 10;
param.show_electrodes        = 1;           
param.show_electrodes_number = 1;

%% ActiveBrain path
current_path = cd;
addpath([current_path,'/activeBrain']);
addpath([current_path,'/activeBrain/geometry']);
addpath([current_path,'/activeBrain/transformModel']);
clear current_path;


%% Brain model calculation
opengl neverselect
fprintf('> Calculating brain model for patient %s ...\n', param.subject_name);
%Load brain model
load (['./AMC/brain_models/AMC',param.subject_name,'_brain.mat']);

%% Calculate Electrode Contributions
if (~strcmp(param.kernel,'linear') || param.param ~= 10 || param.cutoff ~= 10),
    % model is prepared for linear kernel with param 10 and cutoff 10.  
    [ vcontribs ] = electrodesContributions( cortex, tala, param.kernel, param.param, param.cutoff);
end


%% Display options
viewstruct.what2view    = {'brain','activations'};
if strcmp(param.electrodes_pos,'right')
    viewstruct.viewvect     = [90, 0];
    viewstruct.lightpos     = [150, 0, 0];
elseif strcmp(param.electrodes_pos,'left')
    viewstruct.viewvect     = [270, 0];
    viewstruct.lightpos     = [-150, 0, 0];
end

viewstruct.material     = 'dull';
viewstruct.enablelight  = 1;
viewstruct.enableaxis   = 0;
viewstruct.lightingtype = 'gouraud';

cmapstruct.basecol          = [0.7, 0.7, 0.7];
cmapstruct.fading           = true;
cmapstruct.enablecolormap   = true;
cmapstruct.enablecolorbar   = true;
cmapstruct.color_bar_ticks  = 4;

ix = 1;


%% Brain activation plot
fprintf('\t>Ploting correlations on the brain ...\n');
cmapstruct.cmap = colormap('Jet'); close(gcf); %because colormap creates a figure
cmapstruct.ixg2 = floor(length(cmapstruct.cmap) * 0.15);
cmapstruct.ixg1 = -cmapstruct.ixg2;

%Random activation
tala.activations = rand(length(tala.electrodes),1);

cmapstruct.cmin = min(tala.activations);
cmapstruct.cmax = max(tala.activations);
figure;
set(gca,'visible','on');
%figure('visible', 'on');
set(gcf,'Color','w');

activateBrain(cortex, vcontribs, tala, ix, cmapstruct, viewstruct );

plotBalls(tala.trielectrodes, 'r',1.0);

%plotElNums(tala.trielectrodes*1.1,1:size(tala.trielectrodes,1),12,'k');

colorbar off
%zoom(1.4);




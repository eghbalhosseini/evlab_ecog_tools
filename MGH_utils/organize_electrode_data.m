% organize_electrode_data.m
function organize_electrode_data(source_dir, output_dir)
    if nargin<2, error('Usage: organize_electrode_data(source_dir, output_dir)'); end
    if ~exist(output_dir,'dir'), mkdir(output_dir); end

    % Gather localization files (CSV & TSV)
    loc_dir = fullfile(source_dir,'brainstorm_localizations_cat12_fedorenko_nov_2024');
    loc_list = [dir(fullfile(loc_dir,'*.csv')); dir(fullfile(loc_dir,'*.tsv'))];
    loc_list = loc_list(~startsWith({loc_list.name},'.'));

    % Gather anatomical label files (CSV & TSV)
    label_dir = fullfile(source_dir,'brainstorm_localizations');
    label_list = [dir(fullfile(label_dir,'*.csv')); dir(fullfile(label_dir,'*.tsv'))];
    label_list = label_list(~startsWith({label_list.name},'.'));

    fprintf('Found %d localization files and %d label files\n', numel(loc_list), numel(label_list));
    processed={}; failed={};
    for i=1:numel(loc_list)
        fname=loc_list(i).name; fpath=fullfile(loc_list(i).folder,fname);
        subj=extract_subject_id(fname);
        if isempty(subj), warning('No ID in %s',fname); failed{end+1}=fname; continue; end

        % Convert TSV→CSV for localization
        [~,n,ext]=fileparts(fname);
        if strcmpi(ext,'.tsv')
            Tloc=readtable(fpath,'FileType','text','Delimiter','\t');
            csv_loc=fullfile(loc_list(i).folder,[n '.csv']);
            writetable(Tloc,csv_loc);
            fpath=csv_loc;
        end

        % Find & convert label file if needed
        lbl=find_label_file(label_list,subj);
        lbl_path='';
        if ~isempty(lbl)
            [~,m,le]=fileparts(lbl.name);
            orig=fullfile(lbl.folder,lbl.name);
            if strcmpi(le,'.tsv')
                Tlbl=readtable(orig,'FileType','text','Delimiter','\t');
                csv_lbl=fullfile(lbl.folder,[m '.csv']);
                writetable(Tlbl,csv_lbl);
                lbl_path=csv_lbl;
            else
                lbl_path=orig;
            end
        else
            fprintf('Warning: No label file for %s\n',subj);
        end

        out_sub=fullfile(output_dir,subj); if ~exist(out_sub,'dir'), mkdir(out_sub); end

        try
            % Load localization CSV
            Loc=readtable(fpath);
            CT=table(Loc.x,Loc.y,Loc.z,'VariableNames',{'x','y','z'});
            CT.name=Loc.name;
            writetable(CT,fullfile(out_sub,'electrodes_coords.csv'));

            % Load label CSV
            if ~isempty(lbl_path)
                Label=readtable(lbl_path);
                if ~strcmp(Label.Properties.VariableNames{1},'label')
                    Label.Properties.VariableNames{1}='label';
                end
                writetable(Label,fullfile(out_sub,'electrodes_labels.csv'));
            else
                Label=table();
            end

            % Combine
            combined.coords=double([CT.x,CT.y,CT.z]);
            combined.names=CT.name;
            if istable(Label)&&height(Label)==height(CT)
                if ismember('mni_linear_x',Label.Properties.VariableNames)
                    combined.mni_coords=[Label.mni_linear_x,Label.mni_linear_y,Label.mni_linear_z];
                end
                if ismember('label',Label.Properties.VariableNames)
                    combined.labels=Label.label;
                end
                if ismember('HCPMMP1_label_1',Label.Properties.VariableNames)
                    combined.HCP_label=Label.HCPMMP1_label_1;
                end
                if ismember('HCPMMP1_weight_1',Label.Properties.VariableNames)
                    combined.HCP_weight=Label.HCPMMP1_weight_1;
                end
            end

            save(fullfile(out_sub,'electrodes_combined.mat'),'combined','-v7.3');
            processed{end+1}=subj;
            fprintf('Processed %s\n',subj);
        catch ME
            warning('Failed %s: %s',subj,ME.message);
            failed{end+1}=subj;
        end
    end

    % Summary
    fid=fopen(fullfile(output_dir,'subjects_summary.txt'),'w');
    fprintf(fid,'Processed subjects (%d):\n',numel(processed));
    for k=1:numel(processed)
        fprintf(fid,'  %d. %s\n',k,processed{k});
    end
    fclose(fid);
    fprintf('Done: %d success, %d failed\n',numel(processed),numel(failed));
end

function s=extract_subject_id(fn)
    tk=regexp(fn,'EM\d+','match'); s='';
    if ~isempty(tk), s=tk{1}; end
end

function f=find_label_file(list,subj)
    f=[]; for k=1:numel(list)
        if contains(list(k).name,subj), f=list(k); return; end
    end
end

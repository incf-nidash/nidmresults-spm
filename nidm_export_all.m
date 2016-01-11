function nidm_export_all(path_to, out_path)
    if ~isvarname('out_path')
        out_path = '';
    end

    spm('defaults','fmri');
    spm_jobman('initcfg');

    files = dir(path_to);
    subdirs = files([files.isdir]);
    for i = 1:numel(subdirs)
        dname = subdirs(i).name;
        if ~strcmp(dname, 'spm_explicit_mask')
            if strcmp(dname, 'spm_full_example001')
                % For SPM full example 001 we use already exported peaks 
                % and clusters list to get exactly the same graph
                load('nidm_example001.mat');
                SPM.swd=pwd;
                spm_results_nidm(SPM,xSPM,TabDat);
            else
                % TODO: read config.json instead to check for software used
                if strncmpi(dname,'spm_',4)
                    if ~strcmp(dname, 'ground_truth')
                        disp(dname)
                        nidm_export(fullfile(path_to, dname), out_path)
                    end
                end
            end
        end
    end
end
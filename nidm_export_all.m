function nidm_export_all(path_to)
    spm('defaults','fmri');
    spm_jobman('initcfg');

    files = dir(path_to);
    subdirs = files([files.isdir]);
    for i = 1:numel(subdirs)
        dname = subdirs(i).name;
        if ~strcmp(dname, 'spm_explicit_mask')
        % TODO: read config.json instead        
        if strncmpi(dname,'spm_',4)
            if ~strcmp(dname, 'ground_truth')
                disp(dname)
                nidm_export(fullfile(path_to, dname))
            end
        end
        end
    end
end
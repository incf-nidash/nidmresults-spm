function nidm_export(path_to_script_folder, out_path)
    cwd = pwd;
    cd(path_to_script_folder)
    % Remove previous nidm exports    
    files = dir(path_to_script_folder);
    subdirs = files([files.isdir]);
    for i = 1:numel(subdirs)
        dname = subdirs(i).name;
        if strncmpi(dname,'nidm',4)
            disp(['Removing ' dname])
            rmdir(dname,'s')
        end
        
        nidm_zips = cellstr(strvcat(spm_select('FPList', path_to_script_folder, '\.nidm\.zip$')));
        for j = 1:numel(nidm_zips)
            if ~isempty(nidm_zips{j})
                disp(['Deleting ' nidm_zips{j}])
                delete(nidm_zips{j})
            end
        end
    end
    
    run(fullfile(pwd, 'batch.m'))
    result_batch = matlabbatch(end);
    result_batch{1}.spm.stats.results.spmmat = {fullfile(pwd, 'SPM.mat')};
    result_batch{1}.spm.stats.results.print = 'nidm';    
    spm_jobman('run', result_batch)
    
    unzip('spm_0001.nidm.zip', 'nidm')
    
    if ~isempty(out_path)
        target_dir = fullfile(out_path, ['ex_' spm_file(path_to_script_folder, 'basename')]);
        if isdir(target_dir)
            disp(['Removing ' target_dir])
            rmdir(target_dir,'s')
        end
        movefile('nidm', target_dir)
    end
    
    cd(cwd);
end
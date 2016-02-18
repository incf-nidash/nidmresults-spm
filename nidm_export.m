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
    
    test_name = spm_file(path_to_script_folder, 'filename');
    if strcmp(test_name, 'spm_full_example001')
        % For SPM full example 001 we use already exported peaks 
        % and clusters list to get exactly the same graph
        load(fullfile(path_to_script_folder, 'nidm_example001.mat'));
        SPM.swd=pwd;
        spm_results_nidm(SPM,xSPM,TabDat);
    else    
        run(fullfile(pwd, 'batch.m'))
        result_batch = matlabbatch(end);
        result_batch{1}.spm.stats.results.spmmat = {fullfile(pwd, 'SPM.mat')};
        result_batch{1}.spm.stats.results.print = 'nidm';    
        spm_jobman('run', result_batch)
    end
    
    unzip('spm_0001.nidm.zip', 'nidm')
    
    if ~isempty(out_path)
        test_name = spm_file(path_to_script_folder, 'basename');
        
        target_dir = fullfile(out_path, ['ex_' test_name]);
        if isdir(target_dir)
            disp(['Removing ' target_dir])
            rmdir(target_dir,'s')
        end
        movefile('nidm', target_dir)
        json_file = fullfile(path_to_script_folder, 'config.json');
        copyfile(json_file, ...
                 fullfile(target_dir, 'config.json'));
             
        fname = json_file;
        fid = fopen(fname);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);

        expression = '\[".*"\]';
        gt = regexp(str,expression,'match'); 
        gt = strrep(strrep(strrep(gt{1}, '[', ''), ']', ''), '"', '');
        disp(gt)
        % FIXME: version should be extracted from json        
        gt_file = fullfile(path_to_script_folder, '..', 'ground_truth', '1.2.0', gt);
        
        target_gt_dir = fullfile(out_path, 'ground_truth', spm_file(gt,'path'));
        if isdir(target_gt_dir)
            disp(['Removing ' target_gt_dir])
            rmdir(target_gt_dir,'s')
        end
        mkdir(target_gt_dir)
        copyfile(gt_file, target_gt_dir);
    end
    
    cd(cwd);
end
%fileID = fopen('failed_files.txt','r');
%[failed_files_all] = fscanf(fileID, '%s');
%fclose(fileID);


% change '/' to '\' due to the difference between python and matlab
failed_files_all = strrep(fileread('failed_files.txt'),'/','\');
% replace folder
gap_sym = '\Volumes\behavgenom_archive$';

ini_loc = strfind(failed_files_all,gap_sym);
%ini_loc = regexp(failed_files_all,gap_sym);

file_name = {};

% restore file names to independent cell
for ii = 1:numel(ini_loc)-1
    file_name = [file_name;failed_files_all(ini_loc(ii):ini_loc(ii+1)-2)];
end
file_name = [file_name;failed_files_all(ini_loc(numel(ini_loc)):end)];

% record results
success = [];
% record number of failed files
num_fail = 0;
skip_file = [];
% loop goes through all files in txt
for iif = 131:numel(ini_loc);     
    
    % set current file and result hdf5 file
    cur_file = file_name{iif};
    result_file0 = strrep(cur_file, 'MaskedVideos', 'Results');
    result_file = strrep(result_file0, '.hdf5','_skeletons.hdf5');
    
    % show the progress
    fprintf('%i/%i) %s\n', iif, numel(ini_loc),cur_file)
    
    if ~isempty(regexp(cur_file, '\w*.hdf5', 'ONCE'))
        cur_suc = 0;
        % change group folder name to Z:
        masked_image_file = strrep(cur_file,gap_sym,'Z:');
        skeletons_file = strrep(result_file,gap_sym,'Z:');
        fprintf('%i) %s\n', iif, masked_image_file)
        video_timestamp_time = h5read(skeletons_file, '/timestamp/time');
        if video_timestamp_time(end)> 60*61;
            skip_file = [ skip_file; iif];
            % 
             fileID = fopen('long_files.txt','a');
             fprintf(fileID,'%s ',cur_file);
             fclose(fileID);
            continue;
        end
        try
            % record if it successes, 0 or 1
            
            cur_suc = alignStageMotionFun(masked_image_file,skeletons_file);   
        catch ME
            num_fail = num_fail+1;
        end
         success =[success;cur_suc];
         if cur_suc 
             % write to good file txt
             fileID = fopen('good_files.txt','a');
             fprintf(fileID,'%s ',cur_file);
             fclose(fileID);
         else
             % write to bad file txt
             fileID = fopen('bad_files.txt','a');
             fprintf(fileID,'%s ',cur_file);
             fclose(fileID);
         end
    else
        % write to error file txt
        fileID = fopen('error_files.txt','a');
        fprintf(fileID,'%s ',cur_file);
        fclose(fileID);
        error('file name error');
    end
   
    
end

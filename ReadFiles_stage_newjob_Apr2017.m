%fileID = fopen('failed_files.txt','r');
%[failed_files_all] = fscanf(fileID, '%s');
%fclose(fileID);

path =  'X:\Kezhi\FromAvelino\';
% change '/' to '\' due to the difference between python and matlab
failed_files_all = strrep(fileread([path,'missing_in_old_db.txt']),'/','\');
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
for iif = 1:numel(ini_loc);     
    
    % set current file and result hdf5 file
    cur_file_now = strtrim(file_name{iif});
    result_file = strrep(cur_file_now, '.hdf5','_skeletons.hdf5');
    
%     % set current file and result hdf5 file
%     result_file = strtrim(file_name{iif});
%     cur_file = strrep(result_file,'_skeletons.hdf5', '.hdf5');
%     cur_file_now = strrep(cur_file, 'Results', 'MaskedVideos');
%
%     % use MaskedVideos_old here
%     %cur_file_now = strrep(cur_file, 'MaskedVideos', 'MaskedVideos_old');
%     cur_file_now = strrep(cur_file, 'MaskedVideos', 'MaskedVideos_old');
%     
%     result_file0 = strrep(cur_file_now, 'MaskedVideos', 'Results');
%     %result_file0 = strrep(cur_file_now, 'MaskedVideos_old', 'Results_old');
%     
%     result_file = result_file0;
%     %result_file = strrep(result_file0, '.hdf5','_skeletons.hdf5');
    
%     cur_file = strtrim(file_name{iif});
%     cur_file_now = strrep(cur_file, 'Results', 'MaskedVideos');
%     result_file = strrep(cur_file, '.hdf5','_skeletons.hdf5');
    
    % show the progress
    fprintf('%i/%i) %s\n', iif, numel(ini_loc),cur_file_now)
    
    if ~isempty(regexp(result_file, 'swimming', 'ONCE'))
        fileID = fopen('missing_swimming_files.txt','a');
        fprintf(fileID,'%s ',result_file);
        fclose(fileID);
        continue;
    elseif ~isempty(regexp(cur_file_now, '\w*.hdf5', 'ONCE'))
        cur_suc = 0;
        % change group folder name to Z:
        masked_image_file = strrep(cur_file_now,gap_sym,'Z:');
        skeletons_file = strrep(result_file,gap_sym,'Z:');
        fprintf('%i) %s\n', iif, masked_image_file)
        try
            video_timestamp_time = h5read(skeletons_file, '/timestamp/time');
        catch ME
            num_fail = num_fail+1;
        end
        if video_timestamp_time(end)> 60*61;
            skip_file = [ skip_file; iif];
            % 
             fileID = fopen('missing_long_files.txt','a');
             fprintf(fileID,'%s ',cur_file_now);
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
             fileID = fopen([path,'missing_good_files.txt'],'a');
             fprintf(fileID,'%s ',cur_file_now);
             fclose(fileID);
         else
             % write to bad file txt
             fileID = fopen([path,'missing_bad_files.txt'],'a');
             fprintf(fileID,'%s ',cur_file_now);
             fclose(fileID);
         end
    else
        % write to error file txt
        fileID = fopen([path,'missing_error_files.txt'],'a');
        fprintf(fileID,'%s ',cur_file_now);
        fclose(fileID);
        error('file name error');
    end
   
    
end

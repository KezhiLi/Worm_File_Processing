function [match_score, match_mtx] = cal_match_score(cvs_ind, mov_fra_ind, mov_fra_ind_beg, para_extend, csv_mag, frame_mag)
% 
% 
% 
% 
% 
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 03/02/2016
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.

% 
% debug 
if ((mov_fra_ind_beg == 1)||(mov_fra_ind_beg == 2)) && para_extend == 1
    para_extend
end


% initialize
leng_cvs_ind = size(cvs_ind,1);
leng_csv = size(csv_mag,1);
weights = exp(0:-(2/(leng_cvs_ind-1)):-2);

% match_mtx is used to record all possible penalties for different mm1 and
% mm2
match_mtx = zeros(leng_cvs_ind,4);
match_mtx(:,1) = cvs_ind; 

mov_fra_ind_extend = round((mov_fra_ind(mov_fra_ind_beg:end)-...
    mov_fra_ind(mov_fra_ind_beg))*para_extend + cvs_ind(1));

% parameters to calculate the 'shift_to_left'
shift_weights = exp(-2.4:0.6:0)';
normal_shift_weights = shift_weights/sum(shift_weights);
% only conder nearest 10 peaks to calculate the 'shift_to_left'
shift_consider = zeros(1,length(shift_weights));

% compensate the difference: shift impulse index to left 
shift_to_left = zeros(leng_cvs_ind+1,1);



mov_fra_ind_extend_copy = mov_fra_ind_extend;

for ii = 1:leng_cvs_ind;
    [~,min_ind] = min(abs(mov_fra_ind_extend_copy-cvs_ind(ii)-shift_to_left(ii)));
    match_ind = mov_fra_ind_extend_copy(min_ind);
    match_mtx(ii,2) = match_ind;
    
    % update 'shift_to_left'. it is caluclated based on last shift values
    shift_to_left(ii) = match_mtx(ii,2) -match_mtx(ii,1);
    shift_consider(1:end-1) = shift_consider(2:end);
    %shift_consider(end) = min_fra_ind_match(ii,5)-mask_ind_x(ii,1);
    shift_consider(end) =  shift_to_left(ii);
    shift_to_left(ii+1) = round(shift_consider*normal_shift_weights);
    
    mov_fra_ind_extend_copy(min_ind) = 0;
end

for ii = 1:leng_cvs_ind;
    match_mtx(ii,3) = csv_mag(cvs_ind(ii));
    match_mtx(ii,4) = sum(frame_mag(max(1,match_mtx(ii,2)-3):min(match_mtx(ii,2)+3,leng_csv)));
end

zero_ind = find(match_mtx(:,2)==0);
positive_ind = find(match_mtx(:,2)>0);
%if sum(abs(sort(match_mtx(:,2))-match_mtx(:,2)))>1
if min(zero_ind) < max(positive_ind)
    match_score = 1e10;
else
    % absulute difference/ and/ second derivitive 
    match_diff = match_mtx(:,2)-match_mtx(:,1);
    
    % distance from mm1=1, mm2=20
    mm1_add = mov_fra_ind_beg-1;
    mm2_add = abs(round((para_extend-0.94)/0.003)-20);
    
    % % calculate the match score, given weights, and distance mm1=1, mm2=20
    % match_score = sum(abs(match_diff.*weights'))+(mm1_add/2+ mm2_add); % +sum(abs(match_diff(2:end)-match_diff(1:end-1)));
    match_score = sum(abs(match_diff.*weights'))*(1+(mm1_add+ mm2_add*2)*0.02);
end



function [match_score, match_mtx] = cal_match_score(cvs_ind, mov_fra_ind, mov_fra_ind_beg, para_extend, csv_mag, frame_mag)
% 
% 
% 
% 
% 
% 
% 
% 

% 
% debug 
% if ((mov_fra_ind_beg == 1)||(mov_fra_ind_beg == 2)) && para_extend == 1.09
%     para_extend
% end

% initialize
leng_cvs_ind = size(cvs_ind,1);
leng_csv = size(csv_mag,1);
weights = exp(0:-(2/(leng_cvs_ind-1)):-2);


match_mtx = zeros(leng_cvs_ind,4);
match_mtx(:,1) = cvs_ind; 

mov_fra_ind_extend = round((mov_fra_ind(mov_fra_ind_beg:end)-...
    mov_fra_ind(mov_fra_ind_beg))*para_extend + cvs_ind(1));

% cvs_compare = zeros(leng_csv,1);
% mov_compare = zeros(leng_csv,1);
% cvs_compare(cvs_ind) = 1;
% mov_compare(mov_fra_ind_extend) = 50;
% figure, plot(csv_mag,'r-o');
% hold on , plot(mov_compare,'b-');
% plot (frame_mag,'g-');


mov_fra_ind_extend_copy = mov_fra_ind_extend;

for ii = 1:leng_cvs_ind;
    [~,min_ind] = min(abs(mov_fra_ind_extend_copy-cvs_ind(ii)));
    match_ind = mov_fra_ind_extend_copy(min_ind);
    match_mtx(ii,2) = match_ind;
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
    match_score = sum(abs(match_diff.*weights')); % +sum(abs(match_diff(2:end)-match_diff(1:end-1)));
end



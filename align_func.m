function stage_move_new = align_func(stage_move, central_move, stage_impulse,range1, mask_ind, csv_ind, ii)
% function used to align the impulse stage moving in csv file and the real
% stage moving distance in each time slot according to time stamps
% Input: 
% stage_move: the vector to update the results; at first it is a
% all-0 vector; in each iteration, the function update a bunch of moving
% distance corresponding to 1 impluse stage move in csv
% central_move: the vector that records the moving distance of the central point of the worm
% stage_impulse: a vector that saves the distance of impulse stage moving
% range1: the considering range on left/right side around the corresponding predicted
% time stamp
% mask_ind: frame indexes 
% csv_ind: impluse indexes in the csv file
% ii: the index of iteration, corresponding to the index of stage moving 
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 10/12/2015
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.

stage_move_new = stage_move;

% debug purpose
if ii == 199
ii
end
% initialize
    x_align1 = zeros(2*range1+1,1);
    x_align2 = zeros(2*range1+1,1);
    x_align3 = zeros(2*range1+1,1);
    x_align4 = zeros(2*range1+1,1);

    pulse_cental_ii = mask_ind(ii,1);
    % consider -range1 to range1, in total 2*range1+1 points
    x_consider = central_move(max(1,pulse_cental_ii-range1):min(length(central_move),pulse_cental_ii+range1));
    % true stage moving pixels
    x_stagemove = stage_impulse(csv_ind(ii));
    
    if abs(x_stagemove)>1e-10
        


        % set a threshold to distinguish impulse part and other part
        x_consider_normalized = (x_consider-min(x_consider))/max(x_consider-min(x_consider));
        x_level = graythresh(x_consider_normalized)*max(x_consider-min(x_consider))+min(x_consider);
        % find the central impulse indexes
    %     if x_consider(range1+1,1)>0
    %         impulse_ind_binary = x_consider>x_level;       
    %     else
    %         impulse_ind_binary = x_consider<x_level;    
    %     end
        impulse_ind_binary = abs(x_consider)>abs(x_level);
        % only keep the central impulses, cancel other fake ones, by
        % identifying the 2 closest change points near central point
        impulse_ind_binary_diff = abs(impulse_ind_binary(2:end) - impulse_ind_binary(1:end-1));
        if sum(impulse_ind_binary_diff) == 2
            chg_pt1 = find(impulse_ind_binary_diff,1,'first')+1;
            chg_pt2 = find(impulse_ind_binary_diff,1,'last');
        else
            if impulse_ind_binary(range1+1)<0.5 % or ==0 % means the peak is not at the central location
                % if there exists peak in the x_consider
                if sum(impulse_ind_binary)>=2
                    chg_ind = find(impulse_ind_binary_diff>0.5); % >0
                    [~, sort_chg_ind_ind] = sort(abs(chg_ind-range1));
                    chg_pt1 = chg_ind(sort_chg_ind_ind(1));            
                    if length(sort_chg_ind_ind)>1
                        chg_pt2 = chg_ind(sort_chg_ind_ind(2));
                    else % postive values last to one of the end of the vector
                        if abs(chg_pt1-1)<abs(chg_pt1-2*range1)
                            chg_pt2 = 1;
                        else
                            chg_pt2 = 2*range1;
                        end
                    end
                    if chg_pt1>chg_pt2
                        temp_chg_pt1=chg_pt2;
                        chg_pt2 = chg_pt1;
                        chg_pt1 = temp_chg_pt1;
                    end
                else
                    error_peak = 'peak error'
                    chg_pt1  = range1-5;
                    chg_pt2  = range1+5;
                end
            else
                for jj = 1: range1-3;
                    if impulse_ind_binary(range1+1-jj)<0.5  % or ==0
                        chg_pt1 = range1+2-jj;
                        break;
                    end
                    if jj == range1-3
                        chg_pt1 = 4;
                    end
                end
                for jj = 1: range1-3;
                    if impulse_ind_binary(range1+1+jj)<0.5  % or ==0
                        chg_pt2 = range1+jj;
                        break;
                    end
                    if jj == range1-3
                        chg_pt2 = 2*range1-2;
                    end
                end
            end
        end

        % compare two window size
        % align initial
        x_align1(chg_pt1:chg_pt2) = x_consider(chg_pt1:chg_pt2); 
        % compensate average
        x1_ave = (stage_impulse(csv_ind(ii))-sum(x_align1(chg_pt1:chg_pt2)))/(chg_pt2-chg_pt1+1);
        % align results
        x_align1(chg_pt1:chg_pt2)=x_align1(chg_pt1:chg_pt2)+x1_ave;
        % left worm moving by deleting stage move
        x1_left = x_consider;
        x1_left(chg_pt1:chg_pt2) = -x1_ave;
        var_x1_left = var(x1_left);

        chg_pt1
        chg_pt2
        % align initial, left 1 more, right 1 more
        x_align2(max(1,chg_pt1-1):chg_pt2+1) = x_consider(chg_pt1-1:chg_pt2+1); 
        % compensate average
        x2_ave = (x_stagemove-sum(x_align2(max(1,chg_pt1-1):chg_pt2+1)))/(chg_pt2-chg_pt1+3);
        % align results
        x_align2(max(1,chg_pt1-1):chg_pt2+1)=x_align2(max(1,chg_pt1-1):chg_pt2+1)+x2_ave;
        % left worm moving by deleting stage move
        x2_left = x_consider;
        x2_left(max(1,chg_pt1-1):chg_pt2+1) = -x2_ave;
        var_x2_left = var(x2_left);

        % align initial, only left 1 more  
        x_align3(max(1,chg_pt1-1):chg_pt2) = x_consider(max(1,chg_pt1-1):chg_pt2); 
        % compensate average
        x3_ave = (x_stagemove-sum(x_align3(max(1,chg_pt1-1):chg_pt2)))/(chg_pt2-chg_pt1+2);
        % align results
        x_align3(max(1,chg_pt1-1):chg_pt2)=x_align3(max(1,chg_pt1-1):chg_pt2)+x3_ave;
        % left worm moving by deleting stage move
        x3_left = x_consider;
        x3_left(max(1,chg_pt1-1):chg_pt2) = -x3_ave;
        var_x3_left = var(x3_left);

        % align initial, only right 1 more
        x_align4(chg_pt1:chg_pt2+1) = x_consider(chg_pt1:chg_pt2+1); 
        % compensate average
        x4_ave = (x_stagemove-sum(x_align4(chg_pt1:chg_pt2+1)))/(chg_pt2-chg_pt1+2);
        % align results
        x_align4(chg_pt1:chg_pt2+1)=x_align4(chg_pt1:chg_pt2+1)+x4_ave;
        % left worm moving by deleting stage move
        x4_left = x_consider;
        x4_left(chg_pt1:chg_pt2+1) = -x4_ave;
        var_x4_left = var(x4_left);


        % compare variance, and make sure that for each moving, direction
        % (postive/negtive) is the same
        [var_sort, var_ind] = sort([var_x1_left,var_x2_left,var_x3_left,var_x4_left]);
        align_mtx = [x_align1,x_align2,x_align3,x_align4];

        for nn = 1:4;
            nn_var = var_ind(nn);
            if abs(sum(sign(align_mtx(:,nn_var))))==sum(abs(sign(align_mtx(:,nn_var))))
                stage_move_new(pulse_cental_ii-range1:pulse_cental_ii+range1) = align_mtx(:,nn_var);
                break;
            elseif nn ==4
                stage_move_new(pulse_cental_ii-range1:pulse_cental_ii+range1) = align_mtx(:,var_ind(1));
            end
        end

        % check length
        if length(stage_move_new)>length(stage_move)
            stage_move_new = stage_move_new(1:length(stage_move));
        elseif length(stage_move_new)<length(stage_move)
            stage_move_new = [stage_move_new;zeros(length(stage_move)-length(stage_move_new),1)];
        end
    end
end
%     if (var_x1_left > var_x2_left) && (sum(sign(x_align2))==sum(abs(sign(x_align2)))) 
%         stage_move(pulse_cental_ii-range1:pulse_cental_ii+range1) = x_align2;
%     else
%         stage_move(pulse_cental_ii-range1:pulse_cental_ii+range1) = x_align1;
%     end
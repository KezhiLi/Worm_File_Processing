load(['N:\Andre\results-12-05-10\Laura Grundy\egl-17\e1313\CB1313\on_food\XX\30m_wait\R\tracker_3\2010-07-16___13_06_00\',name_temp,'_features.mat'])
% N:\Andre\results-12-05-10\Laura Grundy\egl-17\e1313\CB1313\on_food\XX\30m_wait\L\tracker_2\2010-07-09___11_43_13\

Nannan_ind = worm.posture.skeleton.x(1,:)>0;
worm_skeleton_x = worm.posture.skeleton.x(:,Nannan_ind);
worm_skeleton_y = worm.posture.skeleton.y(:,Nannan_ind);

tt_step = 10;

%%%%%%

x_ske_cum_left = x_ske_cum(left_fra_ind,:)*(-x_pixel_per_microns);  % x_ske_cum
y_ske_cum_left = y_ske_cum(left_fra_ind,:)*(-y_pixel_per_microns);  % y_ske_cum

Nannan_ind2 = ~isnan(x_ske_cum_left(:,1));
x_ske_cum_left2 = x_ske_cum_left(Nannan_ind2,:)+(worm_skeleton_x(1,1)-x_ske_cum_left(1,1));
y_ske_cum_left2 = y_ske_cum_left(Nannan_ind2,:)+(worm_skeleton_y(1,1)-y_ske_cum_left(1,1));

fig66 = figure(66),plot(x_ske_cum_left2(1,:),y_ske_cum_left2(1,:) ); axis equal; hold on

for tt_1 = 2:tt_step:size(x_ske_cum_left2,1);
    %f ~ismember(tt_1+1, cancel_fra_ind)
%         plot(x_ske_cum(tt_1,:),y_ske_cum(tt_1,:),'g','LineWidth',1.1);
%     else
        plot(x_ske_cum_left2(tt_1,:),y_ske_cum_left2(tt_1,:) );
    %end
end
hold off,


% show frames shown
size(x_ske_cum_left2,1)

%%%%%%%%%%%%%%

fig67 = figure(67),plot(worm_skeleton_x(:,1),worm_skeleton_y(:,1) ); axis equal; hold on

for tt_1 = 2:tt_step:size(worm_skeleton_x,2);
    %if ismember(tt_1+1, cancel_fra_ind)
    if ~(isnan(sum(worm_skeleton_x(:,tt_1))) | isnan(sum(worm_skeleton_y(:,tt_1))))
%         plot(worm.posture.skeleton.x(:,tt_1),worm.posture.skeleton.y(:,tt_1),'g','LineWidth',1.1);
%     else
        plot(worm_skeleton_x(:,tt_1),worm_skeleton_y(:,tt_1) );
    end
end
hold off,

% show frames shown
size(worm_skeleton_x,2)
%%%%%

savefig(fig66,[trajectories_folder,name,'-align_ske.fig']);
savefig(fig67,[trajectories_folder,name,'-real_ske.fig']);

% 
% frame_all = zeros(1, size( worm.posture.skeleton.x,2));
% frame_all1 = frame_all;
% frame_all1(diff_mask_central(cancel_fra_ind,8)) = 1;
% frame_stage = frame_all;
% frame_stage(info.video.annotations.frames==3) = 1.3;
% 
% fig71 = figure(71), plot(frame_all1,'g'); 
% hold on, plot(frame_stage,'r')
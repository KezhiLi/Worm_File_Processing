
clear
clc

temp_text_file = ['temp_text',num2str(nf),'.mat'];

load(temp_text_file );


load([path2,name,'_features.mat']);


% 
% % vr = VideoReader([hdf5_path(1:end-5),'.avi']);
% 
% if ~exist('mask')
%     mask = h5read(hdf5_path, '/mask');
% end
% 
% 
% stage_move = [stage_move_x, stage_move_y];
% 
% for ii = 1:size(mask,3)-1;
%     ii
%     % Y_k = read(vr, ii); % first frame
%     Y_k = mask(:,:,ii);
%     figure(1)
%     if sum(abs(stage_move(ii,:)))<0.01;
%         image((Y_k));
%         title('+++ Find Stage Move +++');
%         pause(0.02);
%     else
%         mask_temp = (Y_k);
%         mask_temp([1:25,end-24:end],:) =  uint8(round(rand(size(mask_temp([1:25,end-24:end],:)))*254));
%         mask_temp(:,[1:25,end-24:end]) =  uint8(round(rand(size(mask_temp(:,[1:25,end-24:end])))*254));
%         image(mask_temp);
%         title('+++ Find Stage Move +++');
%         pause(0.2);
%     end
% end

x_shift = diff_mask_central(:,1);
y_shift = diff_mask_central(:,2);
x_shift(1) = x_shift(1) + worm.posture.skeleton.x(25, 1);
y_shift(1) = y_shift(1) + worm.posture.skeleton.y(25, 1);

x_shift((abs(diff_mask_central(:,9))>0))=0;
x_shift((abs(diff_mask_central(:,10))>0))=0;
y_shift((abs(diff_mask_central(:,9))>0))=0;
y_shift((abs(diff_mask_central(:,10))>0))=0;

x_summ = cumsum(x_shift);%*abs(x_pixel_per_microns);
y_summ = cumsum(y_shift);%*abs(y_pixel_per_microns);

figure, plot(x_summ,y_summ,'r'); axis equal
figure, plot(worm.posture.skeleton.x(25, :), worm.posture.skeleton.y(25, :)); axis equal
load('C:\Kezhi\MyCode!!!\ManualVideos\Check_Align_samples\npr-10 (tm1568) on food R_2010_01_25__14_23_46___4___6_features.mat')

% figure(8),plot(worm.posture.skeleton.y(25, 1), worm.posture.skeleton.x(25, 1)); %axis equal
% figure(9), plot(x_summ(1),y_summ(1));
% for ii = 1:length(x_summ);
%     ii
%     figure(7),image(mask(:,:,ii))
%     figure(8),hold on, 
%     plot( worm.posture.skeleton.x(25, ii), worm.posture.skeleton.y(25, ii),'b.');
%     figure(9),hold on, 
%     plot(x_summ(ii),y_summ(ii),'r.');
%     pause(0.00001);
% end

% x_centr = -diff_mask_central(:,9);
% y_centr = -diff_mask_central(:,10);
% x_centr(1) = x_centr(1) + worm.posture.skeleton.x(25, 1);
% y_centr(1) = y_centr(1) + worm.posture.skeleton.y(25, 1);
% 
% % x_centr_summ = round(cumsum(x_centr)*abs(x_pixel_per_microns));
% % y_centr_summ = round(cumsum(y_centr)*abs(y_pixel_per_microns));
% x_centr_summ = round(cumsum(x_centr));
% y_centr_summ = round(cumsum(y_centr));
% 
% x_centr_summ_adjusted = x_centr_summ - 7.95e3;
% y_centr_summ_adjusted = y_centr_summ - 4.17e4;
% 
% curr_img = uint8(zeros(2600,2500));
% half_x = 320;
% half_y = 240;
% 
% for ii = 1:frame_total-1;
%     ii
%     mask_backadjust = mask(:,:,ii+1);
%     mask_backadjust(mask_backadjust<1)=substitute_intensity;
%     
%     mask_backadjust(1:4,:) = 256;
%     mask_backadjust(end-3:end,:) = 256;
%     mask_backadjust(:,1:4) = 256;
%     mask_backadjust(:,end-3:end) = 256;
%     
%     curr_img((y_centr_summ_adjusted(ii)-half_y):(y_centr_summ_adjusted(ii)+half_y-1),...
%        (x_centr_summ_adjusted(ii)-half_x):(x_centr_summ_adjusted(ii)+half_x-1) ) = (mask_backadjust)';
%     figure(11), imshow(curr_img); 
% end




    

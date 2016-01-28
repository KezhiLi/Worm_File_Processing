function ss = read_csv_data(root,excel_name)
%
%
%
%
%
%
%
% @Copyright Kezhi Li, Imperial College London, Medical Research Council,
% Jan 25, 2016

ss = importdata([root,excel_name]);

% fix the problem when the title is not in the first row of csv file
if size(ss.textdata,1) < 5
    ss_xlsread = xlsread([root,excel_name]);
    NaN_row = find(isnan(ss_xlsread(:,1)),1); 
    if size(NaN_row,1) >0
        %1st_row = ss_xlsread(1,:);
        ss_xlsread(2:NaN_row,:) = ss_xlsread(1:NaN_row-1,:);    
        ss_xlsread(1,:) = NaN;     
    else
        error('csv file has a error');
    end
    ss.textdata = ss_xlsread(:,1:2)*24*3600;
    ss.data = ss_xlsread(2:end,4:5);
end

% fix the problem when there is duplicated items in ss
[~,index] = unique(ss.textdata(:,2),'first');   
if length(index)<size(ss.textdata,1)
    %ss.textdata = ss.textdata(1:end-(size(ss.textdata,1)-length(index)),:);
    ss.textdata(2:length(index),:) = ss.textdata(index(1:end-1),:);
    ss.textdata = ss.textdata(1:end-(size(ss.textdata,1)-length(index)),:);
    ss.data(1:length(index)-1,:) = ss.data(index(1:end-1)-1,:);
    ss.data = ss.data(1:end-(size(ss.data,1)-length(index)+1),:);
    
end





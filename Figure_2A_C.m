%% time series plots (M-TsData)
% plot example voxel ts & prediction
% find example voxel id-perceive
% amp, corr and map
% co: variance explained
% amp: sigma (i.e. tuning width)
% map: eccentricity (i.e. numerosity preference)
% Subject directory

clear;
mrVista 3

clear control_file uncon_file con_file
% unconnect files:
uncon_file =dir(fullfile(sub_dir,uncon_dir,'*gFit-gFit.mat')).name;

% connect files
con_file = dir(fullfile(sub_dir,con_dir,'*gFit-gFit.mat')).name;

% control files
control_file = dir(fullfile(sub_dir,control_dir,'*gFit-gFit.mat')).name;

roi_dir = '\Gray\ROIs\';
filename = dir([sub_dir,roi_dir,'LN*.mat']);
roiName = {};
for filenum = 1:length(filename)
    roiName{filenum} = filename(filenum).name(1:length(filename(filenum).name)-4);
end

condition = {'control','uncon', 'con'};
condfile = {control_file, uncon_file, con_file};
conddir =  {control_dir, uncon_dir, con_dir};

% change dataType to control
condnum = 1;
VOLUME{1} = rmSelect(VOLUME{1}, 1, fullfile(sub_dir,conddir{1,condnum},condfile{1,condnum}));
VOLUME{1} = rmLoadDefault(VOLUME{1},0);
for i = 1:length(roiName)
    [VOLUME{1}, ok] = loadROI(VOLUME{1}, [sub_dir,roi_dir,roiName{i},'.mat'],1,'b',1);
    M = rmPlotGUI(VOLUME{1});
    TsData.(roiName{i}).(condition{condnum}) = M;
    for k = 1:length(TsData.(roiName{i}).(condition{condnum}).coords)
        [TsData.(roiName{i}).(condition{condnum}).prediction(:,k),~,TsData.(roiName{i}).(condition{condnum}).rfParams(k,:),~,blanks] = ...
            rmPlotGUI_makePrediction(TsData.(roiName{i}).(condition{condnum}),[],k);
    end
    thresh = ve.(roiName{i}).(condition{condnum}) > 0.3 & pref.(roiName{i}).(condition{condnum}) >=1.05 & pref.(roiName{i}).(condition{condnum}) <= 7;
    id = find(ve.(roiName{i}).(condition{condnum}) > max(ve.(roiName{i}).(condition{condnum})(thresh))-0.001);
    tmp = reshape(TsData.(roiName{i}).(condition{condnum}).tSeries(:,id(1)),[60 4]);
    for j = 1:60  
        minfo = wstat(tmp(j,:));
        ts.(roiName{i}).(condition{condnum})(j) = minfo.mean;
        sterr.(roiName{i}).(condition{condnum})(j) = minfo.sterr;
    end
    tmp = reshape(TsData.(roiName{i}).(condition{condnum}).prediction(:,id(1)),[60 4]);
    for j = 1:60
        minfo = wstat(tmp(j,:));
        pred.(roiName{i}).(condition{condnum})(j) = minfo.mean;
    end
    % figure;
    % errorbar(TsData.(roiName{i}).(condition{condnum}).x(1:60,1),ts.(roiName{i}).(condition{condnum}),sterr.(roiName{i}).uncon,'ko');
    % hold on; plot(TsData.(roiName{i}).(condition{condnum}).x(1:60,1),pred.(roiName{i}).(condition{condnum}),'b-');
    % axis square; box off; 
    % xlabel('time (sec)'); ylabel('bold signal change (%)')
    % %axis([0 90 -2.5 2.5]);
    % title(sprintf('Sub07 %s ve-%s = %.2f',roiName{i},condition{condnum},ve.(roiName{i}).(condition{condnum})(id(1))))

end

% change dataTYPE to unconnect
condnum = 2;
VOLUME{1} = rmSelect(VOLUME{1}, 1, fullfile(sub_dir,conddir{1,condnum},condfile{1,condnum}));
VOLUME{1} = rmLoadDefault(VOLUME{1},0);
for i = 1:length(roiName)
    [VOLUME{1}, ok] = loadROI(VOLUME{1}, [sub_dir,roi_dir,roiName{i},'.mat'],1,'b',1);
    M = rmPlotGUI(VOLUME{1});
    TsData.(roiName{i}).(condition{condnum}) = M;
    for k = 1:length(TsData.(roiName{i}).(condition{condnum}).coords)
        [TsData.(roiName{i}).(condition{condnum}).prediction(:,k),~,TsData.(roiName{i}).(condition{condnum}).rfParams(k,:),~,blanks] = ...
            rmPlotGUI_makePrediction(TsData.(roiName{i}).(condition{condnum}),[],k);
    end
    thresh = ve.(roiName{i}).(condition{condnum}) > 0.3 & pref.(roiName{i}).(condition{condnum}) >=1.05 & pref.(roiName{i}).(condition{condnum}) <= 7;
    id = find(ve.(roiName{i}).(condition{condnum}) > max(ve.(roiName{i}).(condition{condnum})(thresh))-0.001);
    tmp = reshape(TsData.(roiName{i}).(condition{condnum}).tSeries(:,id(1)),[60 4]);
    for j = 1:60  
        minfo = wstat(tmp(j,:));
        ts.(roiName{i}).(condition{condnum})(j) = minfo.mean;
        sterr.(roiName{i}).(condition{condnum})(j) = minfo.sterr;
    end
    tmp = reshape(TsData.(roiName{i}).(condition{condnum}).prediction(:,id(1)),[60 4]);
    for j = 1:60
        minfo = wstat(tmp(j,:));
        pred.(roiName{i}).(condition{condnum})(j) = minfo.mean;
    end
    % figure;
    % errorbar(TsData.(roiName{i}).(condition{condnum}).x(1:60,1),ts.(roiName{i}).(condition{condnum}),sterr.(roiName{i}).uncon,'ko');
    % hold on; plot(TsData.(roiName{i}).(condition{condnum}).x(1:60,1),pred.(roiName{i}).(condition{condnum}),'b-');
    % axis square; box off; 
    % xlabel('time (sec)'); ylabel('bold signal change (%)')
    % %axis([0 90 -2.5 2.5]);
    % title(sprintf('Sub07 %s ve-%s = %.2f',roiName{i},condition{condnum},ve.(roiName{i}).(condition{condnum})(id(1))))

end
%exportgraphics(gcf,outputfile,'ContentType','vector','Append',true)

% change dataTYPE to connect
condnum = 3;
VOLUME{1} = rmSelect(VOLUME{1}, 1, fullfile(sub_dir,conddir{1,condnum},condfile{1,condnum}));
VOLUME{1} = rmLoadDefault(VOLUME{1},0);
for i = 1:length(roiName)
    [VOLUME{1}, ok] = loadROI(VOLUME{1}, [sub_dir,roi_dir,roiName{i},'.mat'],1,'b',1);
    M = rmPlotGUI(VOLUME{1});
    TsData.(roiName{i}).(condition{condnum}) = M;
    for k = 1:length(TsData.(roiName{i}).(condition{condnum}).coords)
        [TsData.(roiName{i}).(condition{condnum}).prediction(:,k),~,TsData.(roiName{i}).(condition{condnum}).rfParams(k,:),~,blanks] = ...
            rmPlotGUI_makePrediction(TsData.(roiName{i}).(condition{condnum}),[],k);
    end
    thresh = ve.(roiName{i}).(condition{condnum}) > 0.3 & pref.(roiName{i}).(condition{condnum}) >=1.05 & pref.(roiName{i}).(condition{condnum}) <= 7;
    id = find(ve.(roiName{i}).(condition{condnum}) > max(ve.(roiName{i}).(condition{condnum})(thresh))-0.001);
    tmp = reshape(TsData.(roiName{i}).(condition{condnum}).tSeries(:,id(1)),[60 4]);
    for j = 1:60
        minfo = wstat(tmp(j,:));
        ts.(roiName{i}).(condition{condnum})(j) = minfo.mean;
        sterr.(roiName{i}).(condition{condnum})(j) = minfo.sterr;
    end
    tmp = reshape(TsData.(roiName{i}).(condition{condnum}).prediction(:,id(1)),[60 4]);
    for j = 1:60
        minfo = wstat(tmp(j,:));
        pred.(roiName{i}).(condition{condnum})(j) = minfo.mean;
    end
    % figure;
    % errorbar(TsData.(roiName{i}).(condition{condnum}).x(1:60,1),ts.(roiName{i}).(condition{condnum}),sterr.(roiName{i}).uncon,'ko');
    % hold on; plot(TsData.(roiName{i}).(condition{condnum}).x(1:60,1),pred.(roiName{i}).(condition{condnum}),'b-');
    % axis square; box off; 
    % xlabel('time (sec)'); ylabel('bold signal change (%)')
    % %axis([0 90 -2.5 2.5]);
    % title(sprintf('Sub07 %s ve-%s = %.2f',roiName{i},condition{condnum},ve.(roiName{i}).(condition{condnum})(id(1))))  
end
%exportgraphics(gcf,outputfile,'ContentType','vector','Append',true)

% save([sub_dir,'\results\tsparam.mat'],"ts","sterr","pred","TsData")
color = [1 0.5 1; 0.5 0.5 1; 0.5 0 1; 1 0.8 0.5; 0 0.5 0;0.5 1 1; 1 0.5 0.5];
conditions = {'control', 'uncon', 'con'};   % conditions: con - connected; uncon - unconnected
conditionColors = containers.Map({'control', 'uncon', 'con'}, {[0 0 1], [1 0.5 0], [0 1 0]}); 

for i = 1:length(roiName)
    for j = 1:length(conditions)
        condition = conditions{j};
        figure;

        thresh = ve.(roiName{i}).control > 0.3 & pref.(roiName{i}).control >=1.05 & pref.(roiName{i}).control <= 7;
        id = find(ve.(roiName{i}).control > max(ve.(roiName{i}).control(thresh))-0.001);


        tmp_ts = reshape(TsData.(roiName{i}).(condition).tSeries(:,id(1)),[60 4]);
        ts_data = zeros(1, 60);
        sterr_data = zeros(1, 60);
        for k = 1:60
            minfo = wstat(tmp_ts(k,:));
            ts_data(k) = minfo.mean;
            sterr_data(k) = minfo.sterr;
        end

        tmp_pred = reshape(TsData.(roiName{i}).(condition).prediction(:,id(1)),[60 4]);
        pred_data = zeros(1, 60);
        for k = 1:60
            minfo = wstat(tmp_pred(k,:));
            pred_data(k) = minfo.mean;
        end

        current_color = conditionColors(condition);
        errorbar(TsData.(roiName{i}).(condition).x(1:60,1),ts_data,sterr_data,'o','color',current_color);
        hold on;
        plot(TsData.(roiName{i}).(condition).x(1:60,1),pred_data,'-','color',current_color, 'LineWidth', 1.5); 
        % plot(TsData.(roiName{i}).(condition).x(1:60,1),pred_data,'-','color',color(1,1:3), 'LineWidth', 1.5);

        axis square;
        box off;
        xlabel('time (sec)'); ylabel('bold signal change (%)')
        %axis([0 90 -0.5 0.5]);
        title(sprintf('%s - %s, ve = %.2f',roiName{i}, condition, ve.(roiName{i}).(condition)(id(1)))) 
        legend([condition ' raw'],[condition ' fit']);
        hold off; 
    end
end
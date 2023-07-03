clearvars;close all;clc % cleaning up

load(fullfile(pwd,'data','CGC_FigParameters')) % fig parameters load
load(fullfile(pwd,'data','imgStats_exp2'))

PlotFig14a=1;
PlotFig14b=1;

% Define max value for Pellacini's c
maxVal = 0.26;minVal = -0.01;

% A list of observers
observerList = {'AG','AM','AZ','ET','LH','NJ','NZ','OP','YW','SB'};

sessionN = 2; % Number of sessions

for observer = 1:length(observerList)
    for session = 1:sessionN
        load(fullfile(repo_basedir,'data','rawdata','exp2',[observerList{observer},'_session',num2str(session),'.mat']))
        [~,order] = sort(groundtruth.imageN);
        [responses_sorted(session),groundtruth_sorted(session)] = sortData(responses,groundtruth,order);
    end
    GroundTruth = groundtruth_sorted(1);
    humanResponse(observer) = UniteSessions(responses_sorted,2);
end

humanResponse(end+1) = UniteObservers(humanResponse);
humanResponse_allobservers = humanResponse(end);

% calculate correlation across observers
% for target = 1:length(observerList)
%     cnt = 0;
%     target_val = mean(humanResponse(target).Pellacini_c,2);
% 
%     cnt = cnt + 1;
%     rest = mean(humanResponse(10).Pellacini_c,2);
%     corrCoeff(cnt) = corr(target_val,rest);
% 
%     meancorrCoeff(target) = mean(corrCoeff);
% end

lightprobecolor = brewermap(12,'Set3');
shapecolor = brewermap(18,'Set3');

%% compute correlation coefficients between human settings and groundtruth
for observer = 1:length(observerList)
    cnt = 0;
    for var = {'Pellacini_c'}
        cnt = cnt + 1;
        for envType = 1:3
            N = size(humanResponse(observer).(var{1}),2);
            
            % load human response and ground-truth values
            human_all_val = mean(humanResponse(11).(var{1})(3:end,:),2);

            gt_val = GroundTruth.(var{1})(3:end)';
            human_val = mean(humanResponse(observer).(var{1})(3:end,:),2);
            human_val_se = std(humanResponse(observer).(var{1})(3:end,:),[],2)/sqrt(N);
            
            % load human response and ground-truth values
            corrCoeff.(var{1})(observer) = corr(gt_val,human_val);
            corrCoeff_acrossParticipants.(var{1})(observer) = corr(human_val,human_all_val);
        end
    end
end

%% Compute Correlation across participants
for observer = 1:length(observerList)
    for var = {'Pellacini_c'}
        human_val = humanResponse(11).(var{1})(3:end,observer);
             
        human_rest_val = mean(humanResponse(11).(var{1})(3:end,[1:observer-1,observer+1:10]),2);
        human_all_val = mean(humanResponse(11).(var{1})(3:end,:),2);

        corrCoeff_imgStats.acrossParticipants_lb.(var{1})(observer) = corr(human_val,human_rest_val);
        corrCoeff_imgStats.acrossParticipants_ub.(var{1})(observer) = corr(human_val,human_all_val);
    end
end

%% Compute correlation between human settings and image statistics
load(fullfile(repo_basedir,'data','imgStats_exp2.mat'))

% this is just to create string list for stats
temp = fieldnames(imgStats);
for n = 1:9
    list_stats{n} = temp{n};
end

for var = {'Pellacini_c'}
    for observer = 1:length(observerList)
        for statsN = 1:length(list_stats)
            if ~strcmp(list_stats{statsN},'specularMetric')
                statsList = list_stats;
                
                gt_val = GroundTruth.(var{1})(3:end)';

                model_val = imgStats.(statsList{statsN})';
                human_val = humanResponse(11).(var{1})(3:end,observer);

                corrCoeff_imgStats.(var{1})(statsN,observer) = corr(model_val,human_val);
                pcorrCoeff_imgStats.(var{1})(statsN,observer) = partialcorr(human_val,model_val,gt_val);

            end
        end
    end
end

% sorting stats for the display purpose
desiredOrder = {'mean','percentile50','percentile25','percentile75','std','skewness','kurtosis','min','max'};
stats_label = {'cover.','sharp.','contrast','mean','median','Q1','Q3','std','skew.','kurtosis','min','max'};

for statsN = 1:length(desiredOrder)
    orderId(statsN) = find(strcmp(list_stats,desiredOrder(statsN)));
end


cnt = 0;
for var = {'Pellacini_c'}
    cnt = cnt + 1;
    for statsN = 1:length(desiredOrder)
        orderId(statsN) = find(strcmp([list_stats'],desiredOrder(statsN)));

        %orderId(statsN) = find(strcmp(['coverage';'sharpness';'contrast';list_stats'],desiredOrder(statsN)));
    end
end

%% Compute correlation with specular reflection metric (contrast, coverage, sharpness)
human = humanResponse(11);

cutoffList = [0 1 3 5 10 20 40];

gt_val = GroundTruth.Pellacini_c(3:end)';

for metric = {'contrast','coverage','sharpness'}
    model_val = [];
    cnt = 0;
    for cutoff = cutoffList   
        cnt = cnt + 1;

        switch metric{1}
            % sub-band contrast
            case 'contrast'
                for bandN = 1:9 % 1 - 8 are each subband and 9 indicates all summed
                    for observerN = 1:length(observerList)
                        model_val = imgStats.specularMetric.contrast.(['cutoff',num2str(cutoff)])(:,bandN);
                        modelval_record.(metric{1})(:,bandN,cnt) = model_val;

                        human_val = human.Pellacini_c(3:end,observerN);
                        corrCoeff_imgStats.(metric{1})(observerN,bandN,cnt) = corr(human_val,model_val);
                        pcorrCoeff_imgStats.(metric{1})(observerN,bandN,cnt) = partialcorr(human_val,model_val,gt_val);
                    end
                end

            % coverage
            case 'coverage'
                for observerN = 1:length(observerList)
                    model_val = imgStats.specularMetric.coverage.(['cutoff',num2str(cutoff)])';
                    modelval_record.(metric{1})(:,cnt) = model_val;

                    human_val = human.Pellacini_c(3:end,observerN);
                    corrCoeff_imgStats.(metric{1})(observerN,cnt) = corr(human_val,model_val);
                    pcorrCoeff_imgStats.(metric{1})(observerN,cnt) = partialcorr(human_val,model_val,gt_val);
                end

            % sharpness
            case 'sharpness'
                for observerN = 1:length(observerList)
                    model_val = imgStats.specularMetric.sharpness.(['cutoff',num2str(cutoff)])';
                    modelval_record.(metric{1})(:,cnt) = model_val;
                    human_val = human.Pellacini_c(3:end,observerN);
                    corrCoeff_imgStats.(metric{1})(observerN,cnt) = corr(human_val,model_val);
                    pcorrCoeff_imgStats.(metric{1})(observerN,cnt) = partialcorr(human_val,model_val,gt_val);
                end
        end
    end
end

% Compute mean across observers and find the best parameter that
% producess the highest correlation with human settings
for metric = {'contrast','coverage','sharpness'}
    if strcmp(metric{1},'contrast')
        meancorrCoeff_imgStats.(metric{1}) = squeeze(mean(corrCoeff_imgStats.contrast,2,'omitnan'));
    else
        meancorrCoeff_imgStats.(metric{1}) = squeeze(mean(corrCoeff_imgStats.(metric{1}),1,'omitnan'));
    end

    if strcmp(metric{1},'contrast')
        [~,id_temp] = max(abs(meancorrCoeff_imgStats.contrast(:)));
        [row, col] = ind2sub(size(meancorrCoeff_imgStats.contrast),id_temp);
        BestId.(metric{1})(1:2) = [row col];
    else
        [~,BestId.(metric{1})] = max(abs(meancorrCoeff_imgStats.(metric{1})(:)));
    end
end


%% Perform multiple regression using all predictors
contrast_best = modelval_record.contrast(:,BestId.contrast(1),BestId.contrast(2));    
coverage_best = modelval_record.coverage(:,BestId.coverage);    
sharpness_best = modelval_record.sharpness(:,BestId.sharpness);    

% cross-validation by leave-one-out
Model_all = [contrast_best coverage_best sharpness_best];
for stats = {'mean','percentile50','percentile25','percentile75','std','skewness','kurtosis','min','max'}
    Model_all = [Model_all,imgStats.(stats{1})'];
end

standarized_b = zeros(12,10);
corrCoeff_mregress = zeros(10,1);
for observerN = 1:length(observerList)
   targetObs = human.Pellacini_c(3:end,observerN); 
   restObs = mean(human.Pellacini_c(3:end,[1:observerN-1,observerN+1:end]),2); 

   [b,bint,r,rint,stats] = regress(restObs,[ones(size(Model_all,1),1),Model_all]);
   
   sy = std(restObs);
   sx = std([ones(size(Model_all,1),1),Model_all]);
   
   standarized_b(:,observerN) = b(2:end).*sx(2:end)'/sy;
   
   model_val = b(1)+Model_all*b(2:end);
   
   corrCoeff_mregress(observerN) = corr(model_val,targetObs);
end

%% Figure 14 - scatter plot human vs. groundturh for lighting environment
img_lightprobe = zeros(128,256,3,12);
for N = 1:12
    img_lightprobe(:,:,:,N) = imresize(imread(['./data/imgs_lightprobe/LightProbe',num2str(N),'.png']),[128 256]);
end

corrCoeff_lightprobe = zeros(12,11);
slope = zeros(12,11);

for observerN = 1:length(observerList)+1
    for lightprobeN = 1:12
        Id = find((GroundTruth.lightprobeN == lightprobeN));

        gt_val = GroundTruth.Pellacini_c(Id)';
        human_val = mean(humanResponse(observerN).Pellacini_c(Id,:),2);
        human_val_se = std(humanResponse(observerN).Pellacini_c(Id,:),[],2)/sqrt(10);
        corrCoeff_lightprobe(lightprobeN,observerN) = corr(gt_val,human_val);

        slope(lightprobeN,observerN) = 1/regress(gt_val,human_val);

        if observerN == length(observerList)+1
            if observerN == length(observerList)+1
                fig = figure;
                fig.Position = [0,10,1200,600];
            end

            l = refline(slope(lightprobeN,observerN));hold on

            l.Color = 'r';

            e = errorbar(gt_val,human_val,human_val_se,'CapSize',0);hold on;

            e.LineStyle = 'none';e.Color = [0 0 0];
            e.LineWidth = .5;

            line([0 400],[0 400],'Color',[0 0 0],'LineStyle',':','LineWidth',1)

            scatter(gt_val,human_val,30,([77 143 143]/255*0.8).^(1/2.2),'o','filled','MarkerEdgeColor',[.3 .3 .3],'LineWidth',0.3,'MarkerFaceAlpha',1);hold on

            axis square

            coeff = sprintf('%1.2f',mean(corrCoeff_lightprobe(lightprobeN,1:10),2));
            coeff_se = sprintf('%1.2f',std(corrCoeff_lightprobe(lightprobeN,1:10),[],2)/sqrt(10));

            [~,p] = corr(gt_val,human_val);
            if p > 0.1
                mark = [];
            elseif p < 0.01
                mark = '**';
            elseif p < 0.05
                mark = '*';
            end

            slope_temp = sprintf('%1.2f',1/regress(gt_val,human_val));

            t = text(0.03,0.95,[coeff,' ± ',coeff_se,' ',mark],'Color','k','Units','normalized');
            t2 = text(0.03,0.83,slope_temp,'Color','r','Units','normalized');

            t.FontSize = fontsize;t.FontName = 'Arial';
            t2.FontSize = fontsize;t2.FontName = 'Arial';

            ax = gca;
            fig.Units = 'centimeters';
            fig.Position = [10,10,twocolumn/6*0.95,twocolumn/6*0.95];
            fig.Color           = 'w';
            fig.InvertHardcopy  = 'off';

            fig.PaperPosition   = [0,10,6.45,6.45];

            ax.XLim = [minVal maxVal];
            ax.YLim = [minVal maxVal];

            xticks([0 0.24])
            yticks([0 0.24])

            ax.XTickLabel = {'0.00','0.24'};
            ax.YTickLabel = {'0.00','0.24'};

            %xlabel(['Ground Truth - ',variable{1}],'FontWeight', 'Bold');ylabel(['Human - ',variable{1}],'FontWeight', 'Bold');
            xlabel('','FontWeight', 'Bold');ylabel('','FontWeight', 'Bold');

            ax.FontName = 'Helvetica';
            ax.Color = 'w';
            ax.FontSize = fontsize;
            ax.XColor = 'k';ax.YColor = 'k';

            ax.LineWidth = 0.5;
            ax.Units = 'centimeters';
            ax.Position = [0.6 0.6 2.2 2.2];

            axis square

            grid off
            box off
            exportgraphics(fig,fullfile(pwd,'figs',['Fig14a_glossSet_lightprobe',num2str(lightprobeN),'.png']),'Resolution','600')
            close all
        end
    end
end

% unite all images
I = zeros(298,332,3,12);
L = zeros(138,332,3,12);

for lightprobeN = 1:12
    I(:,:,:,lightprobeN) = imresize(double(imread(fullfile(pwd,'figs',['Fig14a_glossSet_lightprobe',num2str(lightprobeN),'.png'])))/255,0.5);
    temp = double(imread(['./data/imgs_lightprobe/LightProbe',num2str(lightprobeN),'.png']))/255; % load image
    L(:,:,:,lightprobeN) = [ones(size(temp,1),size(I,2)-size(temp,2),3),temp];
end

wbar = ones(size(I,1),12,3);
halfbar = ones(size(I,1),6,3);

bar_L = ones(size(L,1),12,3);
halfbar_L = ones(size(L,1),6,3);

I1 = [halfbar,I(:,:,:,1),wbar,I(:,:,:,2),wbar,I(:,:,:,3),wbar,I(:,:,:,4),wbar,I(:,:,:,5),wbar,I(:,:,:,6),halfbar];
I2 = [halfbar,I(:,:,:,7),wbar,I(:,:,:,8),wbar,I(:,:,:,9),wbar,I(:,:,:,10),wbar,I(:,:,:,11),wbar,I(:,:,:,12),halfbar];

L1 = [halfbar_L,L(:,:,:,1)*1.2,bar_L,L(:,:,:,2)*1.2,bar_L,L(:,:,:,3)*1.4,bar_L,L(:,:,:,4)*1.4,bar_L,L(:,:,:,5)*1.2,bar_L,L(:,:,:,6)*1.2,halfbar_L];
L2 = [halfbar_L,L(:,:,:,7),bar_L,L(:,:,:,8)*1.5,bar_L,L(:,:,:,9)*1.2,bar_L,L(:,:,:,10)*1.2,bar_L,L(:,:,:,11)*1.4,bar_L,L(:,:,:,12)*1.2,halfbar_L];

out = [L1;ones(10,size(I1,2),3);I1;ones(50,size(I1,2),3);L2;ones(10,size(I1,2),3);I2];

imwrite(out,fullfile(pwd,'figs','Fig14a_combined.png'))

%% Fig14b - bar graph
corCoeff2 = corrCoeff_imgStats.Pellacini_c';

% get specular metric functions with best parameters
contrast = abs(corrCoeff_imgStats.contrast(:,BestId.contrast(1),BestId.contrast(2)));
coverage = abs(corrCoeff_imgStats.coverage(:,BestId.coverage));
sharpness = abs(corrCoeff_imgStats.sharpness(:,BestId.sharpness));

corCoeff_p = pcorrCoeff_imgStats.Pellacini_c';

%ub = mean(meancorrCoeff);
contrast_pcorrCoeff = abs(pcorrCoeff_imgStats.contrast(:,BestId.contrast(1),BestId.contrast(2)));
coverage_pcorrCoeff = abs(pcorrCoeff_imgStats.coverage(:,BestId.coverage));
sharpness_pcorrCoeff = abs(pcorrCoeff_imgStats.sharpness(:,BestId.sharpness));

mean_corCoeff_partial = [mean(coverage_pcorrCoeff) mean(sharpness_pcorrCoeff) mean(contrast_pcorrCoeff) abs(mean(corCoeff_p(:,orderId)))];
std_corCoeff_partial = [std(coverage_pcorrCoeff) std(sharpness_pcorrCoeff) std(contrast_pcorrCoeff) std(corCoeff_p(:,orderId))];

mean_corCoeff = [mean(coverage) mean(sharpness) mean(contrast) abs(mean(corCoeff2(:,orderId)))];
std_corCoeff = [std(coverage) std(sharpness) std(contrast) std(corCoeff2(:,orderId))];

fig = figure;

ub = corrCoeff_imgStats.acrossParticipants_ub.Pellacini_c;
lb = corrCoeff_imgStats.acrossParticipants_lb.Pellacini_c;

ub_gthuman_correlation = mean2(corrCoeff_lightprobe(:,1:10))'+std(mean(corrCoeff_lightprobe(:,1:10),2))/sqrt(10);
lb_gthuman_correlation = mean2(corrCoeff_lightprobe(:,1:10))'-std(mean(corrCoeff_lightprobe(:,1:10),2))/sqrt(10);

rectangle('Position',[0.03 lb_gthuman_correlation 100 ub_gthuman_correlation-lb_gthuman_correlation],'FaceColor',[220 230 240]/255,'EdgeColor','none');hold on;

rectangle('Position',[0.1 mean(lb) 100 mean(ub)-mean(lb)],'FaceColor',[224 195 249]/255,'EdgeColor','none');hold on;

all_mregress = corrCoeff_mregress;

r2 = refline(0,mean2(corrCoeff_lightprobe(:,1:10)));r2.Color = 'b';r2.LineStyle = '-';
b = bar(mean_corCoeff);hold on
b.FaceColor = 'flat';
b.EdgeColor = 'none';

for n = 1:length(stats_label)
    if n < 4
        b.CData(n,:) = [96 179 179]/255*0.8;
    else
        b.CData(n,:) = [149 179 179]/255;
    end
end

e = errorbar(1:length(desiredOrder)+3,mean_corCoeff,std_corCoeff,'LineStyle','none','CapSize',0);
e.Color = [0 0 0];e.LineWidth = .5;

scatter(1:length(desiredOrder)+3,mean_corCoeff_partial',20,[255 255 255]/255*0.7,'o','filled','MarkerEdgeColor','k');

r_all = line([-0.5 13],[mean(all_mregress) mean(all_mregress)]);r_all.Color = [1 0 0]*0.8;r_all.LineStyle = '--';

axis normal

ax = gca;
fig.Units           = 'centimeters';
fig.Position = [10,10,twocolumn/2,twocolumn/3*0.95];
fig.Color           = 'w';
fig.InvertHardcopy  = 'off';

ax.XLim = [0 length(desiredOrder)+4];
ax.YLim = [0 1.0];

xticks(1:length(stats_label))
yticks([0 0.5 1])
desiredOrder{2} = 'median';
ax.XTickLabel = stats_label;

ax.YTickLabel = {'0.0','0.5','1.0'};
ax.XTickLabelRotation = 45;

xlabel('Statistics','FontWeight', 'Bold');ylabel('Correlation to human settings','FontWeight', 'Bold');

ax.FontName = 'Arial';
ax.Color = 'w';
ax.FontSize = fontsize;
ax.XColor = 'k';ax.YColor = 'k';

ax.LineWidth = 0.5;
ax.Units = 'centimeters';
ax.Position = [0.8 1.17 8.4 4.6];
ticklengthcm(ax,0.0)
grid off
box off
exportgraphics(fig,fullfile(pwd,'figs','Fig14b_bar.pdf'),'ContentType','vector')
close all

%% Figure14b - standarized coefficient

order_coeff = 4:12;
order_coeff =  [2 3 1 order_coeff];

fig = figure;
s_coeff = standarized_b(order_coeff,:);
b = bar(mean(s_coeff,2));hold on
b.FaceColor = 'flat';
b.EdgeColor = 'none';

for n = 1:length(stats_label)
    if n < 4
        b.CData(n,:) = [96 179 179]/255*0.8;
    else
        b.CData(n,:) = [149 179 179]/255;
    end
end

e = errorbar(1:length(desiredOrder)+3,mean(s_coeff,2),std(s_coeff,[],2),'LineStyle','none','CapSize',0);
e.Color = [0 0 0];e.LineWidth = .5;

axis normal

ax = gca;
fig.Units           = 'centimeters';
fig.Position = [10,10,twocolumn/2.5,twocolumn/3*0.95];
fig.Color           = 'w';
fig.InvertHardcopy  = 'off';

ax.XLim = [0 length(desiredOrder)+4];
ax.YLim = [-0.6 0.6];

xticks(1:length(stats_label))
yticks([-0.6 0 0.6])
desiredOrder{2} = 'median';
ax.XTickLabel = stats_label;

ax.YTickLabel = {'-0.6','0.0','0.6'};
ax.XTickLabelRotation = 45;

%xlabel(['Ground Truth - ',variable{1}],'FontWeight', 'Bold');ylabel(['Human - ',variable{1}],'FontWeight', 'Bold');
xlabel('Statistics','FontWeight', 'Bold');ylabel('Standardized weights','FontWeight', 'Bold');

ax.FontName = 'Arial';
ax.Color = 'w';
ax.FontSize = fontsize;
ax.XColor = 'k';ax.YColor = 'k';

ax.LineWidth = 0.5;
ax.Units = 'centimeters';
ax.Position = [0.85 1.2 6.4 4.55];
ticklengthcm(ax,0.0)
grid off
box off

exportgraphics(fig,fullfile('figs','Fig14b_weights.pdf'),'ContentType','vector')
close all

%% Custom function to sort human response data (as the presentation order was randomized in the experiment)
function [responses_sorted, groundtruth_sorted] = sortData(RESPONSES,GROUNDTRUTH,order)
    % function to order all field based on order
    responses_sorted.specularity = RESPONSES.specularity(order);
    responses_sorted.Pellacini_c = RESPONSES.Pellacini_c(order);

    groundtruth_sorted.specularity = GROUNDTRUTH.specularity(order);
    groundtruth_sorted.Pellacini_c = GROUNDTRUTH.c_Pellacini(order);
    groundtruth_sorted.imageN = GROUNDTRUTH.imageN(order);
    groundtruth_sorted.lightprobeN = GROUNDTRUTH.lightprobeN(order);
    groundtruth_sorted.objN = GROUNDTRUTH.objN(order);

end

%% Custom function to sort human response data acorss sessions
function humanResponse = UniteSessions(RESPONSES,sessionN)
    humanResponse.specularity = [];
    humanResponse.Pellacini_c = [];
    
    for session = 1:sessionN
        humanResponse.specularity = [humanResponse.specularity,RESPONSES(session).specularity'];
        humanResponse.Pellacini_c = [humanResponse.Pellacini_c,RESPONSES(session).Pellacini_c'];
    end
end

%% Custom function to sort response data across all human observers
function humanResponse_allobservers = UniteObservers(humanResponse)
    Outputs = [size(humanResponse(1).specularity),length(humanResponse)];
    humanResponse_allobservers.specularity = [];
    humanResponse_allobservers.Pellacini_c = [];
    
    for observer = 1:length(humanResponse)
        humanResponse_allobservers.specularity = [humanResponse_allobservers.specularity,humanResponse(observer).specularity];
        humanResponse_allobservers.Pellacini_c = [humanResponse_allobservers.Pellacini_c,humanResponse(observer).Pellacini_c];
    end
    
    % take average across sessions
    humanResponse_allobservers.specularity = squeeze(mean(reshape(humanResponse_allobservers.specularity,Outputs),2));
    humanResponse_allobservers.Pellacini_c = squeeze(mean(reshape(humanResponse_allobservers.Pellacini_c,Outputs),2));
end

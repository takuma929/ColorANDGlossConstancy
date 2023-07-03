clear;close all;clc; % cleaning up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Script to reproduce all figures for Exp1 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(pwd,'data','CGC_FigParameters')) % fig parameters load
load(fullfile(repo_basedir,'data','imgStats_exp1')) % load pre-computed image-statistics for all images

% specify which figures to Plot
PlotFig6=1;
PlotFig7AND10=1;
PlotFig8=1;
PlotFig9=1;
PlotFig12=1;

% Set value range for lightness, hue, chroma, Pellacini's c, and saturation
% (C*/L*ab). These values are used to set axis ranges for figures.
maxVal = [74 368 26 0.1487 0.75];
minVal = [26 -8 8 0 0.1];

% Indices for each lighting environment
% Note: indices 1 and 2 are control images (symmetric matching).
env(1).Id = 3:3:38; % indices for natural lighting environment
env(2).Id = 4:3:38; % indices for gamut-rotated lighting environment
env(3).Id = 5:3:38; % indices for phase-scrambled lighting environment

% list of observers
observerList = {'SB','KK','JM','MF','ET','IW','TG','AM','OP','CY'};

% label for each lighting environment
envType_label = {'natural','gamutrotation','phasescrambled'};

% list of vars
variableList = {'lightness','hue','chroma','Pellacini_c'};

% number of sessions (three sessions in Experiment 1)
sessionN = 3;  

% RGB values for natural, gamut-rotated, phase-scrambled environments
envType_RGB = [205 151 149;179 192 145;143 184 196]/255;
envType_RGB_human = [205 188 187;188 192 177;179 192 195]/255;

% load human response data and associated groundtruth 
for observer = 1:length(observerList)
    clear responses_sorted groundtruth_sorted
    for session = 1:sessionN
        load(fullfile(repo_basedir,'data','rawdata','exp1',[observerList{observer},'_session',num2str(session),'.mat']))
        [~,order] = sort(groundtruth.imageN);
        [responses_sorted(session),groundtruth_sorted] = sortData(responses,groundtruth,order);
    end
    
    % groundtruth are same regardless of observers and sessions so get load one 
    if observer == 1
        groundTruth = groundtruth_sorted;
    end
    humanResponse(observer) = UniteSessions(responses_sorted,3);
end

humanResponse(end+1) = UniteObservers(humanResponse); % concatenate　all observer data
humanResponse_allobservers = humanResponse(end);

%% compute correlation coefficients between human settings and groundtruth
for observer = 1:length(observerList)
    cnt = 0;
    for var = {'lightness','hue','chroma','Pellacini_c','saturation'}
        cnt = cnt + 1;
        for envType = 1:3
            N = size(humanResponse(observer).(var{1}),2);
            
            % load human response and ground-truth values
            if strcmp(var{1},'specularity')||strcmp(var{1},'lightness')||strcmp(var{1},'chroma')||strcmp(var{1},'Pellacini_c')||strcmp(var{1},'saturation')
                human_all_val = mean(humanResponse(11).(var{1})(env(envType).Id,:),2);

                gt_val = groundTruth.(var{1})(env(envType).Id)';
                human_val = mean(humanResponse(observer).(var{1})(env(envType).Id,:),2);
                human_val_se = std(humanResponse(observer).(var{1})(env(envType).Id,:),[],2)/sqrt(N);
            elseif strcmp(var{1},'hue')
                human_all_val = circ_mean(humanResponse(11).(var{1})(env(envType).Id,:)/180*pi,[],2);

                gt_val = groundTruth.(var{1})(env(envType).Id)';
                human_val = mod(circ_mean(humanResponse(observer).(var{1})(env(envType).Id,:)/180*pi,[],2)/pi*180+360,360);
                human_val_se = mod(real(circ_std(humanResponse(observer).(var{1})(env(envType).Id,:)/180*pi,[],[],2)/pi*180+360),360)/sqrt(N);
            end
            
            % load human response and ground-truth values
            if strcmp(var{1},'hue')
                corrCoeff.(var{1})(envType,observer) = circ_corrcc(deg2rad(gt_val),deg2rad(human_val));
                corrCoeff_acrossParticipants.(var{1})(envType,observer) = circ_corrcc(human_val/180*pi,human_all_val);
            else
                corrCoeff.(var{1})(envType,observer) = corr(gt_val,human_val);
                corrCoeff_acrossParticipants.(var{1})(envType,observer) = corr(human_val,human_all_val);
            end
        end
    end
end

%% Compute Correlation across participants
for observer = 1:length(observerList)
    for envType = 1:3
        for var = variableList
            human_val = humanResponse(11).(var{1})(env(envType).Id,observer);

            if strcmp(var{1},'hue')
                human_rest_val = circ_mean(humanResponse(11).(var{1})(env(envType).Id,[1:observer-1,observer+1:10])/180*pi,[],2);
                human_all_val = circ_mean(humanResponse(11).(var{1})(env(envType).Id,:)/180*pi,[],2);
            
                corrCoeff_imgStats.acrossParticipants_lb.(var{1})(envType,observer) = circ_corrcc(human_val/180*pi,human_rest_val);
                corrCoeff_imgStats.acrossParticipants_ub.(var{1})(envType,observer) = circ_corrcc(human_val/180*pi,human_all_val);
            else                
                human_rest_val = mean(humanResponse(11).(var{1})(env(envType).Id,[1:observer-1,observer+1:10]),2);
                human_all_val = mean(humanResponse(11).(var{1})(env(envType).Id,:),2);
                
                corrCoeff_imgStats.acrossParticipants_lb.(var{1})(envType,observer) = corr(human_val,human_rest_val);
                corrCoeff_imgStats.acrossParticipants_ub.(var{1})(envType,observer) = corr(human_val,human_all_val);
            end
        end
    end
end

%% Compute correlation between human settings and image statistics
load(fullfile(repo_basedir,'data','imgStats_exp1.mat'))

% this is just to create string list for stats
list_stats.lightness = fieldnames(imgStats.diffuseOnly.lightness);
list_stats.hue = fieldnames(imgStats.diffuseOnly.hue);
list_stats.chroma = fieldnames(imgStats.diffuseOnly.chroma);
list_stats.Pellacini_c = fieldnames(imgStats.diffuseOnly.Pellacini_c);

for type = {'diffuseANDspecular','diffuseOnly'}
    for envType = 1:3
        for var = variableList
            for observer = 1:length(observerList)
                for statsN = 1:length(list_stats.(var{1}))
                    statsList = list_stats.(var{1});
                    
                    model_val = imgStats.(type{1}).(var{1}).(statsList{statsN})(env(envType).Id-2)';
                    human_val = humanResponse(11).(var{1})(env(envType).Id,observer);
                    
                    if strcmp(var{1},'hue')
                        corrCoeff_imgStats.human.(type{1}).(var{1})(envType,statsN,observer) = circ_corrcc(deg2rad(mod(model_val+360,360)),deg2rad(mod(human_val+360,360)));
                    else
                        corrCoeff_imgStats.human.(type{1}).(var{1})(envType,statsN,observer) = corr(model_val,human_val);
                    end
                end
            end
        end
    end
end

% sorting stats for the display purpose
for var = {'lightness','chroma','Pellacini_c'}
    desiredOrder.(var{1}) = {'mean','percentile50','percentile25','percentile75','std','skewness','kurtosis','min','max'};
    stats_label.(var{1}) = {'mean','median','Q1','Q3','std','skew.','kurtosis','min','max'};
end

desiredOrder.hue = {'mean','percentile50','percentile25','percentile75','std','skewness','kurtosis'};
stats_label.hue = {'mean','median','Q1','Q3','std','skew.','kurtosis'};

desiredOrder.Pellacini_c = {'coverage','sharpness','contrast','mean','percentile50','percentile25','percentile75','std','skewness','kurtosis','min','max'};
stats_label.Pellacini_c = {'cover.','sharp.','contrast','mean','median','Q1','Q3','std','skew.','kurtosis','min','max'};

cnt = 0;
for var = {'lightness','hue','chroma','Pellacini_c'}
    cnt = cnt + 1;
    for statsN = 1:length(desiredOrder.(var{1}))
        if strcmp(var{1},'Pellacini_c')
            orderId.(var{1})(statsN) = find(strcmp(['coverage';'sharpness';'contrast';list_stats.(variableList{cnt})],desiredOrder.(var{1})(statsN)));
        else
            orderId.(var{1})(statsN) = find(strcmp(list_stats.(variableList{cnt}),desiredOrder.(var{1})(statsN)));
        end
    end
end

%% Compute correlation with specular reflection metric (contrast, coverage, sharpness)
human = humanResponse(11);

cutoffList = [0 1 3 5 10 20 40];

for envType = 1:3
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
                            model_val = imgStats.specularMetric.contrast.(['cutoff',num2str(cutoff)])(env(envType).Id-2,bandN);
                            human_val = human.Pellacini_c(env(envType).Id,observerN);
                            corrCoeff_imgStats.human.(metric{1})(envType,observerN,bandN,cnt) = corr(human_val,model_val);
                        end
                    end

                % coverage
                case 'coverage'
                    for observerN = 1:length(observerList)
                        model_val = imgStats.specularMetric.coverage.(['cutoff',num2str(cutoff)])(env(envType).Id-2)';
                        human_val = human.Pellacini_c(env(envType).Id,observerN);
                        corrCoeff_imgStats.human.(metric{1})(envType,observerN,cnt) = corr(human_val,model_val);
                    end

                % sharpness
                case 'sharpness'
                    for observerN = 1:length(observerList)
                        model_val = imgStats.specularMetric.sharpness.(['cutoff',num2str(cutoff)])(env(envType).Id-2)';
                        human_val = human.Pellacini_c(env(envType).Id,observerN);
                        corrCoeff_imgStats.human.(metric{1})(envType,observerN,cnt) = corr(human_val,model_val);
                    end
            end
        end
    end
end

% Compute mean across observers and find the best parameter that
% producess the highest correlation with human settings
for metric = {'contrast','coverage','sharpness'}
    for envType = 1:3
        if strcmp(metric{1},'contrast')
            meancorrCoeff_imgStats.(metric{1}) = squeeze(mean(corrCoeff_imgStats.human.contrast(envType,:,:,:),2));
        else
            meancorrCoeff_imgStats.(metric{1}) = squeeze(mean(corrCoeff_imgStats.human.(metric{1})(envType,:,:),2));
        end
        
        if strcmp(metric{1},'contrast')
            [~,id_temp] = max(abs(meancorrCoeff_imgStats.contrast(:)));
            [row, col] = ind2sub(size(meancorrCoeff_imgStats.contrast),id_temp);
            BestId.(metric{1})(envType,1:2) = [row col];
        else
            [~,BestId.(metric{1})(envType)] = max(abs(meancorrCoeff_imgStats.(metric{1})(:)));
        end
    end
end

%% Figure 6 - results for control condition
if PlotFig6
disp('Generating Figure 6...')

maxVal_control = [74 360 26 0.1487 0.75];
minVal_control = [26 0 8 0 0.1];

label = {'{\it L^*}','{\it Hue [deg]}','{\it C^*_a_b}','{\it Pellacini''s c}'};

for observer = length(observerList)+1
    cnt = 0;
    for var = {'lightness','hue','chroma','Pellacini_c'}
        fig = figure;

        cnt = cnt + 1;

        N = size(humanResponse(observer).(var{1}),2);

        if strcmp(var{1},'specularity')||strcmp(var{1},'lightness')||strcmp(var{1},'chroma')||strcmp(var{1},'Pellacini_c')||strcmp(var{1},'saturation')
            gt_val = groundTruth.(var{1})(1:2)';
            human_val = humanResponse(observer).(var{1})(1:2,:);
            human_val_se = std(humanResponse(observer).(var{1})(1:2,:),[],2)/sqrt(N);
        elseif strcmp(var{1},'hue')
            gt_val = groundTruth.(var{1})(1:2)';
            human_val = mod(humanResponse(observer).(var{1})(1:2,:)/180*pi/pi*180+360,360);
            human_val_se = mod(real(circ_std(humanResponse(observer).(var{1})(1:2,:)/180*pi,[],[],2)/pi*180+360),360)/sqrt(N);
        end

        control = 1; % plot only results for control image 1
        
        % plot error bar
        e = errorbar(1,mean(human_val(1,:),2),human_val_se(1)*sqrt(10),human_val_se(control)*sqrt(10),'LineStyle','none','CapSize',0);
        e.Color = [.5 .5 .5];e.LineWidth = 1;

        % plot line for ground-truth value
        line([0.6,1.4],[gt_val(control),gt_val(control)],'Color','b');hold on
        
        % record precision value
        precision.(var{1})(control,:) = abs(human_val(control,:)-gt_val(control));
        
        % plot human settings
        scatter(control,mean(human_val(control,:),2),40,([64 71 29]/255).^(1/2.2),'o','filled','MarkerEdgeColor', [0 0 0],'LineWidth',0.5,'MarkerFaceAlpha',1);hold on
        
        % make the figure look nicer
        ax = gca;
        fig.Units = 'centimeters';
        fig.Color = 'w';
        fig.InvertHardcopy = 'off';
        fig.PaperPosition   = [0,10,8.45,8.45];
        fig.Position = [10,10,twocolumn/8*0.6,twocolumn/4*0.95];

        ax.XLim = [0.5 1.5];
        ax.YLim = [minVal_control(cnt) maxVal_control(cnt)];

        xticks([1 2])
        yticks([minVal_control(cnt) maxVal_control(cnt)])

        ax.XTickLabel = {num2str(1) num2str(2)};
        ax.YTickLabel = {num2str(round(minVal_control(cnt)*100)/100) num2str(round(maxVal_control(cnt)*100)/100)};                

        xlabel('','FontWeight', 'Bold');ylabel(label{cnt},'FontWeight', 'Bold');

        ax.FontName = 'Arial';
        ax.Color = [.97 .97 .97];
        ax.FontSize = fontsize;
        ax.XColor = 'k';ax.YColor = 'k';

        ax.LineWidth = 1;
        ax.Units = 'centimeters';
        ax.Position = [0.9 0.3 1.2 4.0];

        grid off
        box off
        
        % save figure
        exportgraphics(fig,fullfile(repo_basedir,'figs',['Fig6_',var{1},'.pdf']),'ContentType','vector')
        close all
    end
end
disp('Done.')
end

%% Figures 7 and 10 - scatter plots for human vs. ground-truth
if PlotFig7AND10
disp('Generating Figure 7 and 10...')

% set numbers not to draw in the plot (just for the sake of visiblity)
notplotNumbers.hue.envType1 = [2,3,5]; 
notplotNumbers.lightness.envType1 = [2,3];
notplotNumbers.chroma.envType1 = [];
notplotNumbers.Pellacini_c.envType1 = 10;
notplotNumbers.saturation.envType1 = [2 9];

notplotNumbers.hue.envType2 = [1,5]; 
notplotNumbers.lightness.envType2 = [];
notplotNumbers.chroma.envType2 = 7;
notplotNumbers.Pellacini_c.envType2 = 8;
notplotNumbers.saturation.envType2 = [];

notplotNumbers.hue.envType3 = 7;
notplotNumbers.lightness.envType3 = [];
notplotNumbers.chroma.envType3 = [];
notplotNumbers.Pellacini_c.envType3 = 5;
notplotNumbers.saturation.envType3 = [5 6 7 12];

for observer = length(observerList)+1 % this plots only average across observers 
    for envType = 1:3
        cnt = 0; % loop counter for var
        for var = {'lightness','hue','chroma','Pellacini_c','saturation'}
            cnt = cnt + 1;
            
            N = size(humanResponse(observer).(var{1}),2);
            
            if strcmp(var{1},'specularity')||strcmp(var{1},'lightness')||strcmp(var{1},'chroma')||strcmp(var{1},'Pellacini_c')||strcmp(var{1},'saturation')
                gt_val = groundTruth.(var{1})(env(envType).Id)';
                human_val = mean(humanResponse(observer).(var{1})(env(envType).Id,:),2);
                human_val_se = std(humanResponse(observer).(var{1})(env(envType).Id,:),[],2)/sqrt(N);
            elseif strcmp(var{1},'hue')
                gt_val = groundTruth.(var{1})(env(envType).Id)';
                human_val = mod(circ_mean(humanResponse(observer).(var{1})(env(envType).Id,:)/180*pi,[],2)/pi*180+360,360);
                human_val_se = mod(real(circ_std(humanResponse(observer).(var{1})(env(envType).Id,:)/180*pi,[],[],2)/pi*180+360),360)/sqrt(N);
            end
            
            % Subtract or add 360 degrees to improve the clarity of the
            if envType == 2 && strcmp(var{1},'hue')
                human_val(3) = human_val(3)+360;
            elseif envType == 3 && strcmp(var{1},'hue')
                human_val(7) = human_val(7)-360;
                human_val(11) = human_val(11)-360;
            end
            
            % Plot figures and make it look nice
            fig = figure;
            
            % error bar
            e = errorbar(gt_val,human_val,human_val_se,'CapSize',0);hold on;
            e.LineStyle = 'none';e.Color = [0.0 0.0 0.0];
            e.LineWidth = 1.0;
            
            % unity line
            line([0 400],[0 400],'Color',[0 0 0],'LineStyle',':','LineWidth',1)
           
            % plot data points
            % imgN < 7 is outdoorm imgN >= 7 is indoor
            for imgN = 1:12
                if imgN < 7
                    % circle symbol for outdoor
                    scatter(gt_val(imgN),human_val(imgN),80,envType_RGB(envType,:),'o','filled','MarkerEdgeColor', [0 0 0],'LineWidth',0.5,'MarkerFaceAlpha',1);hold on
                else
                    % diamond symbol for indoor
                    scatter(gt_val(imgN),human_val(imgN),80,envType_RGB(envType,:),'d','filled','MarkerEdgeColor', [0 0 0],'LineWidth',0.5,'MarkerFaceAlpha',1);hold on
                end
            end            
            plotNums = setdiff(1:12,notplotNumbers.(var{1}).(['envType',num2str(envType)]));
            
            % plot image number for each data point
            for imageN = plotNums
                text(gt_val(imageN),human_val(imageN),num2str(imageN),'FontSize',6,'FontName','Arial','Color','w','HorizontalAlignment','center')
            end    
            
            % get string for coeff and se
            coeff = sprintf('%1.2f',mean(corrCoeff.(var{1})(envType,:)));
            coeff_se = sprintf('%1.2f',std(corrCoeff.(var{1})(envType,:))/sqrt(10));
            if ~strcmp(var,'saturation')
                coeff_ub = sprintf('%1.2f',mean(corrCoeff_imgStats.acrossParticipants_ub.(variableList{cnt})(envType,:)));
                coeff_ub_se = sprintf('%1.2f',std(corrCoeff_imgStats.acrossParticipants_ub.(variableList{cnt})(envType,:))/sqrt(10));
            end

            % show whether the correlation coeff. was significant
            [~,p] = corr(gt_val,human_val);
            if p > 0.1
                mark = [];
            elseif p < 0.01
                mark = '**';
            elseif p < 0.05
                mark = '*';
            end
            
            % show correlaiton coefficients
            if strcmp(var{1},'chroma') && envType == 2
                t = text(0.025,0.90,[coeff,' ± ',coeff_se,' ',mark],'Color','k','Units','normalized');
            elseif strcmp(var{1},'lightness') && envType == 1
                t = text(0.025,0.85,[coeff,' ± ',coeff_se,' ',mark],'Color','k','Units','normalized');
            else
                t = text(0.025,0.94,[coeff,' ± ',coeff_se,' ',mark],'Color','k','Units','normalized');   
            end
            
            t.FontSize = fontsize;t.FontName = 'Arial';
            axis square

            ax = gca;
            fig.Units = 'centimeters';
            fig.Color = 'w';
            fig.InvertHardcopy = 'off';
            fig.PaperPosition   = [0,10,8.45,8.45];
            fig.Position = [10,10,twocolumn/4*0.95,twocolumn/4*0.95];
            
            % add axis label and limits
            if strcmp(var{1},'hue')
                if envType == 3
                    ax.XLim = [-40 400];
                    ax.YLim = [-40 400];
                else
                    ax.XLim = [minVal(cnt) maxVal(cnt)];
                    ax.YLim = [minVal(cnt) maxVal(cnt)];
                end

                xticks([0 360])
                yticks([0 360])

                ax.XTickLabel = {'0','360'};
                ax.YTickLabel = {'0','360'};
            else
                ax.XLim = [minVal(cnt) maxVal(cnt)];
                ax.YLim = [minVal(cnt) maxVal(cnt)];

                xticks([minVal(cnt) maxVal(cnt)])
                yticks([minVal(cnt) maxVal(cnt)])

                ax.XTickLabel = {num2str(round(minVal(cnt)*100)/100) num2str(round(maxVal(cnt)*100)/100)};
                ax.YTickLabel = {num2str(round(minVal(cnt)*100)/100) num2str(round(maxVal(cnt)*100)/100)};                
            end
            
            % define a position to add precision line
            ypos = ax.YLim(1)+(ax.YLim(2)-ax.YLim(1))/2;
            xpos = ax.XLim(1)+(ax.XLim(2)-ax.XLim(1))*0.05;

            % add precision line
            if ~strcmp(var{1},'saturation')
                if strcmp(var{1},'hue')
                    pr = circ_mean(circ_mean(precision.(var{1})/180*pi,[],2))/pi*180;
                else
                    pr = mean(mean(precision.(var{1}),2));
                end
                line([xpos,xpos],[ypos-pr/2 ypos+pr/2],'Color',([60 60 255]/255).^(1/2.2),'LineWidth',2)
            end
            
            xlabel('','FontWeight', 'Bold');ylabel('','FontWeight', 'Bold');

            ax.FontName = 'Arial';
            ax.Color = [.97 .97 .97];
            ax.FontSize = fontsize;
            ax.XColor = 'k';ax.YColor = 'k';

            ax.LineWidth = 1;
            ax.Units = 'centimeters';
            ax.Position = [0.62 0.6 3.5 3.5];

            grid off
            box off
            
            % save figures
            if strcmp(var{1},'saturation')
                exportgraphics(fig,fullfile(repo_basedir,'figs',['Fig10_',var{1},'_',envType_label{envType},'.pdf']),'ContentType','vector')
            else
                exportgraphics(fig,fullfile(repo_basedir,'figs',['Fig7_',var{1},'_',envType_label{envType},'.pdf']),'ContentType','vector')
            end
            close all
        end
    end
end
disp('Done.')
end

%% Figure 8 - human-groundtruth correlation and human-human correlation 
if PlotFig8
disp('Generating Figure 8...')

% define the offset to avoid overlap
offset = linspace(-0.2,0.2,10)';

% plot figures
for envType = 1:3
    cnt = 0;
    fig = figure;

    for var = {'hue','lightness','chroma','Pellacini_c'}
        cnt = cnt + 1;
        
        cor_gtvshuman = corrCoeff.(var{1})(envType,:);
        cor_humanvshuman = corrCoeff_acrossParticipants.(var{1})(envType,:);

        x1 = (cnt-1)*2+1+(cnt-0.5);
        x2 = x1+1;
        
        % show bar graphs for two types of correlation
        b = bar(x1,mean(cor_gtvshuman),'FaceColor',envType_RGB(envType,:),'EdgeColor','none');hold on
        bar(x2,mean(cor_humanvshuman),'FaceColor',envType_RGB_human(envType,:),'EdgeColor','none');hold on

        for observerN = 1:10
            line([x1+offset(observerN),x2+offset(observerN)],[cor_gtvshuman(observerN)',cor_humanvshuman(observerN)'],'Color',[.5 .5 .5],'LineWidth',.3)
        end
        scatter(ones(10,1)*x1+offset,cor_gtvshuman,20,envType_RGB(envType,:),'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[.3 .3 .3],'LineWidth',.5);
        scatter(ones(10,1)*x2+offset,cor_humanvshuman,20,envType_RGB_human(envType,:),'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[.3 .3 .3],'LineWidth',.5);
    end

    ax = gca;
    fig.Units = 'centimeters';
    fig.Color = 'w';
    fig.InvertHardcopy = 'off';
    fig.PaperPosition   = [0,10,8.45,8.45];
    fig.Position = [10,10,twocolumn/2,twocolumn/5];

    ax.XLim = [0.5 12.5];ax.YLim = [-0.0 1.05];
    xticks([2 5 8 11])
    yticks([0 0.5 1.0])

    ax.XTickLabel = {'hue','lightness','chroma','Pellacini''s c'};
    ax.YTickLabel = {'0.0','0.5','1.0'};

    ylabel('correlation coefficient','FontWeight', 'Bold');

    ax.FontName = 'Arial';
    ax.Color = [.97 .97 .97];
    ax.FontSize = fontsize+1;
    ax.XColor = 'k';ax.YColor = 'k';

    ax.LineWidth = 0.5;
    ax.Units = 'centimeters';
    ax.Position = [0.85 0.4 8.3 3.2];
    ticklengthcm(ax,0)

    grid on
    box off
    ax.XGrid = 'off';

    exportgraphics(fig,fullfile(repo_basedir,'figs',['Fig8_',envType_label{envType},'.pdf']),'ContentType','vector')
end
close all
disp('Done.')
end

%% Figure 9 - interaction between diffuse reflectance and specular reflectance
if PlotFig9
disp('Generating Figure 9...')

clear p

% a*b* coordinates for color temperature from 4000 to 10000 with 500 step
% (black body locus)
ab_bbl = [7.87 35.23;4.45 25.80;2.32 17.78;1.02 10.96;0.28 5.14;...
        -0.11 0.15;-0.26 -4.17;-0.25 -7.92;-0.15 -11.20;0.02 -14.09;0.22 -16.64;0.45 -18.91;0.68 -20.94];

rgb_plot.envType1 = [0.42,0.44,0.43,0.22,0.45,0.56,0.49,0.54,0.49,0.57,0.49,0.32;0.24,0.40,0.58,0.33,0.36,0.52,0.67,0.41,0.64,0.47,0.61,0.39;0.27,0.29,0.56,0.40,0.48,0.45,0.55,0.51,0.65,0.62,0.60,0.24];
rgb_plot.envType2 = [0.67,0.51,0.45,0.67,0.59,0.35,0.55,0.26,0.53,0.54,0.66,0.43;0.51,0.47,0.36,0.48,0.68,0.43,0.48,0.36,0.39,0.58,0.63,0.65;0.44,0.58,0.42,0.43,0.51,0.32,0.61,0.49,0.38,0.46,0.53,0.63];
rgb_plot.envType3 = [0.62,0.25,0.23,0.37,0.59,0.59,0.62,0.59,0.60,0.22,0.74,0.65;0.50,0.53,0.49,0.43,0.63,0.60,0.46,0.54,0.41,0.40,0.62,0.55;0.48,0.55,0.45,0.32,0.73,0.74,0.50,0.47,0.49,0.44,0.63,0.69];

for envType = 1:3
    fig = figure;
    
    chroma = groundTruth.chroma(env(envType).Id);
    hue = groundTruth.hue(env(envType).Id);
    lightness = groundTruth.lightness(env(envType).Id);
    gloss = groundTruth.Pellacini_c(env(envType).Id);
    
    a = cos(hue/180*pi).*chroma;
    b = sin(hue/180*pi).*chroma;
    
    % get human gloss setting
    human_gloss = mean(humanResponse_allobservers.Pellacini_c(env(envType).Id,:),2);
    
    % compute error between ground-truth and human_gloss
    error = gloss' - human_gloss;
    
    axis square
    
    % plot black body locus
    plot(ab_bbl(:,1),ab_bbl(:,2),'Color','k','LineWidth',0.5);hold on;
    
    % plot 36 images - change the circle size according to the error
    for imageN = 1:length(chroma)
        rgb = lab2rgb(lightness(imageN),a(imageN),b(imageN));
        rgb_record.(['envType',num2str(envType)])(:,imageN) = rgb;
        if error(imageN) > 0
            scatter(a(imageN),b(imageN),abs(round(error(imageN)*5000)),rgb_plot.(['envType',num2str(envType)])(:,imageN)','o','filled','MarkerEdgeColor','r','MarkerFaceAlpha',.5);hold on;
        else
            scatter(a(imageN),b(imageN),abs(round(error(imageN)*5000)),rgb_plot.(['envType',num2str(envType)])(:,imageN)','o','filled','MarkerEdgeColor','b','MarkerFaceAlpha',.5);hold on;
        end
    end
    
    % make the figure look nicer
    ax = gca;
    fig.Units = 'centimeters';
    fig.Color = 'w';
    fig.InvertHardcopy = 'off';
    fig.PaperPosition   = [0,10,8.45,8.45];
    fig.Position = [10,10,twocolumn/4*0.95,twocolumn/4*0.95];

    ax.XLim = [-30 30];
    ax.YLim = [-30 30];

    xticks([-30 0 30])
    yticks([-30 0 30])

    ax.XTickLabel = ["-30","0","+30"];
    ax.YTickLabel = ["-30","0","+30"];

    xlabel('a*','FontWeight', 'Bold');ylabel('b*','FontWeight', 'Bold');

    ax.FontName = 'Arial';
    ax.Color = [.97 .97 .97];
    ax.FontSize = fontsize;
    ax.XColor = 'k';ax.YColor = 'k';

    ax.LineWidth = 1;
    ax.Units = 'centimeters';
    ax.Position = [0.85 0.8 3.25 3.25];

    grid on
    box off
    
    % save image
    exportgraphics(fig,fullfile(repo_basedir,'figs',['Fig9_',envType_label{envType},'.pdf']),'ContentType','vector')
    close all
end

for envType = 1:3
    for imageN = 1:36
        %lightness = groundTruth.lightness(env(envType).Id);
        %gloss = groundTruth.Pellacini_c(env(envType).Id);
    
        %human_gloss = mean(humanResponse_allobservers.Pellacini_c(env(envType).Id,:),2);
    end
end

%% Figure 9 -  interaction between diffuse reflectance and specular reflectance (lightness direction)
for envType = 1:3
    fig = figure;
    
    lightness = groundTruth.lightness(env(envType).Id);
    gloss = groundTruth.Pellacini_c(env(envType).Id);
    
    human_gloss = mean(humanResponse_allobservers.Pellacini_c(env(envType).Id,:),2);
    error = gloss' - human_gloss;
        
    for imageN = 1:length(chroma)
        if error(imageN) > 0
            scatter(72-lightness(imageN)+27,0,abs(round(error(imageN)*5000)),rgb_plot.(['envType',num2str(envType)])(:,imageN)','o','filled','MarkerEdgeColor','r','MarkerFaceAlpha',.5);hold on;
        else
            scatter(72-lightness(imageN)+27,0,abs(round(error(imageN)*5000)),rgb_plot.(['envType',num2str(envType)])(:,imageN)','o','filled','MarkerEdgeColor','b','MarkerFaceAlpha',.5);hold on;
        end
    end
    
    ax = gca;
    fig.Units = 'centimeters';
    fig.Color = 'w';
    fig.InvertHardcopy = 'off';
    fig.PaperPosition   = [0,10,8.45,8.45];
    fig.Position = [10,10,twocolumn/5*0.95,twocolumn/12*0.95];

    ax.XLim = [27 72];
    ax.YLim = [-0.2 0.2];

    xticks([27 72])
    yticks([])
    
    ax.XTickLabel = ["72","27"];
    ax.YTickLabel = [];
    ax.YAxis.Visible = 0;
    ax.XTickLabelRotation = 90;
    ax.YAxis.Visible = 0;

    xlabel('lightness','FontWeight','Bold','Rotation',180);
    ylabel('');

    ax.FontName = 'Arial';
    ax.Color = [.97 .97 .97];
    ax.FontSize = fontsize;
    %ticklengthcm(ax,0.0)

    ax.XColor = 'k';ax.YColor = 'k';
    
    ax.Units = 'centimeters';
    ax.Position = [0.15 0.53 3.2 1.0];

    grid off
    box off
    
    % save image
    exportgraphics(fig,fullfile(repo_basedir,'figs',['Fig9_lightness_',envType_label{envType},'.pdf']),'ContentType','vector')
    close all
end
disp('Done...')
end

%% Figure 12 - human vs. image stats

if PlotFig12
disp('Generating Figure 12...')

contrast = zeros(3,10);
coverage = zeros(3,10);
sharpness = zeros(3,10);

for envType = 1:3
    contrast(envType,:) = corrCoeff_imgStats.human.contrast(envType,:,BestId.contrast(envType,1),BestId.contrast(envType,2));
    coverage(envType,:) = corrCoeff_imgStats.human.coverage(envType,:,BestId.coverage(envType));
    sharpness(envType,:) = corrCoeff_imgStats.human.sharpness(envType,:,BestId.sharpness(envType));
end

for envType = 1:3
    for var = variableList

        if strcmp(var{1},'Pellacini_c')
            corCoeff = [coverage(envType,:)' sharpness(envType,:)' contrast(envType,:)' squeeze(corrCoeff_imgStats.human.diffuseANDspecular.(var{1})(envType,:,:))'];
        else
            corCoeff = squeeze(corrCoeff_imgStats.human.diffuseANDspecular.(var{1})(envType,:,:))';
            corCoeff_diffuseOnly = squeeze(corrCoeff_imgStats.human.diffuseOnly.(var{1})(envType,:,:))';            
            mean_corCoeff_diffuseOnly = abs(mean(corCoeff_diffuseOnly(:,orderId.(var{1}))));
        end
        mean_corCoeff = abs(mean(corCoeff(:,orderId.(var{1}))));

        if strcmp(var{1},'Pellacini_c')
            temp = std(corCoeff(:,orderId.(var{1})))/sqrt(10);
        end

        % get upper boundary (ub) and lower boundary (lb) to draw noise ceiling
        ub = squeeze(corrCoeff_imgStats.acrossParticipants_ub.(var{1})(envType,:));
        lb = squeeze(corrCoeff_imgStats.acrossParticipants_lb.(var{1})(envType,:));

        fig = figure;

        mean_gthuman_correlation = mean(corrCoeff.(var{1})(envType,:));            
        ub_gthuman_correlation = mean_gthuman_correlation+std(corrCoeff.(var{1})(envType,:))/sqrt(10);
        lb_gthuman_correlation = mean_gthuman_correlation-std(corrCoeff.(var{1})(envType,:))/sqrt(10);

        rectangle('Position',[0.1 mean(lb) 100 mean(ub)-mean(lb)],'FaceColor',[224 195 249]/255,'EdgeColor','none');hold on;

        b = bar(mean_corCoeff);hold on

        if ~strcmp(var{1},'Pellacini_c')
            for N = 1:length(mean_corCoeff_diffuseOnly)
                scatter(N,mean_corCoeff_diffuseOnly(N),30,[252 0 0]/255,'d','filled','MarkerFaceAlpha',0.5);hold on
            end
        end

        b.FaceColor = 'flat';

        if strcmp(var{1},'Pellacini_c')
            for n = 1:12
                if n < 4
                    b.CData(n,:) = [96 179 179]/255*0.8;
                else
                    b.CData(n,:) = [149 179 179]/255;
                end
            end
        else
            b.FaceColor = [149 179 179]/255;
        end

        b.EdgeColor = 'none';

        % Plotting error bars
        e = errorbar(1:length(desiredOrder.(var{1})),mean_corCoeff,std(corCoeff(:,orderId.(var{1})))/sqrt(10),'LineStyle','none','CapSize',0);
        e.Color = [0 0 0];e.LineWidth = .5;
        
        axis normal

        % Make figures look nicer
        ax = gca;
        fig.Units           = 'centimeters';
        fig.Position = [10,10,twocolumn/3*0.95,twocolumn/4*0.95];

        fig.Color           = 'w';
        fig.InvertHardcopy  = 'off';

        ax.XLim = [0 length(desiredOrder.(var{1}))+1];
        ax.YLim = [0 1.05];

        xticks(1:length(desiredOrder.(var{1})))
        yticks([0 0.5 1])
        ax.XTickLabel = stats_label.(var{1});
        ax.YTickLabel = {'0.0','0.5','1.0'};
        ax.XTickLabelRotation = 45;

        ax.FontName = 'Arial';
        ax.Color = 'w';
        ax.FontSize = fontsize;
        ax.XColor = 'k';ax.YColor = 'k';

        ax.LineWidth = 0.5;
        ax.Units = 'centimeters';
        ax.Position = [0.45 0.9 5.4 3.2];

        ticklengthcm(ax,0.0)
        grid off
        box off

        exportgraphics(fig,fullfile(repo_basedir,'figs',['Fig12_',(var{1}),'_',envType_label{envType},'.png']),'Resolution','600')
    end
end

close all

%% combine panels in Figure 12 and show
clear imgs
imgsize_w = 1375;
imgsize_h = round(imgsize_w*(953/1375));

gap = round(imgsize_h*0.1);
out2 = [];

% 1: lightness,  2:Hue, 3:chroma, 4:gloss

for envType = 1:3
    out = [];
    %for p = [2 1 3 4]
    for var = {'lightness','chroma','Pellacini_c'}
        out = [out,ones(imgsize_h,gap,3)*255,imresize(imread(fullfile(repo_basedir,'figs',['Fig12_',(var{1}),'_',envType_label{envType},'.png'])),[imgsize_h imgsize_w])];
    end
    out2 = [out2;ones(round(gap*2),size(out,2),3)*255;out];
end

figure;imshow(out2)
imwrite(out2, fullfile(repo_basedir,'figs','Fig12_combined.png'))
close all
disp('Done.')
end

%% Custom function to sort human response data (as the presentation order was randomized in the experiment)
function [responses_sorted, groundtruth_sorted] = sortData(RESPONSES,GROUNDTRUTH,order)
    % function to order all field based on order
    responses_sorted.specularity = RESPONSES.specularity(order);
    responses_sorted.Pellacini_c = RESPONSES.Pellacini_c(order);
    responses_sorted.lightness = RESPONSES.lightness(order);
    responses_sorted.chroma = RESPONSES.chroma(order);
    responses_sorted.hue = RESPONSES.hue(order);
    responses_sorted.saturation = responses_sorted.chroma./responses_sorted.lightness;

    groundtruth_sorted.specularity = GROUNDTRUTH.specularity(order);
    groundtruth_sorted.Pellacini_c = GROUNDTRUTH.c_Pellacini(order);
    groundtruth_sorted.lightness = GROUNDTRUTH.lightness(order);
    groundtruth_sorted.chroma = GROUNDTRUTH.chroma(order);
    groundtruth_sorted.hue = GROUNDTRUTH.hue(order);
    groundtruth_sorted.saturation = groundtruth_sorted.chroma./groundtruth_sorted.lightness;
    groundtruth_sorted.imageN = GROUNDTRUTH.imageN(order);
end

%% Custom function to sort human response data acorss sessions
function humanResponse = UniteSessions(RESPONSES,sessionN)
    humanResponse.specularity = [];
    humanResponse.Pellacini_c = [];
    humanResponse.lightness = [];
    humanResponse.chroma = [];
    humanResponse.hue = [];
    humanResponse.saturation = [];
    
    for session = 1:sessionN
        humanResponse.specularity = [humanResponse.specularity,RESPONSES(session).specularity'];
        humanResponse.Pellacini_c = [humanResponse.Pellacini_c,RESPONSES(session).Pellacini_c'];
        humanResponse.lightness = [humanResponse.lightness,RESPONSES(session).lightness'];
        humanResponse.chroma = [humanResponse.chroma,RESPONSES(session).chroma'];
        humanResponse.hue = [humanResponse.hue,RESPONSES(session).hue'];
        humanResponse.saturation = [humanResponse.saturation,RESPONSES(session).saturation'];
    end
end

%% Custom function to sort response data across all human observers
function humanResponse_allobservers = UniteObservers(humanResponse)
    Outputs = [size(humanResponse(1).hue),length(humanResponse)];
    humanResponse_allobservers.specularity = [];
    humanResponse_allobservers.Pellacini_c = [];
    humanResponse_allobservers.lightness = [];
    humanResponse_allobservers.chroma = [];
    humanResponse_allobservers.hue = [];
    humanResponse_allobservers.saturation = [];
    
    for observer = 1:length(humanResponse)
        humanResponse_allobservers.specularity = [humanResponse_allobservers.specularity,humanResponse(observer).specularity];
        humanResponse_allobservers.Pellacini_c = [humanResponse_allobservers.Pellacini_c,humanResponse(observer).Pellacini_c];
        humanResponse_allobservers.lightness = [humanResponse_allobservers.lightness,humanResponse(observer).lightness];
        humanResponse_allobservers.chroma = [humanResponse_allobservers.chroma,humanResponse(observer).chroma];
        humanResponse_allobservers.hue = [humanResponse_allobservers.hue,humanResponse(observer).hue];
        humanResponse_allobservers.saturation = [humanResponse_allobservers.saturation,humanResponse(observer).saturation];
    end
    
    % take average across sessions
    humanResponse_allobservers.specularity = squeeze(mean(reshape(humanResponse_allobservers.specularity,Outputs),2));
    humanResponse_allobservers.Pellacini_c = squeeze(mean(reshape(humanResponse_allobservers.Pellacini_c,Outputs),2));
    humanResponse_allobservers.lightness = squeeze(mean(reshape(humanResponse_allobservers.lightness,Outputs),2));
    humanResponse_allobservers.chroma = squeeze(mean(reshape(humanResponse_allobservers.chroma,Outputs),2));
    humanResponse_allobservers.hue = mod(rad2deg(squeeze(circ_mean(reshape(deg2rad(humanResponse_allobservers.hue),Outputs),[],2)+2*pi)),360);
    humanResponse_allobservers.saturation = squeeze(mean(reshape(humanResponse_allobservers.saturation,Outputs),2));
end

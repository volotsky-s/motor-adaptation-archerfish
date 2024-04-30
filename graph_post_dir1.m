
data_mat = readtable('post_epoch_subject_dir1.csv');

show_ep = [1 2  4 5  ];
legend_ep = {'baseline','pert start','pert middle','pert end', ...
    'washout start','washout end'};

y_S1 = zeros(size(data_mat,1),length(show_ep));
figure('Position',[50 100 1200 500]); hold on
for e = 1:length(show_ep)
    epoch_str = ['bE_' num2str(show_ep(e)) '_'];
    int_str = ['bES_' num2str(show_ep(e)) '_1_'];
    
    y_S1(:,e) = data_mat.b0+data_mat.(epoch_str)+data_mat.bS_1_ + ...
                data_mat.(int_str);
    histogram(y_S1(:,e),-5:0.1:9)
end
xlabel('Posterior Distribution Error')
legend(legend_ep(show_ep))
title('Posterior Distribution: fish 1')

y_S2 = zeros(size(data_mat,1),length(show_ep));
figure('Position',[50 100 1200 500]); hold on
for e = 1:length(show_ep)
    epoch_str = ['bE_' num2str(show_ep(e)) '_'];
    int_str = ['bES_' num2str(show_ep(e)) '_2_'];

    y_S2(:,e) = data_mat.b0+data_mat.(epoch_str)+data_mat.bS_2_ + ...
                data_mat.(int_str);
    histogram(y_S2(:,e),-5:0.1:9)
end
xlabel('Posterior Distribution Error')
legend(legend_ep(show_ep))
title('Posterior Distribution: fish 2')

y_S3 = zeros(size(data_mat,1),length(show_ep));
figure('Position',[50 100 1200 500]); hold on
for e = 1:length(show_ep)
    epoch_str = ['bE_' num2str(show_ep(e)) '_'];
    int_str = ['bES_' num2str(show_ep(e)) '_3_'];

    y_S3(:,e) = data_mat.b0+data_mat.(epoch_str)+data_mat.bS_3_ + ...
                data_mat.(int_str);
    histogram(y_S3(:,e),-5:0.1:9)
end
xlabel('Posterior Distribution Error')
legend(legend_ep(show_ep))
title('Posterior Distribution: fish 3')

y_S123 = zeros(size(data_mat,1),length(show_ep));
figure('Position',[50 100 1200 500]); hold on
for e = 1:length(show_ep)
    epoch_str = ['bE_' num2str(show_ep(e)) '_'];

    y_S123(:,e) = data_mat.b0+data_mat.(epoch_str);
    histogram(y_S123(:,e),-5:0.1:9)
end
xlabel('Posterior Distribution Error')
legend(legend_ep(show_ep))
title('Posterior Distribution: fish 123')

% effect size for posterior values
es12 = (y_S123(:,2)-y_S123(:,1))./data_mat.ySD;
es24 = (y_S123(:,3)-y_S123(:,2))./data_mat.ySD;
es45 = (y_S123(:,4)-y_S123(:,3))./data_mat.ySD;


%%
% hdi for each fish and for all fish 4 epochs
hdi_y_123 = zeros(4,2); 
for e = 1:4
hdi_y_123(e,:) = find_hdi(y_S123(:,e),0.95);
end

hdi_y_1 = zeros(4,2); 
for e = 1:4
hdi_y_1(e,:) = find_hdi(y_S1(:,e),0.95);
end

hdi_y_2 = zeros(4,2); 
for e = 1:4
hdi_y_2(e,:) = find_hdi(y_S2(:,e),0.95);
end

hdi_y_3 = zeros(4,2); 
for e = 1:4
hdi_y_3(e,:) = find_hdi(y_S3(:,e),0.95);
end

% hdi for difference between epochs for all fish
hdi_diff_dir1(1,:) = find_hdi(y_S123(:,2)-y_S123(:,1),0.95);
hdi_diff_dir1(2,:) = find_hdi(y_S123(:,3)-y_S123(:,2),0.95);
hdi_diff_dir1(3,:) = find_hdi(y_S123(:,4)-y_S123(:,3),0.95);

% hdi for comparison between epochs effect size for all fish
hdi_es_dir1(1,:) = find_hdi(es12,0.95);
hdi_es_dir1(2,:) = find_hdi(es24,0.95);
hdi_es_dir1(3,:) = find_hdi(es45,0.95);


%% side 1 - fish 1-3
clear l
figure('position',[50 50 900 600]); hold on

patch([0 2 2 0],[-10 -10 -7 -7],'g','facealpha', 0.3,'EdgeColor','w')
patch([2 6 6 2],[-10 -10 -7 -7],'r','facealpha', 0.3,'EdgeColor','w')
patch([6 10 10 6],[-10 -10 -7 -7],'g','facealpha', 0.3,'EdgeColor','w')

text(0.1, -7.6,'Baseline','fontsize',24)
text(2.8, -7.6,'Perturbation','fontsize',24)
text(6.1, -7.6,'Washout','fontsize',24)
text(2.5, -9,'start','fontsize',24)
text(4.9, -9,'end','fontsize',24)

line([0 10],[0 0], 'color','k')
line([2 2],[-10 10], 'color','k')
line([6 6],[-10 10], 'color','k')

for e = 1:4
     bar(e*2-1,mean(y_S123(:,e)),1, 'facecolor','w')

    l(4) = line([e*2-0.6 e*2-0.6],[hdi_y_123(e,1) hdi_y_123(e,2)],'linewidth',5,'color','k');

    l(1) = line([e*2-1.4 e*2-1.4],[hdi_y_1(e,1) hdi_y_1(e,2)],'linewidth',3,'color',[0 0.4470 0.7410]);
    l(2) = line([e*2-1.2 e*2-1.2],[hdi_y_2(e,1) hdi_y_2(e,2)],'linewidth',3,'color',[0.8500 0.3250 0.0980]);
    l(3) = line([e*2-1 e*2-1],[hdi_y_3(e,1) hdi_y_3(e,2)],'linewidth',3,'color',[0.9290 0.6940 0.1250]);

end
line([0 8],[-10 -10],'color','k')
line([0 0],[-10 8],'color','k')

xlim([0 8])
ylim([-10 8])

yticks([-5 0 5])
row1 = {'baseline','perturbation','perturbation', 'washout'};
row2 = {' ','      start', '       end', ' '};
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
ax = gca(); 

xticklabels([])

legend(l,{'Fish 1', 'Fish 2', 'Fish 3', 'Fish 1-3' }, ...
    'box','off','position', [0.8 0.63 0.1 0.3])

set(gca,'position', [0.22 0.1 0.65 0.75],'fontname','helvetica','fontsize',25, ...
    'Box','off','TickDir','out','TickLength',[.001 .001])

ylabel({'Error [mm] ', 'Mean \pm 95% HDI'},'fontname','helvetica','fontsize',36)

% exportgraphics(gcf,'summary_dir1.tiff', 'Resolution',300)


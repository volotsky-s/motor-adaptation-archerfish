
data_mat = readtable('post_epoch_subject_reverse.csv');

show_ep = [1 2  4  5 7  8];
legend_ep = {'baseline','pert start','pert middle','pert end', ...
    'pert II start','pert II middle','pert II end', ...
    'washout start','washout end'};

y_S3 = zeros(size(data_mat,1),length(show_ep));
figure('Position',[50 100 1200 500]); hold on
for e = 1:length(show_ep)
    epoch_str = ['bE_' num2str(show_ep(e)) '_'];
    int_str = ['bES_' num2str(show_ep(e)) '_1_'];
    
    y_S3(:,e) = data_mat.b0+data_mat.(epoch_str)+data_mat.bS_1_ + ...
                data_mat.(int_str);
    histogram(y_S3(:,e),-25:0.1:20)
end
xlabel('Posterior Distribution Error')
legend(legend_ep(show_ep))
title('Posterior Distribution: fish 3')

y_S6 = zeros(size(data_mat,1),length(show_ep));
figure('Position',[50 100 1200 500]); hold on
for e = 1:length(show_ep)
    epoch_str = ['bE_' num2str(show_ep(e)) '_'];
    int_str = ['bES_' num2str(show_ep(e)) '_2_'];

    y_S6(:,e) = data_mat.b0+data_mat.(epoch_str)+data_mat.bS_2_ + ...
                data_mat.(int_str);
    histogram(y_S6(:,e),-25:0.1:25)
end
xlabel('Posterior Distribution Error')
legend(legend_ep(show_ep))
title('Posterior Distribution: fish 6')

y_S7 = zeros(size(data_mat,1),length(show_ep));
figure('Position',[50 100 1200 500]); hold on
for e = 1:length(show_ep)
    epoch_str = ['bE_' num2str(show_ep(e)) '_'];
    int_str = ['bES_' num2str(show_ep(e)) '_3_'];

    y_S7(:,e) = data_mat.b0+data_mat.(epoch_str)+data_mat.bS_3_ + ...
                data_mat.(int_str);
    histogram(y_S7(:,e),-25:0.1:25)
end
xlabel('Posterior Distribution Error')
legend(legend_ep(show_ep))
title('Posterior Distribution: fish 7')

y_S367 = zeros(size(data_mat,1),length(show_ep));
figure('Position',[50 100 1200 500]); hold on
for e = 1:length(show_ep)
    epoch_str = ['bE_' num2str(show_ep(e)) '_'];

    y_S367(:,e) = data_mat.b0+data_mat.(epoch_str)  ;
    histogram(y_S367(:,e),-25:0.1:20)
end
xlabel('Posterior Distribution Error')
legend(legend_ep(show_ep))
title('Posterior Distribution: fish 367')

% effect size for posterior values
es12 = (y_S367(:,2)-y_S367(:,1))./data_mat.ySD;
es24 = (y_S367(:,3)-y_S367(:,2))./data_mat.ySD;
es45 = (y_S367(:,4)-y_S367(:,3))./data_mat.ySD;
es57 = (y_S367(:,5)-y_S367(:,4))./data_mat.ySD;
es78 = (y_S367(:,6)-y_S367(:,5))./data_mat.ySD;


%%


hdi_y_367 = zeros(4,2); 
for e = 1:6
hdi_y_367(e,:) = find_hdi(y_S367(:,e),0.95);
end

hdi_y_3 = zeros(4,2); 
for e = 1:6
hdi_y_3(e,:) = find_hdi(y_S3(:,e),0.95);
end

hdi_y_6 = zeros(4,2); 
for e = 1:6
hdi_y_6(e,:) = find_hdi(y_S6(:,e),0.95);
end

hdi_y_7 = zeros(4,2); 
for e = 1:6
hdi_y_7(e,:) = find_hdi(y_S7(:,e),0.95);
end

% hdi for difference between epochs for all fish
hdi_diff(1,:) = find_hdi(y_S367(:,2)-y_S367(:,1),0.95);
hdi_diff(2,:) = find_hdi(y_S367(:,3)-y_S367(:,2),0.95);
hdi_diff(3,:) = find_hdi(y_S367(:,4)-y_S367(:,3),0.95);
hdi_diff(4,:) = find_hdi(y_S367(:,2)-y_S367(:,1),0.95);
hdi_diff(5,:) = find_hdi(y_S367(:,3)-y_S367(:,2),0.95);

% hdi for comparison between epochs effect size for all fish
hdi_es(1,:) = find_hdi(es12,0.95);
hdi_es(2,:) = find_hdi(es24,0.95);
hdi_es(3,:) = find_hdi(es45,0.95);
hdi_es(4,:) = find_hdi(es57,0.95);
hdi_es(5,:) = find_hdi(es78,0.95);

%%
clear l

figure('position',[20 50 3300 500]); hold on

c1 = [0.8 0.1 0.2];
c2 = [0.1 0.7 0.4]; 
c3 = [0.3 0.1 0.5]; 


patch([0 2.2 2.2 0],[-25 -25 -20 -20],'g','facealpha', 0.3,'EdgeColor','w')
patch([2.2 5.5 5.5 2.2],[-25 -25 -20 -20],'r','facealpha', 0.3,'EdgeColor','w')
patch([5.5 10 10 5.5],[-25 -25 -20 -20],[1 0 0],'facealpha', 0.6,'EdgeColor','w')
patch([10 12 12 10],[-25 -25 -20 -20],'g','facealpha', 0.3,'EdgeColor','w')

text(0.3, -21,'Baseline','fontsize',26)
text(2.6, -21,'Perturbation','fontsize',26)
text(2.6, -23.5,'start','fontsize',26)
text(4.2, -23.5,'end','fontsize',26)
text(5.9, -21,'Reverse Perturbation','fontsize',26)
text(6.5, -23.5,'start','fontsize',26)
text(8.7, -23.5,'end','fontsize',26)
text(10.2, -21,'Washout','fontsize',26)
text(10.5, -23.5,'start','fontsize',26)


line([0 25],[0 0], 'color','k')
line([2.2 2.2],[-25 10], 'color','k')
line([5.5 5.5],[-25 10], 'color','k')
line([10 10],[-25 10], 'color','k')

bar_place = [1 3 4.5 7 9 11];
for e = 1:6
     bar(bar_place(e),mean(y_S367(:,e)),1, 'facecolor','w')

    l(4) = line([bar_place(e)+0.4 bar_place(e)+0.4],[hdi_y_367(e,1) hdi_y_367(e,2)],'linewidth',5,'color','k');

    l(1) = line([bar_place(e)-0.4 bar_place(e)-0.4],[hdi_y_3(e,1) hdi_y_3(e,2)],'linewidth',3,'color',[0.6350 0.0780 0.1840]);
    l(2) = line([bar_place(e)-0.2 bar_place(e)-0.2],[hdi_y_6(e,1) hdi_y_6(e,2)],'linewidth',3,'color',c3);
    l(3) = line([bar_place(e) bar_place(e)],[hdi_y_7(e,1) hdi_y_7(e,2)],'linewidth',3,'color',c2);

end
line([0 12],[-25 -25],'color','k')
line([0 0],[-25 10],'color','k')

xlim([0 12])
ylim([-25 10])

yticks([-20 -10 0 10])
xticks(1:2:11)
row1 = {'baseline','perturbation','perturbation', '    reverse', '    reverse','washout'};
row2 = {'       ' , '      start' , '        end','perturbation', 'perturbation','     '};
row3 = {'       ','          ', '           ','      start', '       end','     '};
labelArray = [row1; row2; row3]; 
tickLabels = strtrim(sprintf('%s\\newline%s\\newline%s\n', labelArray{:}));
ax = gca(); 

xticklabels([])

legend(l,{'Fish 3', 'Fish 6', 'Fish 7','Fish 3,6,7' }, ...
    'box','off', 'position', [0.88 0.33 0.1 0.3]) 

set(gca,'position', [0.15 0.2 0.78 0.7],'fontname','helvetica','fontsize',25, ...
    'Box','off','TickDir','out','TickLength',[.001 .001])

ylabel({'Error [mm] ', 'Mean \pm 95% HDI'},'fontname','helvetica','fontsize',36)

% exportgraphics(gcf,'summary_reverse.tiff', 'Resolution',300)


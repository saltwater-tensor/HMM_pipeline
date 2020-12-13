function [] = within_state(lifetimes,states_to_test,hmmstr)

for st1 = 1:1:3
    for st2 = 1:1:3
        
        r = cell2mat(lifetimes(:,states_to_test(st1))');
%         r(find(r < 100)) = [];
        
        s = cell2mat(lifetimes(:,states_to_test(st2))');
%         s(find(s <100)) = [];
%         
%         r(find(r > 1000)) = [];
%         s(find(s > 1000)) = [];
        
        x = (r');
        y = (s');
        [h_row_lesser_than_col(st1,st2),p_off_lesser_than_on(st1,st2)] = ...
            ttest2(x,y,'VarType','unequal','Tail','left') %row < col
        [h_row_greater_than_col(st1,st2),p_off_greater_than_on(st1,st2)] = ...
            ttest2(x,y,'VarType','unequal','Tail','right') %row > col
    end
end

% multiple corrections
p_off_greater_than_on = reshape(mafdr(p_off_greater_than_on(:)),[3,3]);
fig_handle = figure(200);
h = heatmap(p_off_greater_than_on);
h.XLabel = 'State'; %columns
h.YLabel = 'State'; %rows
h.Title = 'p positive ROW > COLUMN';
% h.XDisplayLabels = {num2str(states_to_test(1)),num2str(states_to_test(2)),num2str(states_to_test(3))};
% h.YDisplayLabels = {num2str(states_to_test(1)),num2str(states_to_test(2)),num2str(states_to_test(3))};

h.XDisplayLabels = {'High Pwr','Comms','Local'};
h.YDisplayLabels = {'High Pwr','Comms','Local'};


saveas(fig_handle,['p_positive_HMM_' hmmstr]);
saveas(fig_handle,['p_positive_HMM_' hmmstr],'png');
close(fig_handle)
% multiple corrections
p_off_lesser_than_on = reshape(mafdr(p_off_lesser_than_on(:)),[3,3]);
fig_handle = figure(200);
h = heatmap(p_off_lesser_than_on);
h.XLabel = 'State'; %columns
h.YLabel = 'State'; %rows
h.Title = 'p negative ROW < COLUMN';
% h.XDisplayLabels = {num2str(states_to_test(1)),num2str(states_to_test(2)),num2str(states_to_test(3))};
% h.YDisplayLabels = {num2str(states_to_test(1)),num2str(states_to_test(2)),num2str(states_to_test(3))};

h.XDisplayLabels = {'High Pwr','Comms','Local'};
h.YDisplayLabels = {'High Pwr','Comms','Local'};

saveas(fig_handle,['p_negative_HMM_' hmmstr]);
saveas(fig_handle,['p_negative_HMM_' hmmstr],'png');
close(fig_handle)

end
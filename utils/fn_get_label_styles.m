function [labels, colors, line_styles] = fn_get_label_styles(var,n_levels)

switch var
    case 'sbj_n'
        labels = {'PFC03','PFC04','PFC05','PFC01'};
    case 'reward_cur'
        labels = {'Current Reward'};
    case 'reward_prv'
        labels = {'Previous Reward'};
    case 'effortS_cur'
        labels = {'Current Effort'};
    case 'effortS_prv'
        labels = {'Previous Effort'};
    case 'pAccept_cur'
        labels = {'Current p(Accept)'};
    otherwise
        labels = {strrep(var,'_',' ')};
end

% Colors from lowest to highest
if contains(var,'reward')
    full_colors = [201,148,199;...
                      223,101,176;...
                      231,41,138;...
                      206,18,86;...
                      145,0,63]./256;   % shades of red
elseif contains(var,'effort')
    full_colors = [166,189,219;...
                      116,169,207;...
                      54,144,192;...
                      5,112,176;...
                      3,78,123]./256; % shades of blue
elseif strcmp(var,'sbj_n')
    full_colors = [27, 158, 119;         % teal
                   217, 95, 2;           % burnt orange
                   117, 112, 179;        % purple
                   231, 41, 138]./256;   % magenta
elseif contains(var,'pAccept')
    full_colors = [0 0 0;
                   0 0 0;
                   0 0 0;
                   0.5 0.5 0.5;
                   0.5 0.5 0.5];
else
    full_colors = zeros([5, 3]); % black (power variables)
end

if n_levels==1
    if contains(var,'decision'); color_ix = 1; else, color_ix = 5; end
    line_styles = {'-'};
elseif n_levels==2
    color_ix = [3 5];
    labels = {[labels{1} ' Low'], [labels{1} ' High']};
    line_styles = {'-','--'};
elseif n_levels==3
    color_ix = [1 3 5];
    labels = {[labels{1} ' Low'], [labels{1} ' Middle'], [labels{1} ' High']};
    line_styles = {'-','--',':'};
elseif n_levels==4
    if ~strcmp(var,'sbj_n'); error('use 4 levels for subjects!'); end
    color_ix = 1:4;
    line_styles = {'-','-','-','-'};
elseif n_levels==5
    color_ix = 1:n_levels;
    line_styles = {'-','-','-','-','-'};
    tmp_lab = labels; labels = cell([n_levels 1]);
    for n = 1:n_levels
        labels{n} = [tmp_lab{1} ' Q' num2str(n)];
    end
end
colors = full_colors(color_ix,:);

end
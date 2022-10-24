function Save_all_figures(path)
if ~nargin
    path='.\';
end
% No parameters for this function. The range of handles of figures is
% between 1 to 100. If the saved figure is out of this range, this range
% can be modified.

% count_fig=0;
% if ishandle(1)
%     saveas(1,'1','fig');
%     close 1
%     count_fig=count_fig+1;
% end
% 
% while (i~=0)
% fig=gcf;
% if fig==1
%     close 1
%     break
% else
% saveas(fig,num2str(fig),'fig');
% close gcf
% count_fig=count_fig+1;
% end
% end
% disp([ num2str(count_fig) ' figures have been saved!'])

% No parameters for this function. And The range of handles of figures is
% between 1 to 100. If the saved figure is out of this range, this range
% can be modified.
h = get(0,'children');
if isempty(h)
    disp('There is no current figure that can be saved! ')
    return
else
    for i=1:length(h)
        saveas(h(i).Number,[ path num2str(h(i).Number)], 'fig');
    end
    disp([ num2str(length(h)) ' figures have been saved!'])
end
    
end
function output_txt = myDataTipFunction(obj,event_obj)
% Display data cursor position in a data tip
% obj          Currently not used
% event_obj    Handle to event object
% output_txt   Data tip text, returned as a character vector or a cell array of character vectors

pos = event_obj.Position;


%********* Define the content of the data tip here *********%
% search for the correct data point
    xIdx = find(event_obj.Target.XData == pos(1));
    yIdx = find(event_obj.Target.YData == pos(2));
    % zIdx = find(event_obj.Target.ZData == pos(3));
    % select single data point, and find value
    % idx = intersect(intersect(xIdx,yIdx),zIdx);
    idx = intersect(xIdx,yIdx);
    value = event_obj.Target.CData(idx(1)); 
    % add to the data cursor text


% Display the x and y values:
output_txt = {['X',formatValue(pos(1),event_obj)],...
    ['Y',formatValue(pos(2),event_obj)],...
    ['Val',formatValue(value,event_obj)]};

%***********************************************************%


% If there is a z value, display it:
if length(pos) > 2
    output_txt{end+1} = ['Z',formatValue(pos(3),event_obj)];
end

%***********************************************************%

function formattedValue = formatValue(value,event_obj)
% If you do not want TeX formatting in the data tip, uncomment the line below.
% event_obj.Interpreter = 'none';
if strcmpi(event_obj.Interpreter,'tex')
    valueFormat = ' \color[rgb]{0 0.6 1}\bf';
    removeValueFormat = '\color[rgb]{.25 .25 .25}\rm';
else
    valueFormat = ': ';
    removeValueFormat = '';
end
formattedValue = [valueFormat num2str(value,4) removeValueFormat];

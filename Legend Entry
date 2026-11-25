function [hs, names] = legendEntries(handles, labels)
   hs = gobjects(0);  names = {};
   for k = 1:numel(handles)
       if ~isempty(handles{k}) && isgraphics(handles{k})
           hs(end+1) = handles{k}(1); %#ok<AGROW>
           names{end+1} = labels{k};  %#ok<AGROW>
       end
   end
end

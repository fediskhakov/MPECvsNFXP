function table = tabulate(x)
%TABULATE Frequency table.
%   TABLE = TABULATE(X) takes a vector X and returns a matrix, TABLE.
%   The first column of TABLE contains the unique values of X.  The
%   second is the number of instances of each value.  The last column
%   contains the percentage of each value.  If the elements of X are
%   non-negative integers, then the output includes 0 counts for any
%   integers that are between 1 and max(X) but do not appear in X.
%
%   TABLE = TABULATE(X), where X is a categorical variable, character
%   array, or a cell array of strings, returns TABLE as a cell array.  The
%   first column contains the unique string values in X, and the other two
%   columns are as above.
%
%   TABULATE with no output arguments returns a formatted table
%   in the command window.
%
%   See also PARETO.
   
%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.4 $  $Date: 2011/02/09 19:35:39 $

isnum = isnumeric(x);
if isnum && ~isfloat(x)
    % use of hist() below requires float
    x = double(x);
end
if isnum
   if min(size(x)) > 1,
      error(message('stats:tabulate:InvalidData'));
   end

   y = x(~isnan(x));
else
   y = x;
end

if ~isnum || any(y ~= round(y)) || any(y < 1);
   docell = true;
   [y,yn,yl] = grp2idx(y);
   maxlevels = length(yn);
else
   docell = false;
   maxlevels = max(y);
   %yn = cellstr(num2str((1:maxlevels)'));
end

[counts values] = hist(y,(1:maxlevels));

total = sum(counts);
percents = 100*counts./total;

if nargout == 0
   if docell
      width = max(cellfun('length',yn));
      width = max(5, min(50, width));
   else
      width = 5;
   end
   
   % Create format strings similar to:   '  %5s    %5d    %6.2f%%\n'
   fmt1 = sprintf('  %%%ds    %%5d    %%6.2f%%%%\n',width);
   fmt2 = sprintf('  %%%ds    %%5s   %%6s\n',width);
   fprintf(1,fmt2,'Value','Count','Percent');
   if docell
      for j=1:maxlevels
         fprintf(1,fmt1,yn{j},counts(j),percents(j));
      end
   else
      fprintf(1,'  %5d    %5d    %6.2f%%\n',[values' counts' percents']');
   end
else
   if ~docell
      table = [values' counts' percents'];
   elseif isnum
      table = [yl(:) counts' percents'];
   else
      table = [yn num2cell([counts' percents'])];
   end
end

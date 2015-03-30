function [cstruct, c] = splitcontours(c)
% Split the contour lines in structure array
%
% [cstruct, c] = SPLITCONTOURS(c)
%
% split the contour array returned by
% Matlab built-in function contourc() as a structure array, each
% contains a single contour lines
%
% by Bruno Luong
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/293174


% First pass to count the number of contour lines
l = 0;
idxlabel = 1;
while idxlabel <= size(c,2)
    n = c(2,idxlabel);
    l = l+1;
    idxlabel = idxlabel+n+1;
end
% Allocation
idxremove = zeros(1,l);
cstruct = struct('x', cell(1,l), 'y', cell(1,l));
% Second pass to fill (x,y)
l = 0;
idxlabel = 1;
while idxlabel <= size(c,2)
    l = l+1;
    idxremove(l) = idxlabel;
    n = c(2,idxlabel);
    cstruct(l).x = c(1,idxlabel+(1:n));
    cstruct(l).y = c(2,idxlabel+(1:n));
    idxlabel = idxlabel+n+1;
end
c(:,idxremove) = [];

end % splitcontours
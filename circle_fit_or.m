function [xo yo R variation angles] = circle_fit_or(x,y)
% A function to find the best circle fit (radius and center location) to
% given x,y pairs
% 
% Val Schmidt
% Center for Coastal and Ocean Mapping
% University of New Hampshire
% 2012
% 
% Arguments:
% x:         x coordinates
% y:         y coordinates
%
% Output:
% xo:        circle x coordinate center
% yo:        circle y coordinate center
% R:         circle radius
x = x(:);
y = y(:);

% copies for later
x1 = x;
y1 = y;

% Fitting a circle to the data - least squares style. 
%Here we start with
% (x-xo).^2 + (y-yo).^2 = R.^2
% Rewrite:
% x.^2 -2 x xo + xo^2 + y.^2 -2 y yo + yo.^2 = R.^2
% Put in matrix form:
% [-2x -2y 1 ] [xo yo -R^2+xo^2+yo^2]' = -(x.^2 + y.^2)
% Solve in least squares way...
A = [-2*x -2*y ones(length(x),1)];
x = A\-(x.^2+y.^2);
xo=x(1);
yo=x(2);
R = sqrt(  xo.^2 + yo.^2  - x(3));

% extract orientation of the circle
if size(x1,1)>size(x1,2)
    pts = [x1 y1];
else
    pts = [x1' y1'];
end
for i1=1:size(pts,1)
   pts(i1,3) = atan2(pts(i1,2)-yo,pts(i1,1)-xo);
end

% In some cases, angles can vary in non-monotonous way because of cyclicity
% of angles on the circle. In this case, add 2*pi to points so that we get
% monotonous variation.

% if angles are monotonous, increment should have identical signs
incr = diff(pts(:,3));

if not(sum(sign(incr))==length(incr) || sum(sign(incr))==-length(incr))
    for i1=1:size(pts,1)
       if pts(i1,3)<0
           pts(i1,3) = pts(i1,3) + 2*pi;
       end
    end
end

% finally, extract the orientation
tempa = sort(pts(:,3),'ascend');
tempd = sort(pts(:,3),'descend');

variation = NaN;
if isequal(tempa,pts(:,3))
    % clockwise
    variation = 1;
else
    if isequal(tempd,pts(:,3))
        % counter clockwise
        variation = -1;
    end
end

angles = pts(:,3);

end



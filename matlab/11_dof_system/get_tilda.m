function x_tilda = get_tilda(x)
% given y: x_tilda*y = cross(x,y)
% [Stoneking, page 3 bottom]
%
    x_tilda = [0, -x(3), x(2); ... 
               x(3), 0, -x(1); ...
               -x(2), x(1), 0];

end
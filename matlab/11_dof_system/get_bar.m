function x_bar = get_bar(x)
% given y: x_bar*y = cross(x, cross(x,y))
% [Stoneking, page 3 bottom]
%
   x_bar = [-x(2)^2 - x(3)^2, x(1)*x(2),        x(1)*x(3);...
             x(2)*x(1),      -x(3)^2-x(1)^2,    x(2)*x(3);...
             x(3)*x(1),       x(3)*x(2),        -x(1)^2-x(2)^2];
end
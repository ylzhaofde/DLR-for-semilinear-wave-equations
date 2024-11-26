function [A0,B0] = initials(x,y)
    [x,y] = meshgrid(x,y);
    A0 = 2*sin(3*pi*x).*sin(3*pi*y);
    B0 = -sin(3*pi*x).*sin(3*pi*y);
    
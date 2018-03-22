function h = circle(x,y,r)
hold on;
for i=1:length(x)
    th = 0:pi/50:2*pi;
    xunit = r(i) * cos(th) + x(i);
    yunit = r(i) * sin(th) + y(i);
    h = fill(xunit, yunit,'k');
end
end
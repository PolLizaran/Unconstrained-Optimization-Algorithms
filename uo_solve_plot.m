% Start [uo_solve_plot] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F.-Javier Heredia https://gnom.upc.edu/heredia
function [] = uo_solve_plot(f, xk, gk, xylim, iplot,fs)
%
% Plots a funtion f(x) over R^2 together with the
% sequence {x^k} generated by an optimization algorithm:
%
% f     : name of objective function.
% xk    : matrix with rows storing the iterated points {x^k}.
% gk    : matrix with rows storing the gradient at the iterated points {g(x^k)}.
% xylim : row vector [xmin, xmax, ymin,ymax] defining the limits of the
%       : ploting area. If [0,0,0,0] the ploting ranges are calculated
%       : to accomodate the whole sequence {x^k}.
% iplot : If =1, level curves, {x^k}, {g^k}.
%       : If =2, level curves, {x^k}, {g^k} and surface.
% fs    : if fs>0, sets the font size of the legends to "fs" points.
%       : Any value <=0 is ignored.
%
cla;
hold off;
grid on;
box on;
if fs > 0 set(gca,'Fontsize',fs); end
hold on;
% Plot area
xmin=xylim(1);
xmax=xylim(2);
ymin=xylim(3);
ymax=xylim(4);
a=.3;
max_xk1   = max(xk(1,:));
min_xk1   = min(xk(1,:));
w_xk1 = max_xk1-min_xk1;
max_xk1   = max_xk1 + a*w_xk1;
min_xk1   = min_xk1 - a*w_xk1;
w_xk1 = max_xk1-min_xk1;

max_xk2   = max(xk(2,:));
min_xk2   = min(xk(2,:));
w_xk2 = max_xk2-min_xk2;
max_xk2   = max_xk2 + a*w_xk2;
min_xk2   = min_xk2 - a*w_xk2;
w_xk2 = max_xk2-min_xk2;

difwidth = w_xk1-w_xk2;
if difwidth > 0
    max_xk2 = max_xk2 + difwidth/2;
    min_xk2 = min_xk2 - difwidth/2;
elseif difwidth < 0
    max_xk1 = max_xk1 - difwidth/2;
    min_xk1 = min_xk1 + difwidth/2;
end

if xmax~=xmin
    if xmax < max_xk1 xmax = max_xk1; end
    if xmin > min_xk1 xmin = min_xk1; end
else
    xmax = max_xk1;
    xmin = min_xk1;
end
if ymax~=ymin
    if ymax < max_xk2 ymax = max_xk2; end
    if ymin > min_xk2 ymin = min_xk2; end
else
    ymax = max_xk2;
    ymin = min_xk2;
end
%
set(gcf,'Color','w');  
set(findobj('Type','line'),'LineWidth',4.0)
x=linspace(xmin,xmax);
y=linspace(ymin,ymax);
[X,Y]=meshgrid(x,y);
[m,sx]=size(x);
[m,sy]=size(y);
Z=zeros(sy,sx);
for i=1:sx
    for j=1:sy
        Z(j,i)=f([x(i);y(j)]);
    end
end
zmax=max(max(Z));
zmin=min(min(Z));
zmin=zmin-0.5*(zmax-zmin);
axis([xmin xmax ymin ymax zmin inf]);
if iplot==1
    contour(X,Y,Z)
elseif iplot==2
    view(45, 20);
    surfc(X,Y,Z)
else
    fprintf('[uo_plot] error: iplot = %i, <> 1,2\n',iplot);
end
xlabel('x_1');ylabel('x_2');zlabel('f(x)');title('');
[n,iter] = size(xk);
zl=zlim;
z(1:iter)=zl(1);
for i=1:iter
    zf(i)=f(xk(:,i))+0.0*abs(f(xk(:,i)));
end
plot3( xk(1,:)',  xk(2,:)', z(:), ':ok','LineWidth',3)
%plot3( xk(1,:)',  xk(2,:)', zf(:), ':oy','LineWidth',6)
if iplot==2
    for j=1:iter-1
        xx = linspace(xk(1,j),xk(1,j+1))';
        yy = linspace(xk(2,j),xk(2,j+1))';
        for i=1:size(xx,1) ff(i) = f([xx(i);yy(i)]); end
        plot3(xx',yy',ff, 'y','LineWidth',4)
        plot3([xk(1,j),xk(1,j)],[xk(2,j),xk(2,j)],[z(j),zf(j)],'--k')
    end
end
plot3([xk(1,iter),xk(1,iter)],[xk(2,iter),xk(2,iter)],[z(iter),zf(iter)],'--k')
quiver3( xk(1,1:iter)', xk(2,1:iter)', z(:), gk(1,1:iter)', gk(2,1:iter)', zeros(iter,1), 'LineWidth',2,'Color','r','Autoscale','on')
% End [uo_solve_plot] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

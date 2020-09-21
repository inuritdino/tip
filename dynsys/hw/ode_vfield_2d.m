function ode_vfield_2d(f,np,domain)
% Generate and plot the vector field of the 2D ODE given the RHS
% calculating function.
% USAGE:
%     ODE_VFIELD_2D(F,NP,DOMAIN)
%
% F: RHS calculating function, taking 2d vector x = [x(1); x(2)] and
% producing 2d output vector. Alternatively, if the system is linear F can
% be a 2x2 matrix such that x' = F*x.
%
% NP: number of initial points to generate. Default is 10.
%
% DOMAIN: space domain detalization. Default DOMAIN is -2:.2:2.


if(nargin < 3 || isempty(domain))
    domain = -2:.2:2;
    if(nargin < 2)
        num_poi = 10;
    else
        num_poi = np;
    end
else
    num_poi = np;
end
domain = linspace(min(domain),max(domain),20);

mesh_space = domain;
[x,y] = meshgrid(mesh_space);
% Preallocating. NOTE x and y have the same size
u = zeros(size(x));
v = zeros(size(x));
for ii=1:numel(x)
    if(isa(f,'function_handle'))
        dy = f([x(ii);y(ii)]);
        u(ii) = dy(1);
        v(ii) = dy(2);
    else
        dy = f*[x(ii);y(ii)];
        u(ii) = dy(1);
        v(ii) = dy(2);
    end
end

%num_poi = 10;
startx = (max(mesh_space)-min(mesh_space))*rand(1,num_poi) + min(mesh_space);
starty = (max(mesh_space)-min(mesh_space))*rand(1,num_poi) + min(mesh_space);
%startx = 10;
%starty = 10;
% startx = startx(abs(startx) > 0.5);
% starty = starty(abs(starty) > 0.5);
% min_num_el = min(numel(startx),numel(starty));
% startx = startx(1:min_num_el);
% starty = starty(1:min_num_el);

% Find steady state
% if(isa(f,'function_handle'))
%     barX = fsolve(f,...
%         [min(domain) + (max(domain) - min(domain))*rand;...
%         min(domain) + (max(domain) - min(domain))*rand]);
% else
%     barX = [0; 0];
% end

clf;
streamline(x,y,u,v,startx,starty);
grid on;
xlabel('X');
ylabel('Y');
hold on;
% Plot the steady state
%plot(barX(1),barX(2),'or','MarkerSize',10,'Linewidth',2);
plot(startx,starty,'*g','MarkerSize',6);
quiver(x,y,u,v,'-k');
hold off;
axis tight equal;
end
function [T,Y] = lorenz( sig,r,b,x0,tf )
% Simulates the Lorenz system using ode45 integrator.
% Lorenz system:
%                x' = sigma*(y - x)
%                y' = r*x - y - x*z
%                z' = x*y - b*z
% USAGE:
%     [T,Y] = LORENZ(SIG,R,B,X0,TF)
%
% SIG: sigma parameter.
% R: r parameter.
% B: b parameter.
% X0: vector (1x3) of the initial conditions (default is [0 2 0]).
% TF: final time until which to simulate (default is 100).
% T: the output time vector
% Y: the matrix with each column being a timeseries of the corresponding
% variable of the system.
%
% SEE ALSO: ode45

if(nargin < 5)
    tf = 100;
    if(nargin < 4)
        x0 = [0 2 0];
    end
end

% Euler method
% dt = 0.005;
% t = 0;
% ii = 1;
% n = floor(tf/dt);
% Y = zeros(n,3);
% T = zeros(n,1);
% Y(1,:) = x0';
% while(t <= tf)
%     T(ii) = t;
%     Y(ii+1,1) = Y(ii,1) + dt*sig*(Y(ii,2) - Y(ii,1));
%     Y(ii+1,2) = Y(ii,2) + dt*(r*Y(ii,1) - Y(ii,2) - Y(ii,1)*Y(ii,3));
%     Y(ii+1,3) = Y(ii,3) + dt*(Y(ii,1)*Y(ii,2) - b*Y(ii,3));
%     ii = ii + 1;
%     t = t + dt;
% end

[T,Y] = ode45(@lorenz_rhs,[0 tf],x0);
%plot(T,Y(:,1));
plot3(Y(:,1),Y(:,2),Y(:,3),'k');grid on;
xlabel('X');ylabel('Y');zlabel('Z');
set(gca,'FontSize',25);

    function dy = lorenz_rhs(t,y)
        dy = zeros(3,1);
        dy(1) = sig * (y(2) - y(1));
        dy(2) = r * y(1) - y(2) - y(1)*y(3);
        dy(3) = y(1)*y(2) - b*y(3);
    end

end


function x = run_map(f,x0,ni)
% Run a 1D discrete map defined by function handle F.
% USAGE:
%       X = RUN_MAP(F,X0,NI)
% F: the function, e.g. f = @(x)(l*x(1-x))
% X0: initial value (default 0)
% NI: number of iterations (default 10)
% X: the output time series.
if(nargin < 3)
    ni = 10;
    if(nargin < 2)
        x0 = 0;
    end
end
x = zeros(ni+1,1);
x(1) = x0;
diag = [];
for ii = 1:ni
    x(ii+1) = f(x(ii));
    diag = cat(1,diag,[x(ii) x(ii); x(ii) x(ii+1)]);
end

LOW = max(x) - min(x);
if(LOW > 0.001)
    DX = 0.1*LOW;
else
    DX = 0.5;
end
LOW = min(x) - DX;
HI = max(x) + DX;

subplot(121);
plot(0:ni,x,'o-');grid on;
xlabel('Time');
ylabel('X');
ylim([LOW HI]);
set(gca,'FontSize',20);

subplot(122);
xplot = linspace(LOW,HI,300);
plot(xplot,f(xplot),'-r');hold on;grid on;
text(LOW+(HI-LOW)/2,f(LOW+(HI-LOW)/2),'Y = F(X)','FontSize',15);
plot(xplot,xplot,'-m');
text(LOW+(HI-LOW)/2,LOW+(HI-LOW)/2,'Y = X','FontSize',15);
plot(diag(1,1),diag(1,2),'*g');
plot(diag(:,1),diag(:,2),'-');
plot(diag(end,1),diag(end,2),'*k');hold off;
xlim([LOW HI]);
ylim([LOW HI]);
%legend('y = f(x)','y = x','Initial point','Steps','Last point','Location','NorthWest');
set(gca,'FontSize',20);
axis equal;
end
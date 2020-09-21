% Integration of Hodgkin-Huxley equations with Euler method
% Katri Holm
% modified from http://users.cs.dal.ca/~tt/fcns/fcns_programs/spikes/hh.m

clear;

% maximal conductances (in units of mS/cm^2);
g(1) = 36;    % K
g(2) = 120;   % Na
g(3) = 0.3;   % leak

% battery voltages (in mV);
E(1) = -12;     % K
E(2) = 115;     % Na
E(3) = 10.613;  % V_leak

% Initialization of some variables
I_app = 0;
V = 0;
Cm = 1; % ?F/cm2
x = zeros(1,3); % x includes n, m, h
x(3)=1;

% Time stuff
dt=0.01; % time step for integration
t_stop = 20; % stop time

t_rec=0; % variable for recording

% Integration with Euler method
% start time -30 s --> system in steady-state at 0 s.
for t=-30:dt:t_stop

    % external current applications 
    if t >= 2 && t < 2.5   
        I_app = 5;
    end
            
    if t >= 10 && t < 10.5
        I_app = 20;
        
    end
    
    if t == 2.5 || t == 10.5
        I_app = 0;
    end

    % alpha functions used by Hodgkin-and Huxley
    alpha(1) = (10-V)/(100*(exp((10-V)/10)-1)); % alpha_n
    alpha(2) = (25-V)/(10*(exp((25-V)/10)-1));  % alpha_m
    alpha(3) = 0.07*exp(-V/20);                 % alpha_h

    % beta functions used by Hodgkin-and Huxley
    beta(1) = 0.125*exp(-V/80); % beta_n
    beta(2) = 4*exp(-V/18);     % beta_m
    beta(3) = 1/(exp((30-V)/10)+1); % beta_h

    % Tau_x and x_inf are defined with alpha and beta
    tau = 1./(alpha+beta);
    x_inf = alpha.*tau;

    % integration with Euler method
    x = x + ((x_inf - x)./tau).*dt;

    % calculate actual conductances g with given n, m, h
    gnmh(1) = g(1)*x(1)^4;      % gK
    gnmh(2) = g(2)*x(2)^3*x(3); % gNa
    gnmh(3) = g(3);             % gleak

    % Ohm's law
    I = ???

    % update voltage of membrane
    V = V + dt*((-sum(I) + I_app) / Cm);

    % record some variables after equilibration (t = 0) for plotting 
    if t >= 0;
        t_rec = t_rec+1;
        time(t_rec) = t;
        v(t_rec) = V;
        iapp(t_rec) = I_app;
        n(t_rec) = x(1);
        m(t_rec) = x(2);
        h(t_rec) = x(3);
        gK(t_rec) = gnmh(1);
        gNa(t_rec) = gnmh(2);
        IK(t_rec) = I(1);
        INa(t_rec) = I(2);
        IKNaleak(t_rec) = sum(I);

    end
end

figure;
subplot(6,1,1:2);
plot(time, v);
title('Membrane voltage (mV)');

subplot(6,1,3);
plot(time, n, time, m, time, h);
title('Gating variables');
legend('n(t)', 'm(t)', 'h(t)')

subplot(6,1,4);
plot(time, gK, time, gNa);
title('Conductances (mS/cm^2)');
legend('g_{Na}(t)', 'g_{K}(t)');

subplot(6,1,5);
plot(time, IK, time, INa, time, IKNaleak);
title('Currents (\muA/cm^2)');
legend('I_K(t)', 'I_{Na}(t)', 'I_K+I_{Na}+I_{leak}');

subplot(6,1,6);
plot(time, iapp);
title('Applied current (\muA/cm^2)');
xlabel('time (s)');
axis([0 max(time) -5 max(iapp)+5]);



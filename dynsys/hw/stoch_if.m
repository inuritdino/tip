 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Simulation of (leaky) integrate-and-fire neuron
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function stoch_if
 clear; clf;

 % Parameters of the model
 T     = 450;                       % Final simulation time
 dt    = 10^-2;                     % Time step
 N     = floor(T / dt) + 1;         % Number of points
 Tau_m = 10;                        % Membrane time constant
 u     = 2;                         % Initial membrane potential
 t     = 0:dt:T;                    % Time vector
 A_I_ext = 8.0;                     % Amplitude of the external input
 I_ext = A_I_ext * ones(1,N);       % Constant external input
%  f     = 0.02;                     % Frequency (Hz)
%  I_ext = 6 + sin(2*pi*t*f);        % Fluctuating external input
 R     = 2;                         % Resistance
 u_threshold = 12;                  % Firing threshold
 u_res = u;                         % Reset membrane potential

 U     = zeros(1,N);                % Vector for membrane potentials
 S     = zeros(1,N);                % spikes
 U(1)  = u;                         % Initial state
 sigma = 0;                         % Variance 
 
 % Integration with Euler method
 for i = 2:N;
     % Check the treshold
     if(u >= u_threshold)
         u = u_res;
         S(i - 1) = 1;
     end;
     % Integration
     u = (1 - (1 / Tau_m) * dt) * u + ((R * I_ext(i)) / Tau_m) * dt + sigma * randn *sqrt(dt);
     % Record
     U(i) = u;
 end

 subplot('position',[0.13 0.25 1-0.26 0.5])
   plot(t,U);
   hold on; plot([0 T],[u_threshold u_threshold],'--');
   axis([0 T 0 u_threshold+1])
   ylabel('Membrane potential u(t)')

 subplot('position',[0.13 0.8 1-0.26 0.1])
   plot(t,S,'.','markersize',20);
   axis([0 T 0.5 1.5])
   set(gca,'xtick',[],'ytick',[])
   ylabel('spikes')
   
 subplot('position',[0.13 0.1 1-0.26 0.08]);
   plot(t,I_ext);
    ylabel('I_{ext}');
    xlabel('time [\tau]');

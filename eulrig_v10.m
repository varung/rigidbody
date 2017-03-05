function [ output_args ] = rigid2( input_args )
%RIGID2 Simulates rigid body motion
%  displays plot showing the location of the angular velocity vector
%  in a frame fixed in the rotating body

I1 = 3; I2 = 3; I3 = 4.9;
I = [I1; I2; I3;]; % we only need diagonal terms

% set up initial conditions
W = [ 1; 0; 5;]; 

% main loop
tmax = 1000;
dt = .01;
k = ones( 3, 4 );

circ_t = linspace( 0, 2*pi);
circ_x = cos(circ_t);
circ_y = sin(circ_t);

% plot analytical answer
plot( circ_x, circ_y );
axis( [-5 5 -5 5 -5 5 ] );
axis square;
xlabel('X');
ylabel('Y');
zlabel('Z');
view(2);
Wobj = line('color','r','linestyle','.','markersize',15,...
    'erase','xor','xdata',[],'ydata',[]);

for t = 1:tmax

    % t doesn't matter now, but I was considering having a 
    % time dependent torque
    k(:,1) = dt * getWdot( W, t, I );
    k(:,2) = dt * getWdot( W+.5*k(:,1), t+.5*dt, I);
    k(:,3) = dt * getWdot( W+.5*k(:,2), t+.5*dt, I );
    k(:,4) = dt * getWdot( W+k(:,3), t + dt, I );
    
    W = W + (1/6) * ( k * [1;2;2;1]); 
    
    % plot
    hold on
    set(Wobj,'xdata',W(1),'ydata',W(2) );
    hold off
    drawnow
end

function [Wdot] = getWdot( W, T, I )

Wdot=[ (I(2) - I(3))/I(1) * W(2)*W(3);
    (I(3) - I(1))/I(2) * W(3)*W(1);
    (I(1) - I(2))/I(3) * W(1)*W(2); ];

function [ output_args ] = teul_1( testNum )
%TEUL_1 runs a demo of eulrig to illustrate comparison between analytical
%and true answer
% Tests:
% 1. Two Equal Moments, Zero Torque
%

switch testNum
    case {1,'TwoEqualMoments'}
        subplot(1,2,1);
        subplot(1,2,2);
        I = [2 2 4]';
        W = [1 0 1]';
        % plot analytical answer
        omega = (I(3) - I(1))/I(1) * W(3)
        circ_t = linspace( 0, 2*pi/omega);
        circ_x = W(1) * cos(omega*circ_t);
        circ_y =  W(1) * sin(omega*circ_t);
        circ_z = repmat( W(3), size(circ_x) ); 
        time_to_orbit = 2*pi / omega
        plot3( circ_x, circ_y, circ_z );
        axis square;
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        view(3);
        drawnow
        pause
        [NumSol] = eulrig( I, W);
        Xact = W(1) * cos( omega*NumSol(:,1) );
        Yact = W(1) * sin( omega*NumSol(:,1));
        Zact = W(3) * ones(size(Yact));
        ActSol = [ Xact Yact Zact];
        figure;
        hold on;
        diff = ActSol - NumSol(:,[2,3,4]); 
        plot( diff(:,1), diff(:,2) );
        
        
    case {2,'stable'},
        subplot(1,2,1);
        subplot(1,2,2);
        I = [2 8 4]';
        W = [0 1 .1]';
        % plot analytical answer
        omega = (I(3) - I(1))/I(1) * W(3)
        circ_t = linspace( 0, 2*pi/omega);
        circ_x = W(1) * cos(omega*circ_t);
        circ_y = W(1) * sin(omega*circ_t);
        circ_z = repmat( W(3), size(circ_x) ); 
        time_to_orbit = 2*pi / omega
        plot3( circ_x, circ_y, circ_z );
        axis square;
        grid on;
        
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        view(3);
        eulrig( I, W)
   case {3,'unstable'},
        subplot(1,2,1);
        subplot(1,2,2);
        I = [2 8 4]';
        W = [0 .1 1]';
        % plot analytical answer
        omega = (I(3) - I(1))/I(1) * W(3)
        circ_t = linspace( 0, 2*pi/omega);
        circ_x = W(1) * cos(omega*circ_t);
        circ_y = W(1) * sin(omega*circ_t);
        circ_z = repmat( W(3), size(circ_x) ); 
        time_to_orbit = 2*pi / omega
        plot3( circ_x, circ_y, circ_z );
        axis square;
        grid on;
        
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        view(3);
        eulrig( I, W)
    case {4, 'pureaxis'}
        I = [8 8 8]';
        W = [0 0 1]';
        torque(I, W, [ 0 0 0]' );
    case {5, 'top' }
        I = [2 2 8]';
        W = [0 0 3]';
        T = [0 0 8]';
        torque( I, W, T );
    case {6, 'toptf'}
        I = [2 2 8]';
        W = [.5 0 1]';
        eulrig( I, W );
end
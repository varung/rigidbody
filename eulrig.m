function [ SavedW ] = eulrig( I_in, W_in  )
%RIGID2 Simulates rigid body motion
%  displays plot showing the location of the angular velocity vector
%  in a frame fixed in the rotating body

I1 =I_in(1); I2 = I_in(2); I3 = I_in(3);
I = [I1; I2; I3;]; % we only need diagonal terms

% set up initial conditions
W = W_in; 
Q = (angToQuat( pi/4, [1 0 0] ) )';

% main loop
tmax = 1000;
dt = .1;
k = ones( 3, 4 );

% init display
X = zeros(1,3);
k = 5;
ns = 2^k-1; 
c = hadamard(2^k);
[xd,yd,zd]=ELLIPSOID(X(1),X(2),X(3),I1,I2,I3,ns);
[md nd] = size(xd);

subplot(1,2,1);
subplot(1,2,2);
Wmarker = line( W(1), W(2), W(3), 'Marker','.', 'MarkerSize',5);
L = diag(I) * W;
Lmarker = line( L(1), L(2), L(3), 'Marker','.', 'MarkerSize',10, 'Color', [1 0 0]);

str_title = sprintf('Time: 0');
h_title = title(str_title);
maxd=0;
SavedW = zeros( tmax, 4 );
SavedTopNum = zeros( tmax, 4 );
for t = 1:tmax
    
    % R-K4 Integrate W to find Wnew in body frame
    k_q(:,1) = dt * getQdot( W, Q );
    k_w(:,1) = dt * getWdot( W, t, I );
    
    k_q(:,2) = dt * getQdot( W+.5*k_w(:,1), Q+.5*k_q(:,1) );
    k_w(:,2) = dt * getWdot( W+.5*k_w(:,1), t+.5*dt, I);
    
    k_q(:,3) = dt * getQdot( W+.5*k_w(:,2), Q+.5*k_q(:,2) );
    k_w(:,3) = dt * getWdot( W+.5*k_w(:,2), t+.5*dt, I );
    
    k_q(:,4) = dt * getQdot( W+k_w(:,3), Q+k_q(:,3) );
    k_w(:,4) = dt * getWdot( W+k_w(:,3), t + dt, I );
    
    Q = Q + (1/6) * ( k_q * [1;2;2;1]);
    W = W + (1/6) * ( k_w * [1;2;2;1]);
    
    % renormalize q
    Q = Q ./ sqrt(dot(Q,Q)); 
    SavedW( t, : ) = [ t*dt, W' ];
    
    % plot only every 100 time steps
    if( mod(t,100) )
        % show sphere
        hold on
        subplot(1,2,1);
        points = [ xd(:), yd(:), zd(:) ];
        R = quatToMat( Q );
        points = (R * points')';
        x = reshape( points(:,1), md, nd );
        y = reshape( points(:,2), md, nd );
        z = reshape( points(:,3), md, nd );
        surf(x,y,z,c);
        colormap([1  1  0; 0  1  1])
        axis( [-10 10 -10 10 -10 10] );
        Xlabel('X');
        Ylabel('Y');
        Zlabel('Z');
        axis square;

        subplot(1,2,2);
        L = diag(I) * W;
        L = L ./ norm(L);
        L = R * L;
        TopNum = R * [ 0 0 1]'; 
        SavedTopNum(t,:) = [t*dt, TopNum' ];
        %set(Wmarker,'xdata',W(1),'ydata',W(2), 'zdata', W(3) );
        set(Lmarker,'xdata',L(1),'ydata',L(2), 'zdata', L(3) );
        set(Wmarker,'xdata',SavedTopNum(:,2),'ydata',SavedTopNum(:,3), 'zdata', SavedTopNum(:,4) );
            
        str_title = sprintf('Time: %3.3f', t*dt );
        set(h_title, 'String', str_title );
        maxd = max( 1.1 * max( [abs(W); abs(L)] ), maxd );
        axis([-maxd maxd -maxd maxd -maxd maxd]);
        Xlabel('X_bdy');
        Ylabel('Y_bdy');
        axis square;
        grid on;
        view(3);
        hold off
%       pause(.05);
        drawnow;
    end
end

disp ( 'Done!');

% (dL/dT)inl = (dL/dT)bdy + cross( W, L)
% 0 = IWdot + cross(W, L )
% L = I*W
% I*Wdot = cross( IW, W ) 
% Wdot = I / cross( IW, W );

function [Wdot] = getWdot( W, T, I )
Imat = diag( I);
Imat_inv = diag( 1 ./ I );
Wdot = Imat_inv * cross( Imat * W, W ) ;


% Quaternion subroutines
% Quaternions represent orientation in 4 variables rather than 9 which is
% what keeping the rotation matrix would do
% Moreover, quaternions avoid gimbal lock which can spoil integation with
% euler angles

% converts angular velocity + quaternion to derivative of quaternion
function [qdot] = getQdot( w, q )
% w is defined in the body frame, so transform it to the world frame
R = quatToMat( q );
w_inl = R * w;   

% calculate derivative of quaternion due to angular velocity
qdot = .5 * quatMult( ([0 w_inl']), q');
qdot = qdot';
% converts angle, axis to quaternion
function [quat] = angToQuat( angle, axis )
quat(1) = cos( angle / 2);
quat(2:4) = axis * sin(angle/2);

% converts quaternion to rotation matrix
% q = [s v] where s = scalar, v = vector
% [cos( theta/2) sin(theta/2)*u] where u = unit axis about which 
% the rotation takes place, theta = angle of rotation about axis
function [rot_m] = quatToMat( q )
s = q(1); vx = q(2); vy = q(3); vz = q(4);
rot_m = [ (1 - 2*vy*vy - 2*vz*vz), (2*vx*vy - 2*s*vz), (2*vx*vz + 2*s*vy);
          2*vx*vy + 2*s*vz, 1 - 2*vx*vx - 2*vz*vz, 2*vy*vz - 2*s*vx;
          2*vx*vz - 2*s*vy, 2*vy*vz + 2*s*vx, 1 - 2*vx*vx - 2*vy*vy];

% multiply quaternions together
% qr = q1 * q2
function [qr] = quatMult(q1, q2)
s1 = q1(1); v1 = q1(2:4);
s2 = q2(1); v2 = q2(2:4);
qr = [ s1*s2 - dot(v1,v2), s1*v2 + s2*v1 + cross( v1,v2) ];
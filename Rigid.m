function [ output_args ] = Rigid()
%RIGID Summary of this function goes here
%  Detailed explanation goes here

% inertia tensor, diagonal matrix since we can set body axes to principal
% axes
IA = 4; IB = 4; IC = 8;

I_body = eye(3);
I_body(1,1) = IA;
I_body(2,2) = IB;
I_vody(3,3) = IC;
I_body_inv = inv(I_body)

% state variables
X = [0 0 0]; % initial position
Q = angToQuat( 0, [0 0 1] );
P = [0 0 0]; % intial linear momentum
W = [1 0 5 ]; % set initial ang velocity
L = (I_body * W')';

N = 1000;
dt = .01;
stuff = zeros( N, 4 );

% init display
k = 5;
ns = 2^k-1; 
c = hadamard(2^k);
[xd,yd,zd]=ellipsoid(X(1),X(2),X(3),IA,IB,IC,ns);
[md nd] = size(xd);
for time = 1:N
    stuff(time,:) = Q;
    R = quatToMat( Q );
    % find Inertia matrix in the inertial frame
    % I_inv is similar to I_body_inv, just defined in a different basis
    I_inv = R * I_body_inv * R';
    
    % L = I*W, since L is constant, just keep L as the state variable
    % and solve for W:
    W = I_inv * L';

    % update each state variable
    % no torque, no velocity
    qdot = getQdot( W', Q );

    % euler
    Q = Q + qdot*dt;
    
    % renormalize
    Q = Q ./ sqrt(dot(Q,Q)); 
    
    % show sphere
    hold on
    subplot(1,2,1);
    points = [ xd(:), yd(:), zd(:) ];
    points = (R * points')';
    x = reshape( points(:,1), md, nd );
    y = reshape( points(:,2), md, nd );
    z = reshape( points(:,3), md, nd );
    surf(x,y,z,c);
    colormap([1  1  0; 0  1  1])
    axis( [-10 10 -10 10 -10 10] );
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis square;
    subplot(1,2,2);
    Wbody = R' * W; % rotate from world to body coordinates
    line( Wbody(1), Wbody(2), Wbody(3), 'Marker','+', 'MarkerSize',10);
    axis([-5 5 -5 5 -5 5]);
    xlabel('X_bdy');
    ylabel('Y_bdy');
    view(2);
    axis square;
    hold off
    drawnow;
end


% Quaternion subroutines
% Quaternions represent orientation in 4 variables rather than 9 which is
% what keeping the rotation matrix would do
% Moreover, quaternions avoid gimbal lock which can spoil integation with
% euler angles

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

% converts angular velocity + quaternion to derivative of quaternion
function [qdot] = getQdot( w, q )
qdot = .5 * quatMult( [0 w], q);

% multiply quaternions together
% qr = q1 * q2
function [qr] = quatMult(q1, q2)
s1 = q1(1); v1 = q1(2:4);
s2 = q2(1); v2 = q2(2:4);
qr = [ s1*s2 - dot(v1,v2), s1*v2 + s2*v1 + cross( v1,v2) ];

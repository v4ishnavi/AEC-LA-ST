%Given Constants
m1=10;
m2=5;
l1=0.2;
l2=0.1;
g=9.81;

%Initial State Variables
q1=0;
q2=0;
q3=0.1;
q4=0.1;
q3_dot = 0;
q4_dot = 0;
init_q=[q1,q2, q3, q4, q3_dot, q4_dot];

%PI controller gains
kp=750;
ki=1;
gain=[kp;ki];

%Final Angles
qf=[0;0];

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
[time,s] = ode45(@(t,q) pid(t,q,qf,gain), [0, 30], init_q,options);

q_final=zeros(length(time),1);

figure;
subplot(2,1,1);
plot(time,s(:,3));
axis([0,30,-1,1]);
xlabel("Time");
ylabel("q_1");
hold on;
plot(time,q_final,'g');
hold off;

subplot(2,1,2);
plot(time,s(:,4),'r');
axis([0,30,-1,1]);
xlabel("Time");
ylabel("q_2");
hold on;
plot(time,q_final,'g');
hold off;

sgtitle("Joint Angle vs Time - PI Controller");

function [q_double] = pid(t,q,qf,gains)
    %Given Constants
    m1=10;
    m2=5;
    l1=0.2;
    l2=0.1;
    g=9.81;

   %Error
    error = [qf(1) - q(3);
         qf(2) - q(4)];

    f1 = gains(1) * error(1) + gains(2)*q(1);
    f2 = gains(1) * error(2) + gains(2)*q(2);
    f = [f1; f2];

    %Mass Matrix
    M11 = (m1 + m2) * (l1^2) + m2 * l2*(l2 + 2 * l1 * cos(q(4)));
    M22 = m2 * l2^2;
    M12 = m2*l2*(l2+l1*cos(q(4)));
    M21 = M12;
    M = [M11, M12; M21, M22];

    %Coriolis Matrix
    c11 = -m2 * l1 * l2 * sin(q(4)) * q(5);
    c12 = -m2 * l1 * l2 * sin(q(4)) * (q(5) + q(6));
    c21 = 0;
    c22 = m2 * l1 * l2 * sin(q(4)) * q(6);
    C = [c11, c12; c21, c22];

    %Gravitational Matrix
    G1 = m1*l1*g*cos(q(3)) + m2*g*(l2*cos(q(3)+q(4)) + l1*cos(q(3)));
    G2 = m2 * l2 * g * cos(q(3)+q(4));
    G = [G1;G2];

    dq = [q(5); q(6)];
    
    m_inverse = inv(M);
    W = ((m_inverse) * (-C * dq - G)) + f;

    %Next state
    q_double=zeros(4,1);
    q_double(1)=qf(1)-q(3);
    q_double(2)=qf(2)-q(4);
    q_double(3)=q(5);
    q_double(4)=q(6);
    q_double=[q_double;W];
end


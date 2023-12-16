%Given Constants
m1=10;
m2=5;
l1=0.2;
l2=0.1;
g=9.81;

%Initial State Variables
q1=0.1;
q2=0.1;
q1_dot = 0;
q2_dot = 0;
init_q=[q1,q2, q1_dot, q2_dot];

%PD controller Gains
kp=100;
kd=100;
gain=[kp;kd];

%Final Angles
qf=[0;0];

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
[time,s] = ode45(@(t,q) pd(t,q,qf,gain), [0, 30], init_q,options);

figure;
subplot(2,1,1);
plot(time,s(:,1));
xlabel("Time");
ylabel("q_1");

subplot(2,1,2);
plot(time,s(:,2),'r');
xlabel("Time");
ylabel("q_2");

sgtitle("Joint Angle vs Time - PD Controller");

function [q_double] = pd(t,q,qf,gains)
    %Given Constants
    m1=10;
    m2=5;
    l1=0.2;
    l2=0.1;
    g=9.81;

    %Error
    error = [qf(1) - q(1);
         qf(2) - q(2)];

    f1 = gains(1) * error(1) - gains(2) * q(3);
    f2 = gains(1) * error(2) - gains(2) * q(4);


    %Mass Matrix
    M11 = (m1 + m2) * (l1^2) + m2 * l2*(l2 + 2 * l1 * cos(q(2)));
    M22 = m2 * l2^2;
    M12 = m2*l2*(l2+l1*cos(q(2)));
    M21 = M12;
    M = [M11, M12; M21, M22];

    %Coriolis Matrix
    c11 = -m2 * l1 * l2 * sin(q(2)) * q(3);
    c12 = -m2 * l1 * l2 * sin(q(2)) * (q(3) + q(4));
    c21 = 0;
    c22 = m2 * l1 * l2 * sin(q(2)) * q(4);
    C = [c11, c12; c21, c22];

    %Gravitational Matrix
    G1 = m1*l1*g*cos(q(1)) + m2*g*(l2*cos(q(1)+q(2)) + l1*cos(q(1)));
    G2 = m2 * l2 * g * cos(q(1)+q(2));
    G = [G1;G2];

    f = [f1; f2];
    dq = [q(3); q(4)];
    
    m_inverse = inv(M);
    W = ((m_inverse) * (-C * dq - G)) + f;

    %Next state
    q_double=zeros(2,1);
    q_double(1)=q(3);
    q_double(2)=q(4);
    q_double=[q_double;W];
end


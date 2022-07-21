%% WECC 3 machine 9 bus system
m1 = [2, 2, 2];
m2 = [8, 8, 8];
d = [1, 1, 1]; 
Y = 0.1*[0.8455 - 2.9883i   0.2871 + 1.5129i   0.2096 + 1.2256i
   0.2871 + 1.5129i   0.4200 - 2.7239i   0.2133 + 1.0879i
   0.2096 + 1.2256i   0.2133 + 1.0879i   0.2770 - 2.3681i];

L = -(Y - real(Y));
L = L - diag(diag(L));
L = imag((L - diag(sum(L,2)))); % B prime matrix without removing slack


%% Model
rand('seed',1)
ngen = 3;
d = 2*ngen;
Q = 0.2*eye(2 * ngen);
R = eye(2 * ngen);

alpha = 0.6;

wb = 376.9911;
Ts = 0.1;
T = 300;

t = 1;
numiter = 300;
while t <= numiter
    disp(t)
    if(mod(t,4)==1|mod(t,4)==2)
        m = randsample([2, 2.05, 2.1, 2.15, 2.2], 1); %randsample([2, 2.05, 2.1, 2.15, 2.2], 1); %2
    else
        m = randsample([8, 8.05, 8.1, 8.15, 8.2], 1); %8; %8; %+0.01*rand(1); %7.95 + 0.1*rand(1); %8
    end
    M = diag([m, m, m]);
    D = diag(d);

    % Build A matrix
    A_c = [zeros(ngen), wb .* eye(ngen);
        -inv(M) * L, -inv(M) * D];
    A = expm(A_c * Ts);

    % Build Bdd matrix
    B_c = [zeros(ngen); -inv(M)];
    B = inv(A_c) * (expm(A_c * Ts) - eye(2 * ngen)) * B_c;

    if(~isnan(B))
        if(t==1)
            A_list = A;
            B_list = B;
            m_list = m;
        else
            A_list = [A_list; A];
            B_list = [B_list; B];
            m_list = [m_list; m];
        end
        t = t+1;
    end
        
end
    
w_mean = zeros(1,6);
w_cov = 1.0000e-014*eye(6);
w = mvnrnd(w_mean,w_cov,T);

%MPC
x3 = zeros(2*ngen, T+1);
x3(:,1) = 0.000005*ones(6,1);
J3 = zeros(T,1);
for t = 1:T-1
    A_t = A_list(d*(t-1)+1:d*t,:);
    B_t = B_list(d*(t-1)+1:d*t,:);
    A_t1 = A_list(d*t+1:d*(t+1),:);
    B_t1 = B_list(d*t+1:d*(t+1),:);

    A_trans = A_t1*A_t;
    B_trans = [B_t1, A_t1*B_t];
    if(mod(t,2)==1)
        cvx_begin quiet sdp
            variable S11(d,d) symmetric
            variable S12(d,d) 
            variable S22(d,d) symmetric
            minimize (trace(Q'*S11+R'*S22))
            [A_trans B_trans]*[S11 S12; S12' S22]*[A_trans'; B_trans'] + w_cov==S11
            [S11 S12; S12' S22] >= 0
            [A_trans B_trans]*[S11 S12; S12' S22]*[A_trans'; B_trans'] <= alpha*S11
        cvx_end

        K = S12'*inv(S11);
        if(isnan(K))
            disp('COCO-LQ NaN');
            K = 0;
        end
        u3 = K*x3(:,t);
        act = u3(4:end);
    else
        act = u3(1:3);
    end

%     K = lqr(A_t, B_t, Q, R);
%     u = K*x3(:,t);
%     act = u;
    
    x3(:, t+1) = A_t*x3(:,t)+B_t*act; 
    if(t==1)
        J3(t,:) = x3(:, t)'*Q*x3(:, t)+act'*R(1:3,1:3)*act;
    else
        J3(t,:) = J3(t-1,:)+x3(:, t)'*Q*x3(:, t)+act'*R(1:3,1:3)*act;
    end
end

% LQR with predictions
x = zeros(2*ngen, T+1);
J = zeros(T,1);
x(:,1) = 0.000005*ones(6,1);
R2 = eye(ngen);
for t = 1:T-1
    A_t = A_list(d*(t-1)+1:d*t,:);
    B_t = B_list(d*(t-1)+1:d*t,:);
    A_t1 = A_list(d*t+1:d*(t+1),:);
    B_t1 = B_list(d*t+1:d*(t+1),:);
    
    x11 = x(1,t);
    x12 = x(2,t);
    x13 = x(3,t);
    x14 = x(4,t);
    x15 = x(5,t);
    x16 = x(6,t);
    if(mod(t,2)==1)
        cvx_begin quiet
            variables x21 x22 x23 x24 x25 x26 x31 x32 x33 x34 x35 x36 u11 u12 u13 u21 u22 u23
            minimize Q(1,1)*(square(x21)+square(x31))+Q(2,2)*(square(x22)+square(x32))...
                     +Q(3,3)*(square(x23)+square(x33))+Q(4,5)*(square(x24)+square(x34))...
                     +Q(5,5)*(square(x25)+square(x35))+Q(6,6)*(square(x26)+square(x36))...
                     +R(1,1)*(square(u11)+square(u21))+R(2,2)*(square(u12)+square(u22))+R(3,3)*(square(u13)+square(u23))
            [x21; x22; x23; x24; x25; x26] == A_t*[x11; x12; x13; x14; x15; x16]+B_t*[u11; u12; u13]
            [x31; x32; x33; x34; x35; x36] == A_t1*[x21; x22; x23; x24; x25; x26]+B_t1*[u21; u22; u23]
        cvx_end

        act = [u11; u12; u13];
    else
        act = [u21; u22; u23];
    end

%     K = lqr(A_t, B_t, Q, R);
%     u = K*x3(:,t);
%     act = u;
    
    x(:, t+1) = A_t*x(:,t)+B_t*act;   
    if(t==1)
        J(t,:) = x(:, t)'*Q*x(:, t)+act'*R(1:3,1:3)*act;
    else
        J(t,:) = J(t-1,:)+x(:, t)'*Q*x(:, t)+act'*R(1:3,1:3)*act;
    end
end

% offline optimal
P_list = zeros(6*T,6);
for t = T:-1:1
    if(t==T)
        P_list((t-1)*6+1:t*6,:) = Q;
    else
        A_t = A_list(d*(t-1)+1:d*t,:);
        B_t = B_list(d*(t-1)+1:d*t,:);
        P_t1 = P_list((t)*6+1:(t+1)*6,:);
        P = A_t'*P_t1*A_t+Q-A_t'*P_t1*B_t*inv(R(1:3,1:3)+B_t'*P_t1*B_t)*B_t'*P_t1*A_t;
        P_list((t-1)*6+1:t*6,:) = P;
    end
end

x_opt = zeros(2*ngen, T+1);
x_opt(:,1) = 0.000005*ones(6,1);
J_opt = zeros(T,1);

for t = 1:T-1
    A_t = A_list(d*(t-1)+1:d*t,:);
    B_t = B_list(d*(t-1)+1:d*t,:); 
    P_t1 = P_list((t)*6+1:(t+1)*6,:);
    x_optt = x_opt(:, t);
    K_optt = inv(R(1:3,1:3)+B_t'*P_t1*B_t)*B_t'*P_t1*A_t;
    
    u_optt = -K_optt*x_optt;
    
    x_opt(:, t+1) = A_t*x_opt(:,t)+B_t*u_optt;  
    if(t==1)
        J_opt(t,:) = x_opt(:, t)'*Q*x_opt(:, t)+u_optt'*R(1:3,1:3)*u_optt;
    else
        J_opt(t,:) = J_opt(t-1,:)+x_opt(:, t)'*Q*x_opt(:, t)+u_optt'*R(1:3,1:3)*u_optt;
    end
end


figure; hold on
x0=10;
y0=10;
width=600;
height=300;
set(gcf,'position',[x0,y0,width,height])
grid on
a1 = plot(60*x_opt(1,:), '-g', 'LineWidth',1.5); M1 = 'Offline Optima';
a2 = plot(60*x(1,:), '-b', 'LineWidth',1); M2 = 'H-horizon Control';
a3 = plot(60*x3(1,:), '-r', 'LineWidth',1); M3 = 'COCO-LQ-Prediction';

yline(0.05,'--','','LineWidth',1);
yline(-0.05,'--','','LineWidth',1);
legend([a1;a2;a3], M1, M2, M3, 'FontSize',16, 'Location','NorthWest')
xlabel('Time')
ylabel('Frequency Deviation [Hz]')
ylim([-0.5, 0.5])
xlim([0, 300])
ax = gca
% Set x and y font sizes.
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
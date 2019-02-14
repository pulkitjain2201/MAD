function [dydt]  = fun(t,y)
global N a
% S --> Ca--> Ts--> Ta--> X 
Bi = 2.22*a/0.21;
D_ag = .198;
e = 0.22;
C = 2.27;
d = 1.11;
M = 1783206.34;
E = 21.72;
p_er = a^2/0.417;
t_sa = 2906.48;
p_esr = (a^2)*1476.14;
t_as = 1.29;
g = 2.15;
p_e = a^2*0.007704;
f = 30.018;
% Intermediate value + values at R0
S(1:N+2) = y(1:N+2);
Ca(1:N+2) = y(N+3:2*N+4);
Ts(1:N+2) = y(2*N+5:3*N+6);
Ta(1:N+2) = y(3*N+7:4*N+8);
X(1:N+2) = y(4*N+9:5*N+10);
dydt = zeros(1,5*N+2);
%%
dydt(1) = y(2);
dydt(N+3) = y(N+4);
dydt(2*N+5) = y(2*N+6);
dydt(3*N+7) = y(3*N+8);
dydt(4*N+9) = y(4*N+10);
%%
if Ts(N+2)< 1.16
    alpha = 2.71;
else
    alpha = (2.71^(-17.23/Ts(N+2)))*2.71^(41250/(8.314*333));
end
dydt(2*N+4) = 1;%boundary value C
dydt(3*N+6) = ((N+1)*Ts(N+1) - Bi)/(N+1-Bi);%boundary value Ts
dydt(4*N+8) = 1; %boundary value Ta
dydt(5*N+10) = D_ag * alpha *X(N+2)*(1-e*(X(N+2)-1))/((C+1)*(d*X(N+2) + (1-e*(X(N+2)-1)))) - M*X(N+2)*(2.71)^(-E/Ts(N+2)); %boundary value X
dydt(N+2) = -e*dydt(5*N+10); %boundary value of S

%%
for i = 1:N
    if Ts(i+1)< 1.16
        alpha = (2.71^(-17.23/Ts(i+1)))*2.71^(41250/(8.314*333));
    else
        alpha = (1.2257-Ts(i+1))*14.4;
    end
    dydt(3*N+7+i) = ((2*(N+1)^2/i)*(Ta(i+2) - Ta(i+1)) + (Ta(i+2) + Ta(i) - 2*Ta(i+1))*(N+1)^2)/p_er - t_sa*(Ta(i+1)-Ts(i+1)); %for Ta
    dydt(4*N+9+i) = D_ag * alpha *X(i+1)*(1-e*(X(i+1)-1))* Ca(i+1)/((C+Ca(i+1))*(d*X(i+1) + (1-e*(X(i+1)-1)))) - M*X(i+2)*(2.71)^(-E/Ts(i+2)); %for X
    dydt(i+1) = -e*dydt(4*N+9+i); %for S
    dydt(2*N+5+i) =  ((2*(N+1)^2/i)*(Ts(i+2) - Ts(i+1)) + (Ts(i+2) + Ts(i) - 2*Ts(i+1))*(N+1)^2)/p_esr + t_as*(Ta(i+1)-Ts(i+1))-g*dydt(4*N+9+i);%for Ts
    dydt(N+3+i)  =  ((2*(N+1)^2/i)*(Ca(i+2) - Ca(i+1)) + (Ca(i+2) + Ca(i) - 2*Ca(i+1))*(N+1)^2)/p_e - f*dydt(4*N+9+i);% For C 
end
dydt
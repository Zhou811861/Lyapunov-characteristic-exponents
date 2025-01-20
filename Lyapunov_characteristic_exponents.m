%  Implemented algorithm presented by Benettin.G et al in 1980 to 
% compute Lyapunov characteristic exponents of 3D dynamical systems. 
% 
%  Computes Lyapunov exponents by evolving perturbations in a system's phase space, 
% using the Jacobian matrix, periodic normalization, and Gram-Schmidt orthogonalization. 
K1=0;K2=0;K3=0;
evolution=zeros(3,steps);
for t=1:steps
    [T,y]=ode45(@lorenz,[0 tstep],y0);
    y0=y(end,:)';
    s=norm(y0([4 7 10]));K1=K1+log(s)/log(2);
    y0([4 7 10])=y0([4 7 10])/s;
    y0([5 8 11])=y0([5 8 11])-(y0([4 7 10])*(y0([5 8 11]))')*y0([4 7 10]);
    s=norm(y0([5 8 11]));K2=K2+log(s)/log(2);
    y0([5 8 11])=y0([5 8 11])/s;
    y0([6 9 12])=y0([6 9 12])-(y0([4 7 10])*(y0([6 9 12])'))*y0([4 7 10])-(y0([5 8 11])*(y0([6 9 12]))')*y0([5 8 11]);
    s=norm(y0([6 9 12]));K3=K3+log(s)/log(2);
    y0([6 9 12])=y0([6 9 12])/s;
    evolution(1,t)=K1/(t*tstep);evolution(2,t)=K2/(t*tstep);evolution(3,t)=K3/(t*tstep);
end
K1=K1/(steps*tstep)
K2=K2/(steps*tstep)
K3=K3/(steps*tstep)
plot(1:steps,evolution,'LineWidth',3)
xlabel('time');ylabel('Lyapunov Exponents');
legend('K1','K2','K3');ylim([-35,15]);
%KY=2-K1/K3 %Kaplan-Yorke dimension

clc,clear,close all

% 1D
% parameters
xmin=-5;
xmax=5;
delta=0.5; % range of trial walk 
num_iter=1000000; % number of iterations inside one VMC
num_equi=500; % 500,after which number begin to calculate the average and variance
out_iter=1; % number of iterations that initialization of VMC
h=0.1; % discretize on x axis, step size,analytic solution
alpha=0.5;

ax=0;
bx=0.2;
cx=1;
tol=0.000001;

% main
% Task 1
% compare the expectation value between the VMC and analytic solution,the
% variance with the book, VMC takes 3s
% alpha=0:0.01:1;
alpha1=0.4:0.01:0.6;
% alpha1=0.1:0.1:1;
alpha2=0.05:0.01:1;
tic
[pos,E_L,E_average,E_variance,E_average_tot,E_variance_tot]=VMC_1d1(alpha,num_equi,xmin,xmax,num_iter,delta,out_iter);
toc
% plot histogram
% for i=1:length(alpha)
% [pos,E_L,E_average,E_variance,E_average_tot,E_variance_tot]=VMC_1d1(alpha(1,i),num_equi,xmin,xmax,num_iter,delta,out_iter);
% % subplot(2,2,i);
% yyaxis left
% histogram(pos(2:end,:));
% xlabel('z');
% ylabel('Number of visits through trial walk')
% title(['α = ',num2str(alpha(1,i))]);
% yyaxis right
% x=xmin:0.001:xmax;
% N=size(x);
% H1=ones(1,N(1,2));
% wave_function=exp(-alpha*x.^2);
% y_ana=exp(-alpha(1,i)*x.^2).^2/sqrt(pi);
% [integral_ana]=cal_integral(y_ana,h,H1);
% plot(x,y_ana,'LineWidth',4);
% ylabel('Probability density')
% end

% plot

% [pos,E_L,E_average,E_variance,E_average_tot, E_variance_tot]=VMC_1d1(alpha,num_equi,xmin,xmax,num_iter,delta,out_iter);
% [E_average_tot_1,E_variance_tot_1,E_expetation_1,E_average_tot_2,E_variance_tot_2,E_expetation_2]=iter_1d(alpha1,alpha2,num_equi,xmin,xmax,num_iter,delta,out_iter,h);
% toc
% figure()
% % plot energy and variance with different alpha
% subplot(1,2,1)
% yyaxis left
% plot(alpha2,E_expetation_2,'LineWidth',2);hold on;
% scatter(alpha2,E_average_tot_2,'filled','blue');
% % legend('Analytic solution of E_L','E_L')
% xlabel('α')
% ylabel('E [ℏω]')         
% yyaxis right
% scatter(alpha2,E_variance_tot_2,'filled');
% % xlabel('α')
% ylabel('Variance')
% 
% subplot(1,2,2)
% yyaxis left
% plot(alpha1,E_expetation_1,'LineWidth',2);hold on;
% scatter(alpha1,E_average_tot_1,'filled','blue');
% % legend('Analytic solution of E_L','E_L')
% xlabel('α')
% ylabel('E [ℏω]')
% yyaxis right
% scatter(alpha1,E_variance_tot_1,'filled');
% % xlabel('α')
% ylabel('Variance')

% Task 2
% use golden search to find the minimum
% tic
% [xmin_gs,fmin_gs]=golden_search(ax,bx,cx,tol,@VMC_1d,num_equi,xmin,xmax,num_iter,delta,out_iter);
% toc

function [E_average_tot_1,E_variance_tot_1,E_expetation_1,E_average_tot_2,E_variance_tot_2,E_expetation_2]=iter_1d(alpha1,alpha2,num_equi,xmin,xmax,num_iter,delta,out_iter,h)
E_average_tot_1=zeros(1,1);E_variance_tot_1=zeros(1,1);E_expetation_1=zeros(1,1);
for i=1:length(alpha1)
[pos,E_L,E_average,E_variance,E_average_tot1,E_variance_tot1]=VMC_1d1(alpha1(1,i),num_equi,xmin,xmax,num_iter,delta,out_iter);
[E_expetation1]=cal_ana(h,xmin,xmax,alpha1(1,i)); % analytic solution
E_expetation_1=[E_expetation_1,E_expetation1];
E_average_tot_1=[E_average_tot_1,E_average_tot1];
E_variance_tot_1=[E_variance_tot_1,E_variance_tot1];
end
E_average_tot_1=E_average_tot_1(1,2:end);
E_variance_tot_1=E_variance_tot_1(1,2:end);
E_expetation_1=E_expetation_1(1,2:end);

E_average_tot_2=zeros(1,1);E_variance_tot_2=zeros(1,1);E_expetation_2=zeros(1,1);
for j=1:length(alpha2)
[pos,E_L,E_average,E_variance,E_average_tot2,E_variance_tot2]=VMC_1d1(alpha2(1,j),num_equi,xmin,xmax,num_iter,delta,out_iter);
[E_expetation2]=cal_ana(h,xmin,xmax,alpha2(1,j)); % analytic solution
E_expetation_2=[E_expetation_2,E_expetation2];
E_average_tot_2=[E_average_tot_2,E_average_tot2];
E_variance_tot_2=[E_variance_tot_2,E_variance_tot2];
end
E_average_tot_2=E_average_tot_2(1,2:end);
E_variance_tot_2=E_variance_tot_2(1,2:end);
E_expetation_2=E_expetation_2(1,2:end);
end

% VMC,one particle 1D,original
function [pos,E_L,E_average,E_variance,E_average_tot,E_variance_tot]=VMC_1d1(alpha,num_equi,xmin,xmax,num_iter,delta,out_iter)
pos=zeros(1,num_iter+2);
E_L=zeros(1,num_iter+2);
E_average=zeros(1,1);
E_variance=zeros(1,1);
for i=1:out_iter
[pos1,E_L1,E_average1,E_variance1]=Quantum_Monte_carlo(num_equi,alpha,xmin,xmax,num_iter,delta);
pos=[pos;pos1];
E_L=[E_L;E_L1];
E_average=[E_average,E_average1];
E_variance=[E_variance,E_variance1];
end
E_average_tot=sum(E_average(1,2:end))/out_iter;
E_variance_tot=sum(E_variance(1,2:end))/out_iter;
end

% VMC,one particle 1D,for golden search
function [E_variance_tot]=VMC_1d(alpha,num_equi,xmin,xmax,num_iter,delta,out_iter)
pos=zeros(1,num_iter+2);
E_L=zeros(1,num_iter+2);
E_average=zeros(1,1);
E_variance=zeros(1,1);
for i=1:out_iter
[pos1,E_L1,E_average1,E_variance1]=Quantum_Monte_carlo(num_equi,alpha,xmin,xmax,num_iter,delta);
pos=[pos;pos1];
E_L=[E_L;E_L1];
E_average=[E_average,E_average1];
E_variance=[E_variance,E_variance1];
end
E_average_tot=sum(E_average(1,2:end))/out_iter;
E_variance_tot=sum(E_variance(1,2:end))/out_iter;
end

% analytic solution
function [E_expetation]=cal_ana(h,xmin,xmax,alpha)
x=xmin:h:xmax;
N=size(x);
H=alpha+x.^2*(1/2-2*alpha^2);
H1=ones(1,N(1,2));
wave_function=exp(-alpha*x.^2);
% for i=1:N(1,2)
[integral1]=cal_integral(wave_function,h,H1); % integral of square absolute wave function 
[integral]=cal_integral(wave_function,h,H); % integral of square absolute of wave function times the hamiltonian
E_expetation=integral/integral1; % analytic solution
% end
end

% integral
% h: step size
% H: hamiltonian
function [integral]=cal_integral(wave_funciton,h,H)
integral=0;
N=size(wave_funciton);
    for i=1:N(1,2)-1
    integral=integral+wave_funciton(1,i)*h*H(1,i);
    end
end

function [init_pos]=random_init_pos(xmin,xmax)
init_pos=xmin+(xmax-xmin)*rand();
end

function [pos_update]=displacement(delta,pos)
pos_update=pos+delta*(rand()-0.5);
end

function [pos,E_L1,E_average,E_variance]=Quantum_Monte_carlo(num_equi,alpha,xmin,xmax,num_iter,delta)
pos=zeros(1,1);
[init_pos]=random_init_pos(xmin,xmax);
pos(1,1)=init_pos;
[pos_update]=displacement(delta,init_pos);
psi_old=exp(-alpha*init_pos^2);
psi_old_2=psi_old^2; % probability density, psi^2
psi_new=exp(-alpha*pos_update^2);
psi_new_2=psi_new^2;
p=psi_new_2/psi_old_2;
acc=rand();
E_L=zeros(1,1);E_L(1,1)=alpha+pos(1,1)^2*(1/2-2*alpha^2);
if acc<p
    pos=[pos,pos_update];
else
    pos=[pos,pos(1,1)]; % save data when position does not change
end
for i=1:num_iter
[pos_update]=displacement(delta,pos(1,end));
psi_old=exp(-alpha*pos(1,end)^2);
psi_old_2=psi_old^2; % probability density, psi^2
psi_new=exp(-alpha*pos_update^2);
psi_new_2=psi_new^2;
p=psi_new_2/psi_old_2;
acc=rand();
   if acc<p
       pos=[pos,pos_update];
   else
       pos=[pos,pos(1,end)];
   end
end
E_L1=alpha+pos.^2*(1/2-2*alpha^2);
E_equi=E_L1(1,num_equi:end);
length_E_equi=size(E_equi);
E_average=sum(E_equi)/length_E_equi(1,2);
E_variance=sum((E_equi-E_average).^2)/length_E_equi(1,2);
end

function [xmin,fmin]=golden_search(ax,bx,cx,tol,func,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6)
R=0.61803399;C=1.0-R; % golden ratios
x0=ax;
x3=cx;
if abs(cx-bx)>abs(bx-ax)
    x1=bx;
    x2=bx+C*(cx-bx);
else
    x2=bx;
    x1=bx-C*(bx-ax);
end
f1=func(x1,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6);
f2=func(x2,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6);
while abs(x3-x0)>tol*(abs(x1)+abs(x2))
    if f2<f1
        x0=x1;
        x1=x2;
        x2=R*x2+C*x3;
        f1=f2;f2=func(x2,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6);
    else
        x3=x2;
        x2=x1;
        x1=R*x1+C*x0;
        f2=f1;
        f1=func(x1,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6);
    end
end
if f1<f2
    xmin=x1;
    fmin=f1;
else
    xmin=x2;
    fmin=f2;
end
end

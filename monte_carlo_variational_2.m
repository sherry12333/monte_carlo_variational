clc,clear,close all

% 2D
% parameters
xmin=-5;
xmax=5;
ymin=-5;
ymax=5;
alpha=0.45;
delta=0.2; % range of trial walk
num_iter=2000; % number of iterations inside one VMC
num_equi=500; % after which number begin to calculate the average and variance
out_iter=500; % number of iterations that set initialization to VMC
h=0.1; % discretize on x axis, step size
% lamda=1;
% lamda=[1,2,5,8];% 5,8
lamda=1;
lamda1=8;

ax=0;
bx=0.2;
cx=1;
tol=0.001;

% % plot
alpha2=0.4:0.01:0.6;
% alpha1=0.1:0.1:1;
alpha1=0.05:0.01:1;
alpha=0.05:0.05:0.95;
tic
E_average_tot=zeros(length(lamda),length(alpha));E_variance_tot=zeros(length(lamda),length(alpha));
E_expetation=zeros(length(lamda),length(alpha));
% for j=1:length(lamda)
for i=1:length(alpha1)
[E_average_tot(1,i),E_variance_tot(1,i)]=VMC_2d(alpha1(1,i),lamda,delta,xmin,xmax,ymin,ymax,out_iter,num_iter,num_equi);
% [E_expetation(j,i)]=cal_ana(h,xmin,xmax,alpha(1,i));
% E_expetation=[E_expetation,E_expetation1];
end
subplot(1,2,1)
yyaxis left
% plot(alpha,E_expetation(j,:),'LineWidth',2);hold on;
scatter(alpha1,E_average_tot(1,:),'filled');
xlabel('α');
ylabel('Energy [ℏω]');
yyaxis right
scatter(alpha1,E_variance_tot(1,:),'filled');
ylabel('Variance');
title(['λ = ',num2str(lamda)]);
% 
E_average_tot=zeros(length(lamda),length(alpha1));E_variance_tot=zeros(length(lamda),length(alpha1));
E_expetation=zeros(length(lamda),length(alpha1));
for i=1:length(alpha1)
[E_average_tot(1,i),E_variance_tot(1,i)]=VMC_2d(alpha1(1,i),lamda1,delta,xmin,xmax,ymin,ymax,out_iter,num_iter,num_equi);
end
subplot(1,2,2)
yyaxis left
% plot(alpha,E_expetation(j,:),'LineWidth',2);hold on;
scatter(alpha1,E_average_tot(1,:),'filled');
xlabel('α');
ylabel('Energy [ℏω]');
yyaxis right
scatter(alpha1,E_variance_tot(1,:),'filled');
ylabel('Variance');
title(['λ = ',num2str(lamda1)])

% end
toc

% figure()
% alpha1=0.4:0.01:0.6;
% % alpha1=0.1:0.1:1;
% alpha2=0.05:0.01:1;
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

xmin_gs=zeros(1,length(lamda));
fmin_gs=zeros(1,length(lamda));
for i=1:length(lamda)
% [E_average_tot]=VMC_2d(alpha,lamda,delta,xmin,xmax,ymin,ymax,out_iter,num_iter,num_equi);
[xmin_gs(1,i),fmin_gs(1,i)]=golden_search(ax,bx,cx,tol,@VMC_2d,lamda(1,i),delta,xmin,xmax,ymin,ymax,out_iter,num_iter,num_equi);
end
% E_variance_tot
% E_average_tot: average of energy
% E_variance_tot: variance of energy
% alpha: variation parameter
% lamda: correlation parameter
% xmin: minimum of value on x axis
% out_iter: number of outer iterations
% num_iter: number of inner iterations
% num_equi: number of iterations to reach equilibrium
function [E_variance_tot]=VMC_2d(alpha,lamda,delta,xmin,xmax,ymin,ymax,out_iter,num_iter,num_equi)
E_average=zeros(1,1);E_variance=zeros(1,1);
for j=1:out_iter
pos1=zeros(2,1);pos2=zeros(2,1);
[init_pos1,init_pos2]=random_init_pos(xmin,xmax,ymin,ymax); % initialize positions
pos1(:,end)=init_pos1;pos2(:,end)=init_pos2;
E=zeros(1,1);
r_new=sqrt((pos1(1,end)-pos2(1,end))^2+(pos1(2,end)-pos2(2,end))^2); % calculate distance between the particles
% calculate the second derivatives of x1,x2,y1,y2
second_deri_x1=-1+lamda/(1+alpha*r_new)^2/r_new-lamda*(1+3*alpha*r_new)*(pos1(1,end)-pos2(1,end))^2/(1+alpha*r_new)^3/r_new^3+(-pos1(1,end)+lamda*(pos1(1,end)-pos2(1,end))/(1+alpha*r_new)^2/r_new)^2;
second_deri_x2=-1+lamda/(1+alpha*r_new)^2/r_new-lamda*(1+3*alpha*r_new)*(pos2(1,end)-pos1(1,end))^2/(1+alpha*r_new)^3/r_new^3+(-pos2(1,end)+lamda*(pos2(1,end)-pos1(1,end))/(1+alpha*r_new)^2/r_new)^2;
second_deri_y1=-1+lamda/(1+alpha*r_new)^2/r_new-lamda*(1+3*alpha*r_new)*(pos1(2,end)-pos2(2,end))^2/(1+alpha*r_new)^3/r_new^3+(-pos1(2,end)+lamda*(pos1(2,end)-pos2(2,end))/(1+alpha*r_new)^2/r_new)^2;
second_deri_y2=-1+lamda/(1+alpha*r_new)^2/r_new-lamda*(1+3*alpha*r_new)*(pos2(2,end)-pos1(2,end))^2/(1+alpha*r_new)^3/r_new^3+(-pos2(2,end)+lamda*(pos2(2,end)-pos1(2,end))/(1+alpha*r_new)^2/r_new)^2;
% calculate the energy
E(1,1)=(pos1(1,end)^2+pos2(1,end)^2+pos1(2,end)^2+pos2(2,end)^2)/2-(second_deri_x1+second_deri_x2+second_deri_y1+second_deri_y2)/2+lamda/r_new; 
  for i=1:num_iter
    [pos_update1]=displacement(delta,pos1(:,end));
    [pos_update2]=displacement(delta,pos2(:,end));
    r_old=sqrt((pos1(1,end)-pos2(1,end))^2+(pos1(2,end)-pos2(2,end))^2);
    r_new=sqrt((pos_update1(1,1)-pos_update2(1,1))^2+(pos_update1(2,1)-pos_update2(2,1))^2);
    % calculate the wave function
    psi_old=exp(-(pos1(1,end)^2+pos2(1,end)^2+pos1(2,end)^2+pos2(2,end)^2)/2+(lamda*r_old)/(1+alpha*r_old));
    psi_new=exp(-(pos_update1(1,end)^2+pos_update2(1,end)^2+pos_update1(2,end)^2+pos_update2(2,end)^2)/2+(lamda*r_new)/(1+alpha*r_new));
    psi_old_2=psi_old^2;
    psi_new_2=psi_new^2;
    p=psi_new_2/psi_old_2;
    acc=rand();
    if acc<p
        pos1=[pos1,pos_update1];
        pos2=[pos2,pos_update2];
    else
        pos1=[pos1,pos1(:,end)];
        pos2=[pos2,pos2(:,end)];
    end
    r_new=sqrt((pos1(1,end)-pos2(1,end))^2+(pos1(2,end)-pos2(2,end))^2);
    second_deri_x1=-1+lamda/(1+alpha*r_new)^2/r_new-lamda*(1+3*alpha*r_new)*(pos1(1,end)-pos2(1,end))^2/(1+alpha*r_new)^3/r_new^3+(-pos1(1,end)+lamda*(pos1(1,end)-pos2(1,end))/(1+alpha*r_new)^2/r_new)^2;
    second_deri_x2=-1+lamda/(1+alpha*r_new)^2/r_new-lamda*(1+3*alpha*r_new)*(pos2(1,end)-pos1(1,end))^2/(1+alpha*r_new)^3/r_new^3+(-pos2(1,end)+lamda*(pos2(1,end)-pos1(1,end))/(1+alpha*r_new)^2/r_new)^2;
    second_deri_y1=-1+lamda/(1+alpha*r_new)^2/r_new-lamda*(1+3*alpha*r_new)*(pos1(2,end)-pos2(2,end))^2/(1+alpha*r_new)^3/r_new^3+(-pos1(2,end)+lamda*(pos1(2,end)-pos2(2,end))/(1+alpha*r_new)^2/r_new)^2;
    second_deri_y2=-1+lamda/(1+alpha*r_new)^2/r_new-lamda*(1+3*alpha*r_new)*(pos2(2,end)-pos1(2,end))^2/(1+alpha*r_new)^3/r_new^3+(-pos2(2,end)+lamda*(pos2(2,end)-pos1(2,end))/(1+alpha*r_new)^2/r_new)^2;
    % calculate the energy
    E1=(pos1(1,end)^2+pos2(1,end)^2+pos1(2,end)^2+pos2(2,end)^2)/2-(second_deri_x1+second_deri_x2+second_deri_y1+second_deri_y2)/2+lamda/r_new;
    E=[E,E1];
  end
E_equi=E(1,num_equi:end);
E_average1=sum(E_equi)/length(E_equi);
E_variance1=sum((E_equi-E_average1).^2)/length(E_equi);
E_average=[E_average,E_average1];
E_variance=[E_variance,E_variance1];
end
E_average_tot=sum(E_average)/(length(E_average)-1);
E_variance_tot=sum(E_variance)/(length(E_variance)-1);
end

% move to new position
function [pos_update]=displacement(delta,pos)
pos_update(1,1)=pos(1,1)+delta*(rand()-0.5);
pos_update(2,1)=pos(2,1)+delta*(rand()-0.5);
end

% initialize the positions for 2 particles
function [init_pos1,init_pos2]=random_init_pos(xmin,xmax,ymin,ymax)
init_pos1=zeros(1,2);init_pos2=zeros(1,2);
init_pos1(1,1)=xmin+(xmax-xmin)*rand();
init_pos1(1,2)=ymin+(ymax-ymin)*rand();
init_pos2(1,1)=xmin+(xmax-xmin)*rand();
init_pos2(1,2)=ymin+(ymax-ymin)*rand();
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
    integral=integral+wave_funciton(1,i)^2*h*H(1,i);
    end
end

function [xmin,fmin]=golden_search(ax,bx,cx,tol,func,lamda,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6,parameter7,parameter8)
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
f1=func(x1,lamda,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6,parameter7,parameter8);
f2=func(x2,lamda,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6,parameter7,parameter8);
while abs(x3-x0)>tol*(abs(x1)+abs(x2))
    if f2<f1
        x0=x1;
        x1=x2;
        x2=R*x2+C*x3;
        f1=f2;f2=func(x2,lamda,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6,parameter7,parameter8);
    else
        x3=x2;
        x2=x1;
        x1=R*x1+C*x0;
        f2=f1;
        f1=func(x1,lamda,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6,parameter7,parameter8);
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

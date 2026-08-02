% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% par case
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% loglog plot
clear all; clc;
B = 1;
N = 2;

% variable range
t = [3.0:0.001:3.13];

% function definition
f = (3./(4*N*(N+2))) .* ( 1./t.^2 + 2*B^2./sin(B*t).^2 );
f_lower = (1+2*B*t./abs(sin(B*t))).^2./(4*N*(N+2)*t.^2);


v5_2to3 = load("upper_t3to3dot13.mat");
lower = load("lower_t3to3dot13.mat");
lower = lower.lower_res;

specialstate=[8.38285,9.69421,11.3423,13.4531,16.2174,19.9359,25.1039,32.586,44.0018,62.6829,96.4105,167.053,357.536,1240.24];

% plot
figure('Position',[1, 1, 604, 423]);
% plot(t, f, 'LineWidth', 2); hold on;
loglog(t, f_lower, 'LineWidth',2,'Color','k'); hold on;
loglog([3.0:0.01:3.13],specialstate,'diamond','LineStyle','none','MarkerSize', 8,'Color','green'); hold on;

loglog([3.0:0.01:3.13], ...
     [v5_2to3.upper_res], ...
     'o',  ...
     'MarkerSize',8,'Color','red'); hold on;

loglog([3.0:0.01:3.13],lower,'+','LineStyle','none','MarkerSize', 8,'Color','green'); hold on;


xlim([3.0 3.13]);

% scatter([1.5,1.8,2.1,2.4,2.7,3.0],v7.upper_res(2,5:10),'filled','diamond');
% scatter([1.5,1.8,2.1,2.4,2.7,3.0],v9.upper_res(2,5:10),'filled','pentagram');

lgd = legend('Analytical lower bound','Heuristic state','$B^{\left( i \right)}_+$','$B^{\left( i \right)}_-$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
grid on;
ax = gca;                 
ax.GridLineStyle = '--';  

%%
clear all; clc;
B = 1;
N = 2;

% variable range
t = [3.0:0.01:3.13];

% function definition
f = (3./(4*N*(N+2))) .* ( 1./t.^2 + 2*B^2./sin(B*t).^2 );
% f_lower = (1+2*B*t./abs(sin(B*t))).^2./(4*N*(N+2)*t.^2);


v5_2to3 = load("upper_t3to3dot13.mat");

specialstate=[8.38285,9.69421,11.3423,13.4531,16.2174,19.9359,25.1039,32.586,44.0018,62.6829,96.4105,167.053,357.536,1240.24];

% plot
figure('Position',[1, 1, 604, 423]);
loglog(t, f, 'LineWidth', 2); hold on;
loglog([3.0:0.01:3.13],specialstate,'diamond','LineStyle','none','MarkerSize', 8,'Color','green'); hold on;

loglog([3.0:0.01:3.13], ...
     [v5_2to3.upper_res(1,:)], ...
     'o',  ...
     'MarkerSize',8,'Color','red'); hold on;

xlim([3.0 3.13]);

% scatter([1.5,1.8,2.1,2.4,2.7,3.0],v7.upper_res(2,5:10),'filled','diamond');
% scatter([1.5,1.8,2.1,2.4,2.7,3.0],v9.upper_res(2,5:10),'filled','pentagram');

lgd = legend('Heuristic state in Ref. [22]','Heuristic state in Ref. [42]','$B^{\left( i \right)}_+$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
grid on;
ax = gca;                
ax.GridLineStyle = '--';

%%
clear all; clc;
len = 20;
xaxis = linspace(0,1,len+1); yaxis = linspace(1,0,len+1);
[X, Y] = meshgrid(xaxis,yaxis);
data = load("colormap_data_20.mat");

figure;
data.data(data.data == 0) = NaN;
t = 3; N = 2;
analytical = (1+2*t/abs(sin(t)))^2/(4*N*(N+2)*t^2);
data.data = data.data - analytical;
h = imagesc(xaxis, yaxis, data.data);

set(h, 'AlphaData', ~isnan(data.data));
set(gca, 'YDir', 'normal');
set(gca, 'XTick', linspace(0,1,11));
set(gca, 'YTick', linspace(0,1,11));


axis equal;
axis([0 1 0 1]); 
colorbar;
xlabel('$\theta_1$',"Interpreter","latex"); ylabel('$\theta_2$',"Interpreter","latex");
%%
% open figure
clear all; clc;
fig = openfig('path','reuse');
set(fig,'Units','pixels');
set(fig,'Position',[1 1 604 423]);
inspect(fig);
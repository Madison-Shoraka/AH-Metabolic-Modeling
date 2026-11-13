weight=0.055115566;%input("Patient weight in pounds: ");
gender='M';%input("Patient gender (M/F): ",'s');
%Widmark factor
switch gender
    case 'F'
        rho = 0.55; 
    case 'M'
        rho=0.68;
end 

episodes=input("Number of drinking episodes: ");
drink=zeros(episodes,1);
time_between=zeros(episodes-1,1);

for i = 1:episodes
    drink(i) = input("Number of drinks consumed for episode " + string(i) + ":");
    if i>1
        time_between(i-1) = input("Time between episode " + string(i-1) + " and episode " + string(i) + " (in min):");
    end
end

perturb=zeros(episodes,1);
perturb(1)=100;
for j = 2:episodes
    perturb(j) = perturb(j-1)+time_between(j-1);
end


mM_drink= ((drink*14)./((weight*453.592)*rho))*100*217; %widmark2009erik, ~14g/drink, 453.592 g conversion from lbs, 217: BAC=>mM

health_state = input("AST levels (U/L):");
base_line = input("AST normal lab range (dependends on lab):");
health_state=health_state/base_line;
%healthy<=1
%mild: 1-3
%moderate: 3-10
%severe>=10
 %interpolated values from rate_estimations.m
GSH_vel_func = [1.9918    2.0117    2.0316    2.0515    2.0714    2.0858    2.1003    2.1147    2.1291    2.1225    2.1158    2.1091    2.1024];
Met_vel_func = [0.4154    0.4129    0.4104    0.4079    0.4054    0.4067    0.4080    0.4094    0.4107    0.4004    0.3900    0.3797    0.3694];
Urea_vel_func = [0.7047    0.6836    0.6626    0.6415    0.6205    0.5991    0.5776    0.5562    0.5348    0.5554    0.5761    0.5967    0.6173];
FABs_vel_func = [0.8881    0.8918    0.8955    0.8991    0.9028    0.9084    0.9141    0.9197    0.9254    0.9217    0.9181    0.9145    0.9109];
general_cubic_poly=@(x,p1,p2,p3,p4) p1*x^3+p2*x^2+p3*x+p4;
if health_state>=10
    decay=.005; %capturing permeability of the cell/ ability to transport materials accross cell boundary when necessary
    y1=.5; %yeild coefficient: oxidative stress to be detoxed by GSH per mM of alcohol  proportional to CYP2E1 activation
    y2=.25; %yeild coefficient: mM acetyl-CoA produced per mM of alcohol (14g) proportional to acetate shuttle deactivation
    y4=1;
    velocities_IC= [GSH_vel_func(13) Met_vel_func(13) Urea_vel_func(13) FABs_vel_func(13)];
elseif health_state<=1
    decay=0.00; %capturing permeability of the cell/ ability to transport materials accross cell boundary when necessary
    y1=.6; %yeild coefficient: oxidative stress to be detoxed by GSH per mM of alcohol proportional to CYP2E1 activation
    y2=.4; %yeild coefficient: mM acetyl-CoA produced per mM of alcohol proportional to acetate shuttle deactivation
    y4=1;
    velocities_IC= [GSH_vel_func(1) Met_vel_func(1) Urea_vel_func(1) FABs_vel_func(1)];
else
    %linearly interpolated based on two above states at health_state=6 and health_state=40
    decay=-0.0009*health_state+0.0109; 
    y1=general_cubic_poly(health_state,0.0015,-0.0333,0.304,-0.175);
    y2=general_cubic_poly(health_state,0.0037,-0.0604,0.4246,0.1321);
    y4=1;
    velocities_IC=[general_cubic_poly(health_state,-0.0002,0.0012,0.025,1.9659) general_cubic_poly(health_state, -0.0002,0.003,-0.0126,0.4253) general_cubic_poly(health_state,0.0006,-0.0061,-0.0113,0.7215) general_cubic_poly(health_state, -0.0002,0.002,-0.0011 ,0.8874)];

end


   





T=max(perturb)+4000;%+2880;

 %% BASE PARAMETERS

    %standard parameters
tau=1000000;
y3=1; %Time delay (min) between alcohol metabolism and acetyl-CoA production
dt=.01;
time=0:dt:T-dt;
sigma=4;
time_to_detox=60;
velocities_IC= [2.0737    0.4002    0.6193    0.9068]; %GSH,MET,Urea,FABs
l=[1 0 0 0 %metabolite = GSH_m
   0 1 0 0 %metbolite = Met_c
   0 0 1 0 %metabolite = ammonia_m
   0 0 1 1 %metabolite = a-CoA_m
   0 0 0 1]; %metabolite=malonyl-CoA_m
type=[1 1 0 0 0]; %1=produced, 0=used
concentration_IC=[7 .08 0 0 0];
critical_concentration_vec=[5.25,.04,1,1,.3];
critical_concentrations=[1 0 0 0 %metabolite = GSH_m
    0 1 0 0 %metbolite = MET_c
    0 0 1 0 %metabolite = ammonia_m
    0 0 1 1 %metabolite = a-CoA_m
    0 0 0 1]; %metabolite = malonyl-CoA_m 
critical_concentrations=critical_concentration_vec'.*critical_concentrations;
max_concentration=[concentration_IC(1) concentration_IC(2) 100 100 100];
beta = 7; % sharpness parameter
alcohol_met=-.05;
Adj_mat=[0     7     4     2
     5     0     2     1
     6     5     0     4
     2     1     3     0];

response_rate=.01;
    


% Run model
[concentration,velocities,alcohol,H2O2]=Patient_specific_function(y4,tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
                 
%plot parameters
figure
hold on
plot(time, concentration(:,1),'LineWidth',3,'Color',[1, 0, 0])
plot(time, concentration(:,2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(time, concentration(:,3),'LineWidth',3,'Color',[0,0,1])
plot(time, concentration(:,4),'LineWidth',3,'Color',[1,1,0])
plot(time, concentration(:,5),'LineWidth',3,'Color',[0,1,1])
xlabel('Time (min)');
ylabel('Shifted and Scaled Concentration');
title('Concentration Over Time');
legend(["[GSH]_m","[MET]_c","[ammonia]_m","[a-CoA]_m","[malonyl-CoA]_m"])

%plot parameters
figure
hold on
plot(time,H2O2(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
xlabel('Time (min)');
ylabel('Shifted Concentration \mu M');
title('[H_2O_2]_m Over Time');


grid on;
hold off
figure
hold on
plot(time, velocities(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(time, velocities(:,2)./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(time,velocities(:,3)./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(time, velocities(:,4)./velocities_IC(4),'LineWidth',3,'Color',[1, 1, 0])
xlabel('Time (min)');
ylabel('velocities');
title('velocities Over Time');
legend(["GSH","MET","Urea","FABs"])
grid on;
hold off

grid on;
hold off
figure
hold on
plot(time, alcohol,'LineWidth',3,'Color',[1, 0, 0])
xlabel('Time (min)');
ylabel('alcohol(mM)');
title('alcohol Over Time');
grid on;
hold off

"GSH net velocity change"
sum(velocities(:,1)-velocities_IC(1))

"MET net velocity change"
dt*sum(velocities(:,2)-velocities_IC(2))

"Urea net velocity change"
dt*sum(velocities(:,3)-velocities_IC(3))

"FABs net velocity change"
dt*sum(velocities(:,4)-velocities_IC(4))


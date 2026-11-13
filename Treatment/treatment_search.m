load("param_vals.mat")
index_failure_case1 = find(met_failure_NO_oxidative_failure1+oxidative_failure_AND_met_failure1);
%% standard parameters
y3=1; %Time delay (min) between alcohol metabolism and acetyl-CoA production
tau=10000;
dt=.005;
T=25000;
time=0:dt:T-dt;
nt=length(time);
sigma=5;
time_to_detox=60;
response_rate=.01;
l=[1 0 0 0 %metabolite = GSH_m
   0 1 0 0 %metbolite = SAM_c
   0 0 1 0 %metabolite = ammonia_m
   0 0 1 1 %metabolite = a-CoA_m
   0 0 0 1]; %metabolite=malonyl-CoA_m
type=[1 1 0 0 0]; %1=produced, 0=used
beta = 7; % sharpness parameter
alcohol_met=-.05;
Adj_mat=[0     7     4     2
     5     0     2     1
     6     5     0     4
     2     1     3     0];
velocities_IC= [2.0737    0.4002    0.6193    0.9068]; %GSH,MET,Urea,FABs


concentration_IC=[10 .08 0 0 0];
critical_concentration_vec=[7.5,.04,1,2,.3];
critical_concentrations=[1 0 0 0 %metabolite = GSH_m
    0 1 0 0 %metbolite = SAM_c
    0 0 1 0 %metabolite = ammonia_m
    0 0 1/2 1/2 %metabolite = a-CoA_m
    0 0 0 1]; %metabolite = malonyl-CoA_m 
critical_concentrations=critical_concentration_vec'.*critical_concentrations;
max_concentration=[concentration_IC(1) concentration_IC(2) 100 100 100];

mM_drink=17.4; %mM alcohol per drink
perturb=[200];

%% Set up grid
decay_min=.001;
decay_max=.02;
y1_min=.01;
y1_max=.75;
y2_min=.1;
y2_max=1;

dimension=3;
m=10;
resolution=(m+1)^dimension;


decay_vec=decay_min:(decay_max-decay_min)/m:decay_max;
y1_vec=y1_min:(y1_max-y1_min)/m:y1_max;
y2_vec=y2_min:(y2_max-y2_min)/m:y2_max;

[decay_list,y1_list,y2_list]=meshgrid(decay_vec,y1_vec,y2_vec);
decay_list=reshape(decay_list,resolution,1);
y1_list=reshape(y1_list,resolution,1);
y2_list=reshape(y2_list,resolution,1);

epsilon_t=200;
epsilon_t=T/dt-epsilon_t/dt;

len=length(index_failure_case1);

total_concentration_NOtreatement=zeros(len,5);
total_velocities_NOtreatement=zeros(len,4);
total_H2O2_NOtreatment=zeros(len,1);
oxidative_stress_NOtreatment=zeros(len,1);
metabolic_decomp_NOtreatment=zeros(len,4);

total_concentration_1treatement=zeros(len,5);
total_velocities_1treatement=zeros(len,4);
total_H2O2_1treatment=zeros(len,1);
oxidative_stress_1treatment=zeros(len,1);
metabolic_decomp_1treatment=zeros(len,4);

total_concentration_2treatement=zeros(len,5);
total_velocities_2treatement=zeros(len,4);
total_H2O2_2treatment=zeros(len,1);
oxidative_stress_2treatment=zeros(len,1);
metabolic_decomp_2treatment=zeros(len,4);

total_concentration_3treatement=zeros(len,5);
total_velocities_3treatement=zeros(len,4);
total_H2O2_3treatment=zeros(len,1);
oxidative_stress_3treatment=zeros(len,1);
metabolic_decomp_3treatment=zeros(len,4);



parfor i=1:len
    j=index_failure_case1(i);
    decay=decay_list(j)
    y1=y1_list(j)
    y2=y2_list(j) 
    %% no treatment
    treatment_time=-1;
    [concentration_normal,velocities_normal,alcohol,H2O2_normal,effect,total_oxidative_stress_normal]=treatment_cortocosteriods(treatment_time,tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
         
    total_concentration_NOtreatement(i,:)=dt*sum(concentration_normal)-T*concentration_IC;
    total_velocities_NOtreatement(i,:)=dt*sum(velocities_normal)-T*velocities_IC;
    total_H2O2_NOtreatment(i)=dt*sum(H2O2_normal)-0;
    oxidative_stress_NOtreatment(i)=dt*sum(total_oxidative_stress_normal);
    metabolic_decomp_NOtreatment(i,:)=dt*sum(abs(velocities_normal(epsilon_t:end,:)-velocities_IC));
    
    
    %% treatment 1
    time_IC=2880;
    number_treatments=7;
    treatment_time1=time_IC:1440:number_treatments*1440+time_IC;
    [concentration_treated1,velocities_treated1,alcohol,H2O2_treated1,effect1,total_oxidative_stress_treatment1]=treatment_cortocosteriods(treatment_time1,tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
     
    total_concentration_1treatement(i,:)=dt*sum(concentration_treated1)-T*concentration_IC;
    total_velocities_1treatement(i,:)=dt*sum(velocities_treated1)-T*velocities_IC;
    total_H2O2_1treatment(i)=dt*sum(H2O2_treated1)-0;
    oxidative_stress_1treatment(i)=dt*sum(total_oxidative_stress_treatment1);
    metabolic_decomp_1treatment(i,:)=dt*sum(abs(velocities_treated1(epsilon_t:end,:)-velocities_IC));
    
    
    %% treatment 2
    time_IC=4320;
    number_treatments=7;
    treatment_time2=time_IC:1440:number_treatments*1440+time_IC;
    [concentration_treated2,velocities_treated2,alcohol,H2O2_treated2,effect2,total_oxidative_stress_treatment2]=treatment_cortocosteriods(treatment_time2,tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
     
    total_concentration_2treatement(i,:)=dt*sum(concentration_treated2)-T*concentration_IC;
    total_velocities_2treatement(i,:)=dt*sum(velocities_treated2)-T*velocities_IC;
    total_H2O2_2treatment(i)=dt*sum(H2O2_treated2)-0;
    oxidative_stress_2treatment(i)=dt*sum(total_oxidative_stress_treatment2);
    metabolic_decomp_2treatment(i,:)=dt*sum(abs(velocities_treated2(epsilon_t:end,:)-velocities_IC));
    
    
    %% treatment 3
    time_IC=5760;
    number_treatments=7;
    treatment_time3=time_IC:1440:number_treatments*1440+time_IC;
    [concentration_treated3,velocities_treated3,alcohol,H2O2_treated3,effect3,total_oxidative_stress_treatment3]=treatment_cortocosteriods(treatment_time3,tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
     
    total_concentration_3treatement(i,:)=dt*sum(concentration_treated3)-T*concentration_IC;
    total_velocities_3treatement(i,:)=dt*sum(velocities_treated3)-T*velocities_IC;
    total_H2O2_3treatment(i)=dt*sum(H2O2_treated3)-0;
    oxidative_stress_3treatment(i)=dt*sum(total_oxidative_stress_treatment3);
    metabolic_decomp_3treatment(i,:)=dt*sum(abs(velocities_treated3(epsilon_t:end,:)-velocities_IC));

end 

%% Data processing
met_decomp_0treatment=max(metabolic_decomp_NOtreatment')';
met_decomp_1treatment=max(metabolic_decomp_1treatment')';
met_decomp_2treatment=max(metabolic_decomp_2treatment')';
met_decomp_3treatment=max(metabolic_decomp_3treatment')';

total_decomp=[met_decomp_0treatment,met_decomp_1treatment,met_decomp_2treatment,met_decomp_3treatment];
total_oxidative_stress=[oxidative_stress_NOtreatment,oxidative_stress_1treatment,oxidative_stress_2treatment,oxidative_stress_3treatment];
treatment_list=["None","1","2","3"];
types=["MET AND OX", "~MET AND OX", "MET AND ~OX"];
sorted_treatment_met_decomp=zeros(len,4);
defined_types=strings(len,1);

for i=1:len
    j=index_failure_case1(i);
    if oxidative_failure_NO_met_failure1(j)==1
        defined_types(i)=types(2);
    elseif oxidative_failure_AND_met_failure1(j)==1
        defined_types(i)=types(1);
    else
        defined_types(i)=types(3);
    end 
    [v_sorted, sorted_treatment_met_decomp(i,:)] = sort(total_decomp(i,:));
    
end 



total=[sorted_treatment_met_decomp,defined_types];
sum_oxidative_failure_AND_met_failure1=sum(oxidative_failure_AND_met_failure1);
sum_oxidative_failure_NO_met_failure1=sum(oxidative_failure_NO_met_failure1);
sum_met_failure_NO_oxidative_failure1=sum(met_failure_NO_oxidative_failure1);
% categorize

indicies_treatment_at_least_one_success=[];
indicies_treatment_allfailure=[];
indicies_treatment_at_least_one_failure=[];
indicies_treatment_worse=[];


for i=1:len
    if met_decomp_0treatment(i)>.1
        if met_decomp_1treatment(i) <=.1 || met_decomp_2treatment(i) <=.1|| met_decomp_3treatment(i) <=.1
            indicies_treatment_at_least_one_success=[indicies_treatment_at_least_one_success,i];
        else
            indicies_treatment_allfailure=[indicies_treatment_allfailure,i];
        end 
        if met_decomp_1treatment(i) >.1 || met_decomp_2treatment(i) >.1|| met_decomp_3treatment(i) >.1
        indicies_treatment_at_least_one_failure=[indicies_treatment_at_least_one_failure,i];
        end 
    end 
    if met_decomp_0treatment(i)<=.1 && (met_decomp_1treatment(i) >.1 || met_decomp_2treatment(i) >.1|| met_decomp_3treatment(i) >.1)
        indicies_treatment_worse=[indicies_treatment_worse,i];
    end 
    
end


param_that_stab_eventually_no_treatment=(met_decomp_0treatment<=.1);
param_unstab=find((met_decomp_0treatment>.1));
length_stable=sum(param_that_stab_eventually_no_treatment);
length_no_stable=length(met_decomp_0treatment)-length_stable;
total_unstable_without_treatment=total(find(1-param_that_stab_eventually_no_treatment),:);

treatment_improvement=zeros(length_no_stable,3);
treatment_not_responsive=zeros(length_no_stable,3);
for i=1:length_no_stable
    sorted=double(total_unstable_without_treatment(i,1:4));
    idx1 = find(sorted == 1, 1);
    idx2 = find(sorted == 2, 1);
    idx3 = find(sorted == 3, 1);
    idx4 = find(sorted == 4, 1);
    if idx1>idx2
        treatment_improvement(i,1)=1;
    else 
        treatment_not_responsive(i,1)=1;
    end 
    if idx1>idx3
        treatment_improvement(i,2)=1;
    else 
        treatment_not_responsive(i,2)=1;
    end
    if idx1>idx4
        treatment_improvement(i,3)=1;
    else 
        treatment_not_responsive(i,3)=1;
    end
end




everything=[y1_list(index_failure_case1(param_unstab)),y2_list(index_failure_case1(param_unstab)),decay_list(index_failure_case1(param_unstab)),total_decomp(param_unstab),total_H2O2_NOtreatment(param_unstab),total_H2O2_1treatment(param_unstab),total_H2O2_2treatment(param_unstab),total_H2O2_3treatment(param_unstab)];
%treatment_stats
disp('Percent of cases that will eventually regain metabolic stability (~ 17 days) ')
disp(length_stable/length(met_decomp_0treatment))

disp('Percent of cases that will not regain metabolic stability (~ 17 days)')
disp(length_no_stable/length(met_decomp_0treatment))

disp('Percent unstable cases that that respond to treatment 1')
disp(sum(treatment_improvement(:,1))/length_no_stable)
disp('Percent unstable cases that that respond to treatment 2')
disp(sum(treatment_improvement(:,2))/length_no_stable)
disp('Percent unstable cases that respond to treatment 3')
disp(sum(treatment_improvement(:,3))/length_no_stable)
disp('Average percent of treatment responsiveness')
disp((.8276+.7241+.8103)/3)



disp('Percent unstable cases that that regain stability treatment 1')
disp(sum(met_decomp_1treatment(param_unstab)<=.1)/length_no_stable)
disp('Percent unstable cases that that regain stability treatment 2')
disp(sum(met_decomp_2treatment(param_unstab)<=.1)/length_no_stable)
disp('Percent unstable cases that regain stability treatment 3')
disp(sum(met_decomp_3treatment(param_unstab)<=.1)/length_no_stable)



disp('Of the cases that will regain stability...')
disp('Percent that treatment would stop the regain of stability')
disp(sum(ismember(find(param_that_stab_eventually_no_treatment),indicies_treatment_worse))/length_stable)





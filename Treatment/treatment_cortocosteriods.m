function [concentration,velocities,alcohol,ROS,effect,total_oxidative_stress]=treatment_cortocosteriods(treatment_time,tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb)
                 
%metabolic model of alcohol consumption
%Cycles included: (1) GSH (2) MET (3) Urea (4)FABs
%Mets tracked: (1) GSH_m (2) SAMe_c (3) ammonia (4) acetyl-CoA (5)manoyl-CoA
%(6) H2O2
           
steriod=0;

nt=T/dt;
effect=zeros(nt,1);
total_oxidative_stress=zeros(nt,1);

num_cycles=length(Adj_mat);

capacity=zeros(nt,num_cycles);
capacity(1,:)=ones(1,num_cycles);
num_mets=5; 
concentration=zeros(nt,num_mets); 
concentration(1,:)=concentration_IC;

velocities=zeros(nt,num_cycles);
velocities(1,:)=velocities_IC;
k_p=[sum((l(1,:).*velocities(1,:))), sum((l(2,:).*velocities(1,:))),sum((l(3,:).*velocities(1,:))),sum((l(4,:).*velocities(1,:))),sum((l(5,:).*velocities(1,:)))];




ROS=zeros(nt,1);
alcohol=zeros(nt,1);
acetate=zeros(nt,1);
perturb=perturb./dt;
treatment_time=treatment_time./dt;
total_cap=velocities(1,:);
    Tv_per_metabolite=[sum((l(1,:).*total_cap)), sum((l(2,:).*total_cap)),sum((l(3,:).*total_cap)),sum((l(4,:).*total_cap)),sum((l(5,:).*total_cap))];
    capacity=[(l(1,:).*total_cap)./Tv_per_metabolite(1); %capacity for cycle i to metabolize met j
    (l(2,:).*total_cap)./Tv_per_metabolite(2)
    (l(3,:).*total_cap)./Tv_per_metabolite(3)
    (l(4,:).*total_cap)./Tv_per_metabolite(4)
    (l(5,:).*total_cap)./Tv_per_metabolite(5)];

for j=2:nt
    
    redistribution=((1-type)'.*concentration(j-1,:)'.*capacity)+type'.*(max_concentration-concentration(j-1,:))'.*capacity;
    z = max((1-type)'.*((concentration(j-1,:)'.*capacity)./critical_concentrations)+...
        type'.*((min(sigma,(critical_concentrations./(concentration(j-1,:)'.*capacity))).*((concentration(j-1,:)'.*capacity)< critical_concentrations))...
        +((max_concentration'-concentration(j-1,:)'.*capacity)./(max_concentration'-critical_concentrations)).*((concentration(j-1,:)'.*capacity)> critical_concentrations).*((concentration(j-1,:)'.*capacity)<= max_concentration')));
    demand=max(redistribution./l./time_to_detox+abs(k_p)'.*capacity./l);
   
    w = exp(beta*(z - max(z))); % stable softmax
    w = w./sum(w); %weight
    b = velocities(j-1,:);  % start allotment b_i
    available_percent = min(...
    (1-type)'.*(1-(concentration(j-1,:)'.*capacity)./critical_concentrations).*((concentration(j-1,:)'.*capacity)<critical_concentrations)+...
    type'.*(1.*((concentration(j-1,:)'.*capacity)>max_concentration')+((((concentration(j-1,:)'.*capacity)-critical_concentrations)./(max_concentration'-critical_concentrations)).*((concentration(j-1,:)'.*capacity)<=max_concentration').*((concentration(j-1,:)'.*capacity)>critical_concentrations))));
    s = available_percent .* b;    % shareable pool s_i
    
    delta = max(0, demand - b);     % deficit delta_i
    R = w .* delta;    % priority-weighted request intensity R_i
    
    rowSums = sum(Adj_mat,2);      
    P = Adj_mat ./ rowSums;          

    % request matrix req_{i,j} = R_i * p_{i,j}
    REQ = (R.') .* P;                          % each row i scaled by R_i

    % donor-normalized allocation from donor j to receiver i:
    % X(i,j) = s_j * REQ(i,j) / sum_i REQ(i,j)
    eps0 = 1e-12; %small coeff to not divide by 0
    colSums = sum(REQ,1) + eps0;                                
    X = (REQ ./ colSums) .* s;       % multiply each column j by s_j

    %  net change in b due to flows
    inflow  = sum(X,2).';     % receiver sums over donors
    outflow = sum(X,1);      % donor sums to receivers
    allotted = b + inflow - outflow;    % updated allotment b_i
    % --- end continuous redistribution ---
    
     %alcohol update rule
    alcohol(j) = max(0,alcohol(j-1) + dt*(alcohol_met)); % Update alcohol consumption (-.017 drink/min)
    
     if alcohol(j-1)>0
        net_diff_alc=dt*(alcohol_met);
    else
        net_diff_alc=0;
    end
    
    %velocities update rule
    velocities(j,:) = velocities(j-1,:)+dt*response_rate.*(allotted-velocities(j-1,:)); 



    %concentrations update rules
    tv_vel=sum((-velocities(j-1,:).*l)')+k_p-decay;
    % GSH
    concentration(j,1) = max(0,min(max_concentration(1),concentration(j-1,1) - dt*tv_vel(1)));
    ROS(j)=ROS(j-1)-y1*net_diff_alc;
    net_diff_ROS=min(ROS(j),dt*tau*concentration(j,1)^2*ROS(j));
    ROS(j)=ROS(j)-net_diff_ROS;
    concentration(j,1)=concentration(j,1)-2*net_diff_ROS;
    total_oxidative_stress(j)=ROS(j)+max(0,1/2*(critical_concentrations(1)-concentration(j,1)));
    %SAM update rule 
    concentration(j,2) = max(0,min(max_concentration(2),concentration(j-1,2) - dt*tv_vel(2))); %SAM
    %ammonia update rule
    concentration(j,3) = max(0,concentration(j-1,3) + dt*tv_vel(3)); %ammonia
    %acetate, A-CoA update rule
    acetate_step=acetate(j-1)-net_diff_alc*y2;
    acetate(j)=max(0,acetate_step-dt*y3);
    net_diff_acetate=acetate_step-acetate(j);
    concentration(j,4) = max(0,concentration(j-1,4) + dt*tv_vel(4)); 
    concentration(j,4)=concentration(j,4)+net_diff_acetate;
    %malonyl-CoA update rule
    concentration(j,5) = max(0,concentration(j-1,5) + dt*tv_vel(5)); %fats
    
    
    steriod=steriod+dt*(-.012)*(steriod>=.012);
    effect(j)=effect(j-1)+dt*(.012)*(steriod>=.012)+dt*(-.0005)*(effect(j-1)>=.0005);
    k_p(1)=velocities_IC(1)-min((.001)*effect(j),.001);

    if ismember(j,treatment_time)==1
        steriod=steriod+1;
    end
    

    if ismember(j,perturb)==1
        alcohol(j)=alcohol(j)+mM_drink;  
    end 

   


   

end



end 

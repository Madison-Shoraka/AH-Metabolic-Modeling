for examples=1:3

    % standard parameters
    y3=1; %Time delay (min) between alcohol metabolism and acetyl-CoA production
    tau=10000;
    dt=.005;
    T=1640+1440;
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
    if examples==1
        Adj_mat=[0     7     4     0
         5     0     0     0
         6     0     0     4
         0     0     3     0];
    elseif examples==2
        Adj_mat=[0     1     0     0
         1     0     0     0
         0     0     0     1
         0     0     1     0];
    elseif examples==3
        Adj_mat=[0     0     0     1
         0     0     1     0
         1     0     0     1
         0     1     0     0];
    end 
    
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
    
    % Set up grid
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
    
    
    % Run cases
    
    %case 1
    total_concentration_case1=zeros(resolution,5);
    total_velocities_case1=zeros(resolution,4);
    total_H2O2_case1=zeros(resolution,1);
    total_GSH_case1=zeros(resolution,1);
    
    
    mM_drink=17.4; %mM alcohol per drink
    perturb=[200]; %minutes from start, can include multiple
    parfor j=1:resolution
            decay=decay_list(j)
            y1=y1_list(j)
            y2=y2_list(j)
            
           [concentration,velocities,alcohol,H2O2]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
            total_concentration_case1(j,:)=dt*sum(concentration)-T*concentration_IC;
            total_velocities_case1(j,:)=dt*sum(velocities)-T*velocities_IC;
            total_H2O2_case1(j)=dt*sum(H2O2)-0;
            
    end 
       
    
    %case 2
    total_concentration_case2=zeros(resolution,5);
    total_velocities_case2=zeros(resolution,4);
    total_H2O2_case2=zeros(resolution,1);
    total_GSH_case2=zeros(resolution,1);
    
    
    mM_drink=17.4; %mM alcohol per drink
    perturb=[200,1640]; %minutes from start, can include multiple (24 hours later)
    parfor j=1:resolution
            decay=decay_list(j)
            y1=y1_list(j)
            y2=y2_list(j)
           [concentration,velocities,alcohol,H2O2]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
            total_concentration_case2(j,:)=dt*sum(concentration)-T*concentration_IC;
            total_velocities_case2(j,:)=dt*sum(velocities)-T*velocities_IC;
            total_H2O2_case2(j)=dt*sum(H2O2)-0;
    end 
    save(sprintf('results_example%d.mat',examples), 'total_concentration_case1', 'total_concentration_case2','total_H2O2_case1','total_H2O2_case2',"total_velocities_case1","total_velocities_case2","total_GSH_case1","total_GSH_case2")

end


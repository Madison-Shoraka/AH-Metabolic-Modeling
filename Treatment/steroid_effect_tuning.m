time_IC=2880;
dose_schedule=1440; % minutes = 24 hours
number_treatments=1;
treatment_time=time_IC:dose_schedule:(number_treatments-1)*dose_schedule+time_IC;
T=15000;
dt=.1;
nt=T/dt;
time=1:nt:T;
steroid=zeros(nt,1);
effect=zeros(nt,1);
alpha1=.012;
alpha2=.0005;
for i=2:nt
    if ismember(i,treatment_time/dt)==1
        steroid=steroid+1;
    end
    steroid(i)=steroid(i-1)+dt*(-alpha1)*(steroid(i-1)>=alpha1);
    effect(i)=effect(i-1)+dt*(alpha1)*(steroid(i-1)>=alpha1)+dt*(-alpha2)*(effect(i-1)>=alpha2);
end 

plot(effect)
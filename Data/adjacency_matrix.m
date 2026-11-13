GSH_L2=readtable('Cycles equations.xlsx',Sheet='GSH_L2');
UREA_L2=readtable('Cycles equations.xlsx',Sheet='Urea_L2');
MET_L2=readtable('Cycles equations.xlsx',Sheet='Met_L2');
FABs_L2=readtable('Cycles equations.xlsx',Sheet='FABs_L2');


GSH_L1_reactants=strrep(rmmissing(GSH_L2.ReactantFromL1Produced, 'MinNumMissing', 1), ' ', '');
GSH_L1_transport= strrep(strrep(strrep(GSH_L1_reactants, ' ', ''),'[c]',''),'[m]','');
GSH_L1_transport=[GSH_L1_transport+"[m]"
    GSH_L1_transport+"[c]"];
UREA_L1_reactants=strrep(rmmissing(UREA_L2.ReactantFromL1Produced, 'MinNumMissing', 1), ' ', '');
UREA_L1_transport= strrep(strrep(strrep(UREA_L1_reactants, ' ', ''),'[c]',''),'[m]','');
UREA_L1_transport=[UREA_L1_transport+"[m]"
    UREA_L1_transport+"[c]"];
MET_L1_reactants=strrep(rmmissing(MET_L2.ReactantFromL1Produced, 'MinNumMissing', 1), ' ', '');
MET_L1_transport= strrep(strrep(strrep(MET_L1_reactants, ' ', ''),'[c]',''),'[m]','');
MET_L1_transport=[MET_L1_transport+"[m]"
    MET_L1_transport+"[c]"];
FABs_L1_reactants=strrep(rmmissing(FABs_L2.ReactantFromL1Produced, 'MinNumMissing', 1), ' ', '');
FABs_L1_transport= strrep(strrep(strrep(FABs_L1_reactants, ' ', ''),'[c]',''),'[m]','');
FABs_L1_transport=[FABs_L1_transport+"[m]"
    FABs_L1_transport+"[c]"];



GSH_L2_reactants=rmmissing(GSH_L2.Reactants, 'MinNumMissing', 1);
GSH_L2_reactants_parsed=[];
for i=1:length(GSH_L2_reactants)
   list=split(GSH_L2_reactants(i),';');
   for j=1:length(list)
    GSH_L2_reactants_parsed=[GSH_L2_reactants_parsed;strrep(list(j), ' ', '')];
   end
end
GSH_L2_reactants_parsed=[GSH_L2_reactants_parsed;GSH_L1_transport];
GSH_L2_reactants_parsed=unique(GSH_L2_reactants_parsed);


UREA_L2_reactants=rmmissing(UREA_L2.Reactants, 'MinNumMissing', 1);
UREA_L2_reactants_parsed=[];
for i=1:length(UREA_L2_reactants)
   list=split(UREA_L2_reactants(i),';');
   for j=1:length(list)
    UREA_L2_reactants_parsed=[UREA_L2_reactants_parsed;strrep(list(j), ' ', '')];
   end
end
UREA_L2_reactants_parsed=[UREA_L2_reactants_parsed;UREA_L1_transport];
UREA_L2_reactants_parsed=unique(UREA_L2_reactants_parsed);

MET_L2_reactants=rmmissing(MET_L2.Reactants, 'MinNumMissing', 1);
MET_L2_reactants_parsed=[];
for i=1:length(MET_L2_reactants)
   list=split(MET_L2_reactants(i),';');
   for j=1:length(list)
    MET_L2_reactants_parsed=[MET_L2_reactants_parsed;strrep(list(j), ' ', '')];
   end
end
MET_L2_reactants_parsed=[MET_L2_reactants_parsed;MET_L1_transport];
MET_L2_reactants_parsed=unique(MET_L2_reactants_parsed);


FABs_L2_reactants=rmmissing(FABs_L2.Reactants, 'MinNumMissing', 1);
FABs_L2_reactants_parsed=[];
for i=1:length(FABs_L2_reactants)
   list=split(FABs_L2_reactants(i),';');
   for j=1:length(list)
    FABs_L2_reactants_parsed=[FABs_L2_reactants_parsed;strrep(list(j), ' ', '')];
   end
end
FABs_L2_reactants_parsed=[FABs_L2_reactants_parsed;FABs_L1_transport];
FABs_L2_reactants_parsed=unique(FABs_L2_reactants_parsed);


%original definition

ADJ_mat=zeros(4,4); %GSH MET UREA FABs
ADJ_mat(1,2)=length(intersect(GSH_L1_reactants,MET_L2_reactants_parsed));
ADJ_mat(1,3)=length(intersect(GSH_L1_reactants,UREA_L2_reactants_parsed));
ADJ_mat(1,4)=length(intersect(GSH_L1_reactants,FABs_L2_reactants_parsed));


ADJ_mat(2,1)=length(intersect(MET_L1_reactants,GSH_L2_reactants_parsed));
ADJ_mat(2,3)=length(intersect(MET_L1_reactants,UREA_L2_reactants_parsed));
ADJ_mat(2,4)=length(intersect(MET_L1_reactants,FABs_L2_reactants_parsed));


ADJ_mat(3,1)=length(intersect(UREA_L1_reactants,GSH_L2_reactants_parsed));
ADJ_mat(3,2)=length(intersect(UREA_L1_reactants,MET_L2_reactants_parsed));
ADJ_mat(3,4)=length(intersect(UREA_L1_reactants,FABs_L2_reactants_parsed));


ADJ_mat(4,1)=length(intersect(FABs_L1_reactants,GSH_L2_reactants_parsed));
ADJ_mat(4,2)=length(intersect(FABs_L1_reactants,MET_L2_reactants_parsed));
ADJ_mat(4,3)=length(intersect(FABs_L1_reactants,UREA_L2_reactants_parsed));


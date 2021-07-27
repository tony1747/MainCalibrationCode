[num,txt,raw]=xlsread('fitted params (PBT138, slow, k_n1 included).xlsx');
choices=txt(:,1);
choices(1,:)=[];
m=size(choices);
m=m(1);
for i=1:m
    choice=choices(i);
    choice=char(choice);
    multiple_binding_solver(choice);
end

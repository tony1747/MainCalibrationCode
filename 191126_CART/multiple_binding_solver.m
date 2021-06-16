function multiple_binding_solver(choice)%for example, choice='H, 1:10, 1'
[num,txt,raw]=xlsread('fitted params(PBT138, slow).xlsx');

row=find(strcmp(txt,choice));
params=num(row,:);
receptor=choice(1);
ratio=extractBetween(choice,' ',',','Boundaries','exclusive');%extract 1:10 out of the choice
if string(ratio{1})=='1:5'
    CARTnum=5;    
elseif string(ratio{1})=='1:10'
    CARTnum=10;   
elseif string(ratio{1})=='1:20'
    CARTnum=20;
end
CARTratio=1/CARTnum;
read_CARTdata()
n=str2double(choice(end));  %choose the index of 'H, 1:10, 1', 
data=selected_cancerdata(n);
x(1)=data(1);
x(2)=data(1)*CARTratio;
x(3:4)=0;
[t,tumVol] = ode23(@(t,y)two_binding_slow(t,y,param), x, time, params );
function multiple_binding_solver(string)%for example, string='H, 1:10, 1'
[num,txt,raw]=xlsread('fitted params.xlsx');

row=find(strcmp(txt,string));
params=num(row,:);
receptor=string(1);
ratio=extractBetween(string,' ',',','Boundaries','exclusive');
if ratio{1}=='1:5'
    CARTnum=5;    
elseif ratio{1}=='1:10'
    CARTnum=10;   
elseif ratio{1}=='1:20'
    CARTnum=20;
end
CARTratio=1/CARTnum;
read_CARTdata()
n=str2double(string(end));
data=selected_cancerdata(n);
x(1)=data(1);
x(2)=data(1)/CARTratio;
x(3:4)=0;
[t,tumVol] = ode23(@(t,y)two_binding_slow(t,y,param), x, time, params );
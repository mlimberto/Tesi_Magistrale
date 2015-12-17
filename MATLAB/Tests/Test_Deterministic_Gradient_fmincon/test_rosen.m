fun = @rosenbrockwithgrad;
x0 = [-1,2];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [-2,-2];
ub = [2 , 2];
nonlcon = [];

options = optimoptions('fmincon','GradObj','on');

x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
4784project
===========
constants 

gk = 36;
gna = 120;
gl = 0.3;
Ek = -12;
Ena = 115;
El = 10.6;
Vrest = -70;
Cm = 1.0;



%initializing 
t =0; %time in milisecond seconds 
Vm = 0;  %membrane
VecX = [0:100];
VecY = [0];

%grating Variables
am = .1*((25-Vm)/(exp((25-Vm)/10)-1));
Bm = 4*exp(-Vm/18);
an = .01*((10-Vm)/(exp((10-Vm)/10)-1));
Bn = .125*exp(-Vm/80);
ah = .07*exp(-Vm/20);
Bh = 1/(exp((30-Vm)/10)+1);
    
m = am/(am+Bm);
n = an/(an+Bn);
h = ah/(ah+Bh);

I = 0;
Ina = (m^3)*gna*h*(Vm-Ena);
Ik = (n^4)*gk*(Vm-Ek);
Il = gl*(Vm-El);
Iion = I-Ina-Ik-Il;



ss = 1;
dVm = Iion/Cm;
Vm= Vm + ss*dVm;
dm = am*(1-m)-Bm*m;
dn = an*(1-n)-Bn*n;
dh = ah*(1-h)-Bh*h;
m = m + ss*dm;
n = n + ss*dn;
h = h + ss*dh;
       
t = t +1
%loop
while t <= 100 
    %grating Variables
    am = .1*((25-Vm)/(exp((25-Vm)/10)-1));
    Bm = 4*exp(-Vm/18);
    an = .01*((10-Vm)/(exp((10-Vm)/10)-1));
    Bn = .125*exp(-Vm/80);
    ah = .07*exp(-Vm/20);
    Bh = 1/(exp((30-Vm)/10)+1);
    
    


    %currents
    I = 0;
    Ina = (m^3)*gna*h*(Vm-Ena);
    Ik = (n^4)*gk*(Vm-Ek);
    Il = gl*(Vm-El);
    Iion = I-Ina-Ik-Il;
    
    %derivatives
    
    dVm = Iion/Cm;
    
    
    
    %updating stuff
    Vm= Vm + ss*dVm;
    
    dm = am*(1-m)-Bm*m;
    dn = an*(1-n)-Bn*n;
    dh = ah*(1-h)-Bh*h;
    m = m + ss*dm;
    n = n + ss*dn;
    h = h + ss*dh;
    
    VecY = [VecY Vm];
     t = t+1;
     
end
VecY = VecY-70;
VoltageGraph = plot(VecX,VecY,'m-')
title('Voltage during Action potential')
xlabel('Time (mili seconds)')
ylabel('Voltage (milivolts)')l')
xlabel('Time (mili seconds)')
ylabel('Voltage (milivolts)')

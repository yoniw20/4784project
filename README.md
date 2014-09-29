4784project
===========

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%question 1: representing a steady state neuron's voltage and chanel
%conductances over time

%constants 
gk = [36]; %potassium channel conductance
gna = [120]; %Sodium channel conductance
gl = 0.3; %leak channel conductance
Ek = -12; %Voltage across potassium channel
Ena = 115; %Voltage across soduim channel
El = 10.6; %Voltage across 
Vrest = -70; %Membrane resting voltage
Cm = 1.0; %membrane capacitance 



%initializing 
t =0; %time in milisecond seconds 
Vm = 0;  %membrane
VecX = [0:100];
VecY = [0];
gna_vec = (gna);
gk_vec = (gk);
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

%currents
I = 0;
Ina = (m^3)*gna*h*(Vm-Ena);
Ik = (n^4)*gk*(Vm-Ek);
Il = gl*(Vm-El);
Iion = I-Ina-Ik-Il;



ss = .01; %stepsize 

dVm = Iion/Cm; %change in voltage
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
    am = .1.*((25-Vm)/(exp((25-Vm)/10)-1));
    Bm = 4.*exp(-Vm/18);
    an = .01*((10-Vm)/(exp((10-Vm)/10)-1));
    Bn = .125.*exp(-Vm/80);
    ah = .07*exp(-Vm/20);
    Bh = 1./(exp((30-Vm)/10)+1);
    
    


    %currents
    I = 0;
    Ina = (m^3).*gna.*h.*(Vm-Ena);
    Ik = (n^4).*gk.*(Vm-Ek);
    Il = gl.*(Vm-El);
    Iion = I-Ina-Ik-Il;
    
    %derivatives
    
    dVm = Iion./Cm;
    
    
    
    %updating stuff
    Vm= Vm + ss.*dVm;
    
    dm = am*(1-m)-Bm*m;
    dn = an*(1-n)-Bn*n;
    dh = ah*(1-h)-Bh*h;
    m = m + ss*dm;
    n = n + ss*dn;
    h = h + ss*dh;
    
    VecY = [VecY Vm];
     t = t+1;
     
     %Conductance vectors
     gna_vec = [gna_vec gna];
     gk_vec = [gk_vec gk];
     
end

VecY = VecY-70;

figure
VoltageGraph = plot(VecX,VecY,'m-')
title('Voltage during Resting potential(question 1')
axis([0,100,-100,100])
xlabel('Time (mili seconds)')
ylabel('Voltage (milivolts)')

figure
conductance_na=plot(VecX,gna_vec,'b');
hold on
conductance_k= plot(VecX, gk_vec, 'r');
title('Channel Conductances for steady state neuron')
legend([conductance_na, conductance_k], 'conductance for Na+', 'conductance for K+')
axis([0,100,-100,150])
xlabel('Time (miliseconds)')
ylabel('Coductance (mS/cm^2)')

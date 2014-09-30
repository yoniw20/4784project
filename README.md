4784project
===========



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%question 1: representing a steady state neuron's voltage and chanel
%conductances over time

%constants 
gk = [36]; %potassium channel conductance (mS/cm^2)
gna = [120]; %Sodium channel conductance (mS/cm^2)
gl = 0.3; %leak channel conductance (mS/cm^2)
Ek = -12; %Voltage across potassium channel (mV)
Ena = 115; %Voltage across soduim channel(mV)
El = 10.6; %Voltage across leak channel (mV)
Vrest = -70; %Membrane resting voltage (mV)
Cm = 1.0; %membrane capacitance (microFarads/cm^2)



%initializing stuff
t =0; %time in milisecond seconds 
Vm = 0;  %membrane
VecX = [0:.01:100]; %time vector
VecY = [0]; %voltage vector. going to re-agust values at end by subtrating -70
gna_vec = (gna); %cunductance vectors
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
I = 0; %initial current
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
VecY = [VecY Vm];
       
t = t + .01
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
    I =0; %at steady state current = 0
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
     t = t+.01;
     
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%question 2: representing a neuron's voltage and chanel during action
%potential
clear

%constants 
gki = 36; %potassium channel conductance (mS/cm^2)
gnai = 120; %Sodium channel conductance (mS/cm^2)
gl = 0.3; %leak channel conductance (mS/cm^2)
Ek = -12; %Voltage across potassium channel (mV)
Ena = 115; %Voltage across soduim channel(mV)
El = 10.6; %Voltage across leak channel (mV)
Vrest = -70; %Membrane resting voltage (mV)
Cm = 1.0; %membrane capacitance (microFarads/cm^2)



%initializing stuff
t =0; %time in milisecond seconds 
Vm = 0;  %membrane
VecX = [0:.01:100]; %time vector
VecY = [0]; %voltage vector. going to re-adgust values at end by 
%subtrating -70

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
Ina = (m^3)*gnai*h*(Vm-Ena);
Ik = (n^4)*gki*(Vm-Ek);
Il = gl*(Vm-El);
Iion = I-Ina-Ik-Il;



ss = .01; %stepsize 

%deivatives
dVm = Iion/Cm; %change in voltage
Vm= Vm + ss*dVm;
dm = am*(1-m)-Bm*m;
dn = an*(1-n)-Bn*n;
dh = ah*(1-h)-Bh*h;

%updates
m = m + ss*dm;
n = n + ss*dn;
h = h + ss*dh;
gna= m.^3*h*gnai; %updating cunductances
gk= n.^4*gki;

gna_vec = (gna); %cunductance vectors
gk_vec = (gk);
       
t = t +.01
%loop
while t <= 100 
    %grating Variables
    am = .1.*((25-Vm)/(exp((25-Vm)/10)-1));
    Bm = 4.*exp(-Vm/18);
    an = .01*((10-Vm)/(exp((10-Vm)/10)-1));
    Bn = .125.*exp(-Vm/80);
    ah = .07*exp(-Vm/20);
    Bh = 1./(exp((30-Vm)/10)+1);
    
    if t <50 %this is the step pulse for .5 ms
        I = 5 % membrane current = 5 microA/cm^2 
    else
        I = 0 % goes back to steadystate
    end


    %currents
    
    Ina = (m^3).*gnai.*h.*(Vm-Ena);
    Ik = (n^4).*gki.*(Vm-Ek);
    Il = gl.*(Vm-El);
    Iion = I-Ina-Ik-Il;
    
    %derivatives
    
    dVm = Iion./Cm;
    
    
    
    %updating stuff
    Vm= Vm + ss.*dVm;
    
    %derivatives
    dm = am*(1-m)-Bm*m;
    dn = an*(1-n)-Bn*n;
    dh = ah*(1-h)-Bh*h;
    
    
    %updates
    m = m + ss*dm;
    n = n + ss*dn;
    h = h + ss*dh;
    gna= m.^3*h*gnai; %updating cunductances
    gk= n.^4*gki;
    VecY = [VecY Vm];
     t = t+.01;
     
     %Conductance vectors
     gna_vec = [gna_vec gna];
     gk_vec = [gk_vec gk];
     
end

VecY = VecY-70;

%%%%%making graphs
figure
VoltageGraph = plot(VecX,VecY,'m-')
title('Voltage during action potential(question 2)')
axis([0,100,-100,100])
xlabel('Time (mili seconds)')
ylabel('Voltage (milivolts)')

figure
conductance_na=plot(VecX,gna_vec,'b');
hold on
conductance_k= plot(VecX, gk_vec, 'r');
title('Channel Conductances for neuron with pulse current')
legend([conductance_na, conductance_k], 'conductance for Na+', 'conductance for K+')
axis([0,100,-100,150])
xlabel('Time (miliseconds)')
ylabel('Coductance (mS/cm^2)')


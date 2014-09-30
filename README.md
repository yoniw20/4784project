4784project

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%question 1: representing a steady state neuron's voltage and chanel
%conductances over time

%constants 
gk_max = 36; %max potassium channel conductance (mS/cm^2)
gna_max = 120; %max Sodium channel conductance (mS/cm^2)
gl = 0.3; %leak channel conductance (mS/cm^2)
Ek = -12; %Voltage across potassium channel (mV)
Ena = 115; %Voltage across soduim channel(mV)
El = 10.6; %Voltage across leak channel (mV)
Vrest = -70; %Membrane resting voltage (mV)
Cm = 1.0; %membrane capacitance (microFarads/cm^2)

%%%%%%%%%%%initializing stuff
t =0; %time in milisecond seconds 
Vm = 0;  %membrane

%%initializing graph stuff

gk_vec = [];
VecX = [0:.01:100]; %time vector
VecY = [Vm]
ss = .01 %step size = .01ms

%grating Variables
am = .1*((25-Vm)/(exp((25-Vm)/10)-1));
Bm = 4*exp(-Vm/18);
an = .01*((10-Vm)/(exp((10-Vm)/10)-1));
Bn = .125*exp(-Vm/80);
ah = .07*exp(-Vm/20);
Bh = 1/(exp((30-Vm)/10)+1);

%initial activation probabilities 
m = am/(am+Bm);
n = an/(an+Bn);
h = ah/(ah+Bh);


gna= m^3*h*gna_max;
gk= n^4*gk_max;
gna_vec = [gna];
gk_vec = [gk]





while t <=100
    %grating Variables update
    am = .1.*((25-Vm)/(exp((25-Vm)/10)-1));
    Bm = 4.*exp(-Vm/18);
    an = .01*((10-Vm)/(exp((10-Vm)/10)-1));
    Bn = .125.*exp(-Vm/80);
    ah = .07*exp(-Vm/20);
    Bh = 1./(exp((30-Vm)/10)+1);
    
    %current stuff
    I =0; %at steady state current = 0
    Ina = (m^3).*gna_max.*h.*(Vm-Ena);
    Ik = (n^4).*gk_max.*(Vm-Ek);
    Il = gl.*(Vm-El);
    Iion = I-Ina-Ik-Il;
    
    %%%%%%%%%%updates
    dVm = Iion./Cm; %change in voltage over change in time
    Vm= Vm + ss.*dVm; %updating voltage
    
    %change in activation pobababilties over time
    dm = am*(1-m)-Bm*m; 
    dn = an*(1-n)-Bn*n;
    dh = ah*(1-h)-Bh*h;
    %updating activation probabilities 
    m = m + ss*dm;
    n = n + ss*dn;
    h = h + ss*dh;
    
    %updating voltage vector
    VecY = [VecY Vm];
    
    gna= m^3*h*gna_max;
    gk= n^4*gk_max;
    gna_vec = [gna_vec gna];
	gk_vec = [gk_vec gk];
    
    
    t = t +.01 %updating time
end

VecY = VecY-70;

%%%%%%Plotting voltage
figure
VoltageGraph = plot(VecX,VecY,'m-')
title('Voltage during Resting potential(question 1')
axis([0,100,-100,100])
xlabel('Time (mili seconds)')
ylabel('Voltage (milivolts)')


%%%%%%plotting conductances
figure
conductance_na=plot(VecX,gna_vec,'b');
hold on
conductance_k= plot(VecX, gk_vec, 'r');
title('Channel Conductances for steady state neuron')
legend([conductance_na, conductance_k], 'conductance for Na+', 'conductance for K+')
axis([0,100,0,5])
xlabel('Time (miliseconds)')
ylabel('Coductance (mS/cm^2)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%question 2: representing neuron's voltage and chanel
%conductances over time when stimulated with a pulse 

clear
%constants 
gk_max = 36; %max potassium channel conductance (mS/cm^2)
gna_max = 120; %max Sodium channel conductance (mS/cm^2)
gl = 0.3; %leak channel conductance (mS/cm^2)
Ek = -12; %Voltage across potassium channel (mV)
Ena = 115; %Voltage across soduim channel(mV)
El = 10.6; %Voltage across leak channel (mV)
Vrest = -70; %Membrane resting voltage (mV)
Cm = 1.0; %membrane capacitance (microFarads/cm^2)

%%%%%%%%%%%initializing stuff
t =0; %time in milisecond seconds 
Vm = 0;  %membrane

%%initializing graph stuff

gk_vec = [];
VecX = [0:.01:100]; %time vector
VecY = [Vm]
ss = .01 %step size = .01ms

%grating Variables
am = .1*((25-Vm)/(exp((25-Vm)/10)-1));
Bm = 4*exp(-Vm/18);
an = .01*((10-Vm)/(exp((10-Vm)/10)-1));
Bn = .125*exp(-Vm/80);
ah = .07*exp(-Vm/20);
Bh = 1/(exp((30-Vm)/10)+1);

%initial activation probabilities 
m = am/(am+Bm);
n = an/(an+Bn);
h = ah/(ah+Bh);


gna= m^3*h*gna_max;
gk= n^4*gk_max;
gna_vec = [gna];
gk_vec = [gk]




while t<=100
    %grating Variables update
    am = .1.*((25-Vm)/(exp((25-Vm)/10)-1));
    Bm = 4.*exp(-Vm/18);
    an = .01*((10-Vm)/(exp((10-Vm)/10)-1));
    Bn = .125.*exp(-Vm/80);
    ah = .07*exp(-Vm/20);
    Bh = 1./(exp((30-Vm)/10)+1);
    
    %current stuff
    if t <= .5   %pulse for .5 seconds
        I = 5   
    else 
        I = 0 %at steady state current = 0
        
    end 
    Ina = (m^3).*gna_max.*h.*(Vm-Ena);
    Ik = (n^4).*gk_max.*(Vm-Ek);
    Il = gl.*(Vm-El);
    Iion = I-Ina-Ik-Il;
    
    %%%%%%%%%%updates
    dVm = Iion./Cm; %change in voltage over change in time
    Vm= Vm + ss.*dVm; %updating voltage
    
    %change in activation pobababilties over time
    dm = am*(1-m)-Bm*m; 
    dn = an*(1-n)-Bn*n;
    dh = ah*(1-h)-Bh*h;
    %updating activation probabilities 
    m = m + ss*dm;
    n = n + ss*dn;
    h = h + ss*dh;
    
    %updating voltage vector
    VecY = [VecY Vm];
    
    gna= m^3*h*gna_max;
    gk= n^4*gk_max;
    gna_vec = [gna_vec gna];
	gk_vec = [gk_vec gk];
    
    
    t = t +.01 %updating time
end

VecY = VecY-70;

%%%%%%Plotting voltage
figure
VoltageGraph = plot(VecX,VecY,'m-')
title('Voltage during Resting potential(question 1')
axis([0,100,-100,100])
xlabel('Time (mili seconds)')
ylabel('Voltage (milivolts)')


%%%%%%plotting conductances
figure
conductance_na=plot(VecX,gna_vec,'b');
hold on
conductance_k= plot(VecX, gk_vec, 'r');
title('Channel Conductances for steady state neuron')
legend([conductance_na, conductance_k], 'conductance for Na+', 'conductance for K+')
axis([0,100,0,40])
xlabel('Time (miliseconds)')
ylabel('Coductance (mS/cm^2)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%question 3: representing neuron's voltage and chanel
%conductances over time when stimulated by a constant current

clear
%constants 
gk_max = 36; %max potassium channel conductance (mS/cm^2)
gna_max = 120; %max Sodium channel conductance (mS/cm^2)
gl = 0.3; %leak channel conductance (mS/cm^2)
Ek = -12; %Voltage across potassium channel (mV)
Ena = 115; %Voltage across soduim channel(mV)
El = 10.6; %Voltage across leak channel (mV)
Vrest = -70; %Membrane resting voltage (mV)
Cm = 1.0; %membrane capacitance (microFarads/cm^2)

%%%%%%%%%%%initializing stuff
t =0; %time in milisecond seconds 
Vm = 0;  %membrane

%%initializing graph stuff

gk_vec = [];
VecX = [0:.01:100]; %time vector
VecY = [Vm]
ss = .01 %step size = .01ms

%grating Variables
am = .1*((25-Vm)/(exp((25-Vm)/10)-1));
Bm = 4*exp(-Vm/18);
an = .01*((10-Vm)/(exp((10-Vm)/10)-1));
Bn = .125*exp(-Vm/80);
ah = .07*exp(-Vm/20);
Bh = 1/(exp((30-Vm)/10)+1);

%initial activation probabilities 
m = am/(am+Bm);
n = an/(an+Bn);
h = ah/(ah+Bh);


gna= m^3*h*gna_max;
gk= n^4*gk_max;
gna_vec = [gna];
gk_vec = [gk]





while t <=100
    %grating Variables update
    am = .1.*((25-Vm)/(exp((25-Vm)/10)-1));
    Bm = 4.*exp(-Vm/18);
    an = .01*((10-Vm)/(exp((10-Vm)/10)-1));
    Bn = .125.*exp(-Vm/80);
    ah = .07*exp(-Vm/20);
    Bh = 1./(exp((30-Vm)/10)+1);
    
    %current stuff
    I =5; %at steady state current = 0 but now it's stimulated
    Ina = (m^3).*gna_max.*h.*(Vm-Ena);
    Ik = (n^4).*gk_max.*(Vm-Ek);
    Il = gl.*(Vm-El);
    Iion = I-Ina-Ik-Il;
    
    %%%%%%%%%%updates
    dVm = Iion./Cm; %change in voltage over change in time
    Vm= Vm + ss.*dVm; %updating voltage
    
    %change in activation pobababilties over time
    dm = am*(1-m)-Bm*m; 
    dn = an*(1-n)-Bn*n;
    dh = ah*(1-h)-Bh*h;
    %updating activation probabilities 
    m = m + ss*dm;
    n = n + ss*dn;
    h = h + ss*dh;
    
    %updating voltage vector
    VecY = [VecY Vm];
    
    gna= m^3*h*gna_max;
    gk= n^4*gk_max;
    gna_vec = [gna_vec gna];
	gk_vec = [gk_vec gk];
    
    
    t = t +.01 %updating time
end

VecY = VecY-70;

%%%%%%Plotting voltage
figure
VoltageGraph = plot(VecX,VecY,'m-')
title('Voltage during Resting potential(question 1')
axis([0,100,-100,100])
xlabel('Time (mili seconds)')
ylabel('Voltage (milivolts)')


%%%%%%plotting conductances
figure
conductance_na=plot(VecX,gna_vec,'b');
hold on
conductance_k= plot(VecX, gk_vec, 'r');
title('Channel Conductances for steady state neuron')
legend([conductance_na, conductance_k], 'conductance for Na+', 'conductance for K+')
axis([0,100,-100,100])
xlabel('Time (miliseconds)')
ylabel('Coductance (mS/cm^2)')

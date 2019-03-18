%effective lattice constant. 1:6 is discharge, 6:12 is charge. 
aeff4=[
    8.098841292
8.082505276
8.091335852
8.144877527
8.137046808
8.183136774
8.190092491
8.177716293
8.177138676
8.142470778
8.155375217
8.142351682
8.140935712
8.141033592
8.187651256
8.172903312
8.132190681
8.173910569
8.180330397
8.175844893
8.169025791
8.114736231
8.098563875
8.143313028
8.174397507
8.193170026
8.200700122];
%delta corresponding to lattice constant. 1:6 is discharge, 6:12 is charge
exp4ddat=[
0.53844
0.21515
0.1786
0.14955
0.11008
1.39e-04
0.05268
0.11342
0.37208
0.72746
0.86812
0.89536
0.06519
0.55164
0.03544
0.50202
0.07728
1.39e-04
0.05893
0.10702
0.37208
0.72746
0.49925
0.21515
0.16484
0.06519
0.03878];

%josh PXRD data for charge
josh_delta=[
0
0.21307
0.41
0.49
0.64254
0.7384
0.841
0.888];

josh_lat1=[
	8.1782
8.177263
8.16835
8.1405
8.139617
8.146
8.132
8.138];

josh_lat2=[
	8.0904
8.089
8.0747
8.0783];

josh_lat3=[
	8.00195
8.000195
8.00538];

%josh PXRD data for discharge
josh_delta_d=[0.809, 0.548, 0.486, 0.361, 0.295, 0.163, .052];

josh_lat1_d=[8.120433, 8.1286, 8.1283, 8.1344, 8.13,8.156, 8.166];
    
josh_lat2_d=[8.0904,8.0939]; 

josh_lat3_d=[8.0127];


%% Make my own PXRD graph

hold on
%charge data
plot(exp4ddat(7:12),aeff4(7:12),'bs') %my individual particle charge data
plot(exp4ddat(1:6),aeff4(1:6),'kd') %my individual particle discharge
plot(josh_delta,josh_lat1,'rh') %josh phase 1
plot(josh_delta(end-3:end),josh_lat2,'rh')  %josh phase 2
plot(josh_delta(end-2:end),josh_lat3,'rh')  %josh phaes 3
box on
saveas(1,'test.eps','psc2')
%discharge data
hold on
plot(josh_delta_d,josh_lat1_d,'kx') %josh phase 1
plot(josh_delta_d(1:2),josh_lat2_d,'kx') %josh phase 2
plot(josh_delta_d(1),josh_lat3_d,'kx') %%josh phaes 3
plot(exp4ddat(1:6),aeff4(1:6),'r+') %my individual discharge data
xlabel('\delta in Li_{1-\delta}')
ylabel('Lattice Constant (Angstrom)')


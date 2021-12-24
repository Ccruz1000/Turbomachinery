# Function calculates both specific thrust and fuel consumption without having to loop through twice
function [specific_thrust, SFC] = SFC_Calc (Ca, ta, pa, Ai, B, FPR, OPR, TIT, eta_m_high_pressure, eta_m_low_pressure, eta_i, eta_inf_c, eta_inf_f, eta_inf_lpt, eta_inf_hpt, eta_comb, deltap_comb, eta_j_hot, eta_j_cold, LHV)

r=OPR/FPR;
mdot=(pa/(287*ta))*Ai*Ca;
mdot_hot=mdot/(B+1);
mdot_cold=B*mdot_hot;
p01=(pa)*(1+(eta_i*Ca^2)/(2*1005*ta))^(1.4/(1.4-1));
t01=ta+Ca^2/(2*1005);
p02=p01*FPR;
t02=t01*FPR^(0.4/(1.4*eta_inf_f));
p03=p02*r;
t03=t02*r^(0.4/(1.4*eta_inf_c));
p04=p03-deltap_comb;
t04=TIT;
wdot_comp=mdot_hot*1005*(t03-t02);
wdot_hpt=wdot_comp/eta_m_high_pressure;
temp_drop_hpt=wdot_hpt/(mdot_hot*1148);
t05=t04-temp_drop_hpt;
high_pressure_tpr=(t05/t04)^(1.333/(0.333*eta_inf_hpt));
p05=p04*high_pressure_tpr;
wdot_fan=mdot*1005*(t02-t01);
wdot_lpt=wdot_fan/eta_m_low_pressure;
temp_drop_lpt=wdot_lpt/(mdot_hot*1148);
t06=t05-temp_drop_lpt;
low_pressure_tpr=(t06/t05)^(1.333/(0.333*eta_inf_lpt));
p06=p05*low_pressure_tpr;
ambient_to_p06=pa/p06;
critical_to_p06=(1-(1/eta_j_hot)*(0.333/2.333))^(1.333/0.333);
if(ambient_to_p06<critical_to_p06)
p7=p06*critical_to_p06;
t7=t06*(2/(2.333));
Cj_hot=sqrt(1.333*287*t7);
else
p7=pa;
t7=((((p7/p06)^(0.333/1.333)-1)*-eta_j_hot)-1)*(-t06);
Cj_hot= sqrt(2*1148*eta_j_hot*t06(1-(p7/p06)^(0.333/1.333)));
endif
ambient_to_p02=pa/p02;
critical_to_p02=(1-(1/eta_j_cold)*(0.4/2.4))^(1.4/0.4);
if(ambient_to_p02<critical_to_p02)
p8=p02*critical_to_p02;
t8=t02*(2/2.4);
Cj_cold=sqrt(1.4*287*t8);
else
p8=pa;
t8=((((p8/p02)^(0.4/1.4)-1)*-eta_j_cold)-1)*(-t02);
Cj_cold=sqrt(2*1005*eta_j_cold*t02(1-(p8/p02)^(0.4/1.4)));
endif
Ah=mdot_hot/((p7/(287*t7))*Cj_hot);
Ac=mdot_cold/((p8/(287*t8))*Cj_cold);
hot_thrust=mdot_hot*(Cj_hot-Ca)+(p7-pa)*Ah;
cold_thrust=mdot_cold*(Cj_cold-Ca)+(p8-pa)*Ac;
intake_momentum_drag=mdot*Ca;
net_thrust=hot_thrust+cold_thrust
specific_thrust=net_thrust/mdot
f_theoretical=((1148*t04-1005*t03)/(LHV-1148*t04));
f_actual=f_theoretical/eta_comb
fuel_flow=f_actual*mdot_hot
SFC=(fuel_flow/net_thrust)*3600
eta_overall=(1000*Ca)/(SFC*LHV)
 
Endfunction

# Initialize Test Case
LHV=43100000;
eta_m_low_pressure=0.99;
eta_m_high_pressure=0.99;
Ca=221.325;
ta=216.7;
pa=17930;
Ai=0.150453;
eta_i=0.95;
r=7.9;
OPR = 1.7*7.9
eta_inf_c=0.9;
TIT=1220;
eta_inf_lpt=0.9;
eta_inf_hpt=0.9;
B=3.8;
FPR=1.7;
eta_inf_f=0.9;
eta_j_hot=1;
eta_j_cold=1;
eta_comb=0.98;
deltap_comb=12040;

# Find optimal
function [optimal_FPR] = optimal_FPR(Ca, ta, pa, Ai, B, OPR, TIT, eta_m_high_pressure, eta_m_low_pressure, eta_i, eta_inf_c, eta_inf_f, eta_inf_lpt, eta_inf_hpt, eta_comb, deltap_comb, eta_j_hot, eta_j_cold, LHV)
# No need to make for loop to calculate fpr, can just declare fpr like this
fpr = 1.4:0.04:1.8;
# Allows changing the step for fpr from 0.04 without having to change other code if you want a different range   
len = length(fpr)  
# Initialize FPR_Curve so that it has the desired shape. Also saves time if this function needs to be called a bunch 
FPR_Curve = zeros(len, 3);
# For loop to make calculations. Might have to change the SFC_Calc input parameters because I’m not sure what you’re doing here
for(i=1:len)
FPR_Curve(i, 1) = fpr(i);
[x(i,2), x(i,3)] = SFC_Calc (Ca, ta, pa, Ai, B, FPR, OPR, TIT, eta_m_high_pressure, eta_m_low_pressure, eta_i, eta_inf_c, eta_inf_f, eta_inf_lpt, eta_inf_hpt, eta_comb, deltap_comb, eta_j_hot, eta_j_cold, LHV);
endfor
 # Determines maximum value, and row that max value occurs in
[maxval, row] = max(FPR_curve(:,2))
 
optimal_FPR=FPR_curve(row,1);
optimal_specific_thrust = FPR_curve(row,2);
optimal_SFC=FPR_curve(row,3);
Endfunction


from uncertainties import *

lumi = 3.189 * 1e15
f_d = ufloat(0.402, 0.007)
σ_bb = ufloat(295, 29) * 1e-6

BR_B2Dmumu = ufloat(2.6, 1.3) * 1e-9
BR_B2Dstmumu = ufloat(1.4, 0.5) * 1e-8
BR_B2Dpmumu = ufloat(2.5, 0.5) * 1e-12
BR_B2Dpstmumu = ufloat(1.0, 0.5) * 1e-11
BR_B2Dsmumu = ufloat(4.3, 0.5) * 1e-11
BR_B2Dsstmumu = ufloat(2, 0.5) * 1e-10
BR_Bc2Dsmumu = ufloat(1.0, 1.0) * 1e-6
BR_Bc2Dsstmumu = ufloat(5.0, 1.0) * 1e-6

eff_strip = 0.20
eff_geom = 0.15
eff_trig = 0.8

BR_Dst2Dpi = 0.619
BR_D2Kpi = 0.0387
BR_Dp2Kpipi = 0.0913
BR_Dp2Kpipi = 0.0913
BR_Ds2KKpi = 0.0549
BR_Dst2Dpi = 0.619
BR_Dstp2Dpi = 0.677
BR_Dsst2Dspi = 0.058

eff = eff_strip * eff_geom * eff_trig

print('B2Dmumu',    lumi * σ_bb * 2 * f_d * eff * BR_B2Dmumu * BR_D2Kpi)
print('B2Dstmumu',  lumi * σ_bb * 2 * f_d * eff * BR_B2Dstmumu * BR_Dst2Dpi * BR_D2Kpi)
print('B2Dpmumu',   lumi * σ_bb * 2 * f_d * eff * BR_B2Dpmumu * BR_Dp2Kpipi)
print('B2Dpstmumu', lumi * σ_bb * 2 * f_d * eff * BR_B2Dpstmumu * BR_Dstp2Dpi * BR_D2Kpi)
print('B2Dsmumu',   lumi * σ_bb * 2 * f_d * eff * BR_B2Dsmumu * BR_Ds2KKpi * BR_D2Kpi)
print('B2Dsstmumu', lumi * σ_bb * 2 * f_d * eff * BR_B2Dsstmumu * BR_Dsst2Dspi * BR_Ds2KKpi)
 

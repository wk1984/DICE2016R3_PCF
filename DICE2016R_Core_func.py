#!/usr/bin/env python
# coding: utf-8

# In[1]:



import pyomo.environ as pyo
import numpy as np

class DICE2016R_Core:
    
    def eqDDTpfEQ(m, t):
        return m.DDTpf[t] == 200.46*m.TATM[t] +1038.89
    
    def eqALTpfEQ(m, t):
        return m.ALTpf[t] == 0.045*m.DDTpf[t]**0.5
    
    def eqC_PFEQ(m,t):
        return m.C_PF[t] == -68.41*(m.ALTpf[t])**2+ 563*(m.ALTpf[t]) - 13.41
        
    def eqCthawedPFEQ(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.CthawedPF[t] == m.C_PF[t]-m.C_PF[t-1]

    def eqCCumPFEQ(m,t):
        if t == 1:
            return pyo.Constraint.Skip
        else:        
            list6=np.arange(1,t)
            list7=[0]
            for s in list6:
                z=m.CthawedPF[s]
                x=(t-s)*5/70
                e=pyo.exp(x)
                y=0.6*(z)*(1-1/e)
                list7.append(y)
            co2em=sum(list7)
            return m.CCumPF[t] == co2em

    def eqCO2_PFEQ(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.CO2_PF[t] == 3.666*(1-0.023)*(m.CCumPF[t]- m.CCumPF[t-1])

    def eqCH4_PFEQ(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.CH4_PF[t] == 1.333*0.023*(m.CCumPF[t]- m.CCumPF[t-1])

    def eqE_PFEQ(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.E_PF[t] == m.CO2_PF[t] + m.RE * m.CH4_PF[t]
        
    def eqEEQ(m, t):
        return m.E[t] == m.EIND[t] + m.etree[t] + m.E_PF[t]/5
    
    def eqEINDEQ(m, t):
        return m.EIND[t] == m.sigma[t] * m.YGROSS[t] * (1 - m.MIU[t])

    def eqMIUEQ(m, t):
        if t==1:
            return pyo.Constraint.Skip
        else:
            return m.MIU[t] <= m.MIU[t-1]+0.2
        
    def eqCCACCA(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.CCA[t] == m.CCA[t-1] + m.EIND[t-1] * m.tstep / 3.666

    def eqCCATOTEQ(m, t):
        return m.CCATOT[t] == m.CCA[t] + m.cumetree[t]

    def eqFORCEQ(m, t):       
        return m.FORC[t] == m.fco22x/pyo.log(2) * (pyo.log(m.MAT[t]/m.mateq)) + m.forcoth[t]

    def eqDAMFRACEQ(m, t):
        return m.DAMFRAC[t] == m.a1 * m.TATM[t] + m.a2 * m.TATM[t] ** m.a3

    def eqDAMEQ(m, t):
        return m.DAMAGES[t] == m.YGROSS[t] * m.DAMFRAC[t]

    def eqABATEEQ(m, t):
        return m.ABATECOST[t] == m.YGROSS[t] * m.cost1[t] * (m.MIU[t] ** m.expcost2)

    def eqMCABATEEQ(m, t):
        return m.MCABATE[t] == m.pbacktime[t] * m.MIU[t] ** (m.expcost2 - 1)

    def eqCARBPRICEEQ(m, t):
        return m.CPRICE[t] == m.pbacktime[t] * m.MIU[t] ** (m.expcost2 - 1)

    # Climate and carbon cycle =============

    def eqMMAT(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.MAT[t] == m.MAT[t-1] * m.b11 + m.MU[t-1]*m.b21 + (m.E[t-1]*m.tstep/3.666)

    def eqMML(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.ML[t] == m.ML[t-1] * m.b33 + m.MU[t-1] * m.b23

    def eqMMU(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.MU[t] == m.MAT[t-1] * m.b12 + m.MU[t-1] * m.b22 + m.ML[t-1] * m.b32

    def eqTATMEQ(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.TATM[t] == m.TATM[t-1] + m.c1 * (m.FORC[t] - m.fco22x/m.t2xco2*m.TATM[t-1] - m.c3 * (m.TATM[t-1]-m.TOCEAN[t-1]))

    def eqTOCEANEQ(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.TOCEAN[t] == m.TOCEAN[t-1] + m.c4 * (m.TATM[t-1] - m.TOCEAN[t-1])

    # Economic variables ===================

    def eqYGROSSEQ(m, t):
        return m.YGROSS[t] == (m.al[t]*(m.l[t]/1000)**(1-m.gama))*(m.K[t]**m.gama)

    def eqYNETEQ(m, t):
        return m.YNET[t] == m.YGROSS[t] * (1 - m.DAMFRAC[t])

    def eqYY(m, t):
        return m.Y[t] == m.YNET[t] - m.ABATECOST[t]

    def eqCC(m, t):
        return m.C[t] == m.Y[t] - m.I[t]

    def eqCPCE(m, t):
        return m.CPC[t] == 1000*m.C[t]/m.l[t]

    def eqSEQ(m, t):
        return m.I[t] == m.S[t]*m.Y[t]

    def eqKK(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.K[t] <= (1-m.dk)**m.tstep*m.K[t-1]+m.tstep*m.I[t-1]

    def eqRIEQ(m, t):
        if t == 1:
            return pyo.Constraint.Skip
        else:
            return m.RI[t-1] == (1+m.prstp) * (m.CPC[t]/m.CPC[t-1])**(m.elasmu/m.tstep) - 1

    # Utility ================

    def eqCEMUTOTPEREQ(m, t):
        return m.CEMUTOTPER[t] == m.PERIODU[t] * m.l[t] * m.rr[t]

    def eqPERIODUEQ(m, t):
        return m.PERIODU[t] == (m.CPC[t]**(1-m.elasmu)-1)/(1-m.elasmu) - 1

    def eqUTIL(m):
        return m.UTILITY == m.tstep * m.scale1 * pyo.summation(m.CEMUTOTPER) + m.scale2

    def eqOBJ(m):
        return m.UTILITY

    # ==========================================================================================================

    def initialize_parameters(dice):

        dice.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)

        dice.ifopt = pyo.Param(initialize=0, mutable=True,
                               doc='Indicator where optimized is 1 and base is 0')

        dice.tstep = pyo.Param(
            initialize=5, mutable=True, doc='Years per Period')
        dice.t_yr = pyo.Param(dice.t, initialize=2010,
                              mutable=True, doc='Years AD')
        dice.t0 = pyo.Param(initialize=2010, mutable=True, doc='Years onset')

        dice.fosslim = pyo.Param(initialize=6000, mutable=True,
                                 doc='Maximum cumulative extraction fossil fuels (GtC)')
        dice.elasmu = pyo.Param(initialize=1.45, mutable=True,
                                doc='Elasticity of marginal utility of consumption')
        dice.prstp = pyo.Param(initialize=0.015, mutable=True,
                               doc='Initial rate of social time preference per year')
        dice.gama = pyo.Param(initialize=0.300, mutable=True,
                              doc='Capital elasticity in production function')
        dice.pop0 = pyo.Param(initialize=7403, mutable=True,
                              doc='Initial world population 2015 (millions)')
        dice.popadj = pyo.Param(initialize=0.134, mutable=True,
                                doc='Growth rate to calibrate to 2050 pop projection')
        dice.popasym = pyo.Param(
            initialize=11500, mutable=True, doc='Asymptotic population (millions)')
        dice.dk = pyo.Param(initialize=0.100, mutable=True,
                            doc='Depreciation rate on capital (per year)')
        dice.q00 = pyo.Param(initialize=105.5, mutable=True,
                            doc='Initial world gross output 2015 (trill 2010 USD)')
        dice.k00 = pyo.Param(initialize=223, mutable=True,
                            doc='Initial capital value 2015 (trill 2010 USD)')
        dice.a00= pyo.Param(initialize=5.115, mutable=True,
                            doc='Initial level of total factor productivity')
        dice.ga0 = pyo.Param(initialize=0.076, mutable=True,
                             doc='Initial growth rate for TFP per 5 years')
        dice.dela = pyo.Param(initialize=0.005, mutable=True,
                              doc='Decline rate of TFP per 5 years')
        dice.gsigma1 = pyo.Param(
            initialize=-0.0152, mutable=True, doc='Initial growth of sigma (per year)')
        dice.dsig = pyo.Param(initialize=-0.001, mutable=True,
                              doc='Decline rate of decarbonization (per period)')
        dice.eland0 = pyo.Param(initialize=2.6, mutable=True,
                                doc='Carbon emissions from land 2015 (GtCO2 per year)')
        dice.deland = pyo.Param(initialize=0.115, mutable=True,
                                doc='Decline rate of land emissions (per period)')
        dice.e0 = pyo.Param(initialize=35.85, mutable=True,
                            doc='Industrial emissions 2015 (GtCO2 per year)')
        dice.miu0 = pyo.Param(initialize=0.03, mutable=True,
                              doc='Initial emissions control rate for base case 2015')
        dice.mat0 = pyo.Param(initialize=851, mutable=True,
                              doc='Initial Carbon in atmosphere 2015 (GtC)')
        dice.mu0 = pyo.Param(initialize=460, mutable=True,
                             doc='Initial Concentration in upper strata 2015 (GtC)')
        dice.ml0 = pyo.Param(initialize=1740, mutable=True,
                             doc='Initial Concentration in lower strata 2015 (GtC)')
        dice.mateq = pyo.Param(initialize=588, mutable=True,
                               doc='Equilibrium concentration atmosphere  (GtC)')
        dice.mueq = pyo.Param(initialize=360, mutable=True,
                              doc='Equilibrium concentration in upper strata (GtC)')
        dice.mleq = pyo.Param(initialize=1720, mutable=True,
                              doc='Equilibrium concentration in lower strata (GtC)')

        dice.b12 = pyo.Param(initialize=0.12, mutable=True,
                             doc='Carbon cycle transition matrix')
        dice.b23 = pyo.Param(initialize=0.007, mutable=True,
                             doc='Carbon cycle transition matrix')
        dice.b11 = pyo.Param(initialize=0, mutable=True,
                             doc='Carbon cycle transition matrix')
        dice.b21 = pyo.Param(initialize=0, mutable=True,
                             doc='Carbon cycle transition matrix')
        dice.b22 = pyo.Param(initialize=0, mutable=True,
                             doc='Carbon cycle transition matrix')
        dice.b32 = pyo.Param(initialize=0, mutable=True,
                             doc='Carbon cycle transition matrix')
        dice.b33 = pyo.Param(initialize=0, mutable=True,
                             doc='Carbon cycle transition matrix')
        dice.sig0 = pyo.Param(initialize=0, mutable=True,
                              doc='Carbon intensity 2010 (kgCO2 per output 2005 USD 2010)')

        dice.t2xco2 = pyo.Param(initialize=3.1, mutable=True,
                                doc='Equilibrium temp impact (oC per doubling CO2)')
        dice.fex0 = pyo.Param(initialize=0.5, mutable=True,
                              doc='2015 forcings of non-CO2 GHG (Wm-2)')
        dice.fex1 = pyo.Param(initialize=1.0, mutable=True,
                              doc='2100 forcings of non-CO2 GHG (Wm-2)')
        dice.fex_idx = pyo.Param(initialize=17, mutable=True)
        dice.tocean0 = pyo.Param(initialize=0.0068, mutable=True,
                                 doc='Initial lower stratum temp change (C from 1900)')
        dice.tatm0 = pyo.Param(initialize=0.85, mutable=True,
                               doc='Initial atmospheric temp change (C from 1900)')
        dice.c1 = pyo.Param(initialize=0.1005, mutable=True,
                            doc='Climate equation coefficient for upper level')
        dice.c3 = pyo.Param(initialize=0.088, mutable=True,
                            doc='Transfer coefficient upper to lower stratum')
        dice.c4 = pyo.Param(initialize=0.025, mutable=True,
                            doc='Transfer coefficient for lower level')
        dice.fco22x = pyo.Param(initialize=3.6813, mutable=True,
                                doc='Forcings of equilibrium CO2 doubling (Wm-2)')

        dice.a10 = pyo.Param(initialize=0, mutable=True,
                             doc='Initial damage intercept')
        dice.a20 = pyo.Param(initialize=0, mutable=True,
                             doc='Initial damage quadratic term')
        dice.a1 = pyo.Param(initialize=0, mutable=True, doc='Damage intercept')
        dice.a2 = pyo.Param(initialize=0.00236, mutable=True,
                            doc='Damage quadratic term')
        dice.a3 = pyo.Param(initialize=2.00, mutable=True,
                            doc='Damage exponent')

        dice.expcost2 = pyo.Param(
            initialize=2.6, mutable=True, doc='Exponent of control cost function')
        dice.pback0 = pyo.Param(initialize=550, mutable=True,
                               doc='Cost of backstop 2010$ per tCO2 2015')
        dice.gback = pyo.Param(initialize=0.025, mutable=True,
                               doc='Initial cost decline backstop cost per period')
        dice.limmiu = pyo.Param(
            initialize=1.2, mutable=True, doc='Upper limit on control rate after 2150')
        dice.tnopol = pyo.Param(initialize=45, mutable=True,
                                doc='Period before which no emissions controls base')
        dice.cprice00 = pyo.Param(
            initialize=2, mutable=True, doc='Initial base carbon price (2010$ per tCO2)')
        dice.gcprice = pyo.Param(
            initialize=0.02, mutable=True, doc='Growth rate of base carbon price per year')

        dice.scale1 = pyo.Param(initialize=0.0302455265681763,
                                mutable=True, doc='Multiplicative scaling coefficient')
        dice.scale2 = pyo.Param(
            initialize=-10993.704, mutable=True, doc='Additive scaling coefficient')

        dice.l = pyo.Param(dice.t, initialize=dice.pop0,
                           mutable=True, doc='Level of population and labor')
        dice.al = pyo.Param(dice.t, initialize=0, mutable=True,
                            doc='Level of total factor productivity')
        dice.sigma = pyo.Param(
            dice.t, initialize=0, mutable=True, doc='CO2-equivalent-emissions output ratio')
        dice.rr = pyo.Param(dice.t, initialize=0, mutable=True,
                            doc='Average utility social discount rate')
        dice.ga = pyo.Param(dice.t, initialize=0, mutable=True,
                            doc='Growth rate of productivity from')
        dice.forcoth = pyo.Param(dice.t, initialize=0, mutable=True,
                                 doc='Exogenous forcing for other greenhouse gases')
        # dice.gl       = pyo.Param(dice.t, initialize=0, mutable=True, doc = 'Growth rate of labor')
        dice.gcost1 = pyo.Param(
            initialize=0, mutable=True, doc='Growth of cost factor')
        dice.gsig = pyo.Param(dice.t, initialize=0, mutable=True,
                              doc='Change in sigma (cumulative improvement of energy efficiency)')
        dice.etree = pyo.Param(
            dice.t, initialize=0, mutable=True, doc='Emissions from deforestation')
        dice.cumetree = pyo.Param(
            dice.t, initialize=0, mutable=True, doc='Cumulative from land')
        dice.cost1 = pyo.Param(dice.t, initialize=0,
                               mutable=True, doc='Adjusted cost for backstop')
        dice.gfacpop = pyo.Param(dice.t, initialize=0,
                                 mutable=True, doc='Growth factor population')
        dice.pbacktime = pyo.Param(
            dice.t, initialize=0, mutable=True, doc='Backstop price')
        dice.optlrsav = pyo.Param(
            initialize=0, mutable=True, doc='Optimal long-run savings rate used for transversality')
        dice.scc = pyo.Param(dice.t, initialize=0,
                             mutable=True, doc='Social cost of carbon')
        dice.cpricebase = pyo.Param(
            dice.t, initialize=0, mutable=True, doc='Carbon price in base case')
        dice.photel = pyo.Param(dice.t, initialize=0, mutable=True,
                                doc='Carbon Price under no damages (Hotelling rent condition)')
        dice.ppm = pyo.Param(dice.t, initialize=0, mutable=True,
                             doc='Atmospheric concentrations parts per million')
        dice.atfrac = pyo.Param(
            dice.t, initialize=0, mutable=True, doc='Atmospheric share since 1850')
        dice.atfrac2010 = pyo.Param(
            dice.t, initialize=0, mutable=True, doc='Atmospheric share since 2010')
        dice.cumetree[1] = 100
        dice.RE = pyo.Param(initialize=29, mutable=True, doc='Radiative Efficiency')
       
        dice.p2018 = pyo.Param(initialize=1.2, mutable=True, doc='Price level 2018 relative to 2010') 
        return

    def precalculate_parameters(dice):

        dice.q0 = pyo.value(dice.q00*dice.p2018)
        dice.k0 = pyo.value(dice.k00*dice.p2018)
        dice.cprice0 = pyo.value(dice.cprice00*dice.p2018)
        dice.pback = pyo.value(dice.pback0*dice.p2018)
        dice.a0 = pyo.value(dice.a00*dice.p2018**((1-dice.gama)))

        dice.b11 = 1 - pyo.value(dice.b12)
        dice.b21 = pyo.value(dice.b12*dice.mateq/dice.mueq)
        dice.b22 = pyo.value(1 - dice.b21 - dice.b23)
        dice.b32 = pyo.value(dice.b23*dice.mueq/dice.mleq)
        dice.b33 = pyo.value(1 - dice.b32)

        dice.a20 = pyo.value(dice.a2)
        dice.sig0 = pyo.value(dice.e0/(dice.q0 * (1-dice.miu0)))

        dice.al[1] = pyo.value(dice.a0)
        dice.gsig[1] = pyo.value(dice.gsigma1)
        dice.sigma[1] = pyo.value(dice.sig0)
     

        for ttt in dice.t:
            dice.t_yr[ttt] = dice.t0 + dice.t.at(ttt) * dice.tstep
            dice.ga[ttt] = dice.ga0 * pyo.exp(-dice.dela * dice.tstep * (ttt-1))
            dice.etree[ttt] = dice.eland0 * (1-dice.deland) ** (ttt - 1)
            
        for ttt in np.arange(dice.t.first(), dice.t.last()):
            dice.l[ttt+1] = dice.l[ttt] * (dice.popasym/dice.l[ttt])**dice.popadj
            dice.al[ttt+1] = dice.al[ttt] / (1-dice.ga[ttt])
            dice.sigma[ttt+1] = dice.sigma[ttt] *np.exp(pyo.value(dice.gsig[ttt]*dice.tstep))
            dice.gsig[ttt+1] = dice.gsig[ttt] * ((1+dice.dsig)**dice.tstep)
            dice.cumetree[ttt+1] = dice.cumetree[ttt] + dice.etree[ttt] * dice.tstep / 3.666

        for ttt in np.arange(dice.t.first(), dice.t.last()+1):
            dice.pbacktime[ttt] = dice.pback * (1 - dice.gback) ** (ttt-1)
            dice.cost1[ttt] = dice.pbacktime[ttt] * dice.sigma[ttt]/dice.expcost2 / 1000
            dice.rr[ttt] = 1/((1+dice.prstp) ** (dice.tstep*(ttt-1)))
            dice.cpricebase[ttt] = dice.cprice0 * (1 + dice.gcprice) ** (dice.tstep*(ttt-1))
            if ttt <= pyo.value(dice.fex_idx):
                item2 = 1/dice.fex_idx * (dice.fex1-dice.fex0) * (ttt - 1)
            else:
                item2 = dice.fex1 - dice.fex0
            dice.forcoth[ttt] = dice.fex0 + item2

        dice.optlrsav = pyo.value(
            (dice.dk+0.004)/(dice.dk + 0.004 * dice.elasmu + dice.prstp) * dice.gama)

        return

    # ==========================================================================================================

    def initialize_variables(dice):

        dice.MIU = pyo.Var(dice.t, bounds=(0, +np.inf),
                           doc='Emission control rate GHGs')
        dice.FORC = pyo.Var(
            dice.t, doc='Increase in radiative forcing (watts per m2 from 1900)')
        dice.TATM = pyo.Var(dice.t, bounds=(
            0, +np.inf), doc='Increase temperature of atmosphere (degrees C from 1900)')
        dice.TOCEAN = pyo.Var(
            dice.t, doc='Increase temperature of lower oceans (degrees C from 1900)')
        dice.MAT = pyo.Var(dice.t, initialize=dice.mat0, bounds=(
            0, np.inf), doc='Carbon concentration increase in atmosphere (GtC from 1750)')
        dice.MU = pyo.Var(dice.t, bounds=(
            0, np.inf), doc='Carbon concentration increase in shallow oceans (GtC from 1750)')
        dice.ML = pyo.Var(dice.t, bounds=(
            0, np.inf), doc='Carbon concentration increase in lower oceans (GtC from 1750)')
        dice.E = pyo.Var(dice.t, doc='Total CO2 emissions (GtCO2 per year)')
        dice.EIND = pyo.Var(
            dice.t, doc='Industrial emissions (GtCO2 per year)')
        dice.C = pyo.Var(dice.t, bounds=(0, np.inf),
                         doc='Consumption (trillions 2005 US dollars per year)')
        dice.K = pyo.Var(dice.t, bounds=(0, np.inf),
                         doc='Capital stock (trillions 2005 US dollars)')
        dice.CPC = pyo.Var(
            dice.t, doc='Per capita consumption (thousands 2005 USD per year)')
        dice.I = pyo.Var(dice.t, bounds=(0, np.inf),
                         doc='Investment (trillions 2005 USD per year)')
        dice.S = pyo.Var(
            dice.t, doc='Gross savings rate as fraction of gross world product')
        dice.RI = pyo.Var(dice.t, doc='Real interest rate (per annum)')
        dice.Y = pyo.Var(dice.t, bounds=(
            0, np.inf), doc='Gross world product net of abatement and damages (trillions 2005 USD per year)')
        dice.YGROSS = pyo.Var(dice.t, bounds=(
            0, np.inf), doc='Gross world product (trillions 2005 USD per year)')
        dice.YNET = pyo.Var(
            dice.t, doc='Output net of damages equation (trillions 2005 USD per year)')
        dice.DAMAGES = pyo.Var(
            dice.t, doc='Damages (trillions 2005 USD per year)')
        dice.DAMFRAC = pyo.Var(
            dice.t, doc='Damages as fraction of gross output')
        dice.ABATECOST = pyo.Var(
            dice.t, doc='Cost of emissions reductions  (trillions 2005 USD per year)')
        dice.MCABATE = pyo.Var(
            dice.t, doc='Marginal cost of abatement (2005$ per ton CO2)')
        dice.CCA = pyo.Var(
            dice.t, doc='Cumulative industrial carbon emissions (GTC)')
        dice.CCATOT = pyo.Var(dice.t, doc='Total carbon emissions (GtC)')
        dice.PERIODU = pyo.Var(dice.t, doc='One period utility function')
        dice.CPRICE = pyo.Var(
            dice.t, doc='Carbon price (2005$ per ton of CO2)')
        dice.CEMUTOTPER = pyo.Var(dice.t, doc='Period utility')
        dice.UTILITY = pyo.Var(doc='Welfare function')
# ==========================================================================================================
        dice.CthawedPF = pyo.Var(
            dice.t,doc=' the amount of carbon in newly thawed permafrost(GtCO2)')  
        dice.CCumPF = pyo.Var(dice.t,doc='CO2em')
        dice.CO2_PF = pyo.Var(
            dice.t, doc=' the amount of carbondioxide from permafrost(GtCO2)')  
        dice.CH4_PF = pyo.Var(
            dice.t, doc=' the amount of methane from permafrost(GtCH4)')  
        dice.E_PF= pyo.Var(dice.t,doc='CO2 emissions (GtCO2 5 year)') 
# ==========================================================================================================
        dice.DDTpf = pyo.Var(dice.t, bounds=(0, +np.inf),doc='thawing index')
        dice.ALTpf = pyo.Var(dice.t,doc='active layer thickness') 
        dice.C_PF = pyo.Var(dice.t,doc='soil organic carbon in the active layer')         

        return

    def set_bounds(dice):

        for ttt in dice.t:

            dice.CCA[ttt].setub(dice.fosslim)

            dice.K[ttt].setlb(1)
            dice.MAT[ttt].setlb(10)
            dice.MU[ttt].setlb(100)
            dice.ML[ttt].setlb(1000)
            dice.C[ttt].setlb(2)
            dice.TATM[ttt].setub(20)
            dice.TOCEAN[ttt].setlb(-1)
            dice.TOCEAN[ttt].setub(20)
            dice.CPC[ttt].setlb(0.01)
            dice.TATM[ttt].setub(12)

            if ttt in dice.lag10:
                dice.S[ttt].fix(dice.optlrsav)

        dice.CCA[dice.t.first()].fix(400)
        dice.K[dice.t.first()].fix(dice.k0)
        dice.MAT[dice.t.first()].fix(dice.mat0)
        dice.MU[dice.t.first()].fix(dice.mu0)
        dice.ML[dice.t.first()].fix(dice.ml0)
        dice.TATM[dice.t.first()].fix(dice.tatm0)
        dice.TOCEAN[dice.t.first()].fix(dice.tocean0)
        dice.MIU[dice.t.first()].fix(dice.miu0)  
        dice.MIU[dice.t.last()].fix(0)      
        dice.CthawedPF[dice.t.first()].fix(0)
        dice.CCumPF[dice.t.first()].fix(0)
        dice.CO2_PF[dice.t.first()].fix(0)        
        dice.CH4_PF[dice.t.first()].fix(0) 
        dice.E_PF[dice.t.first()].fix(0)    
        return

    # ==========================================================================================================

    def run_DICE(dice, ifopt=0, n_repeat=1, tee=False):

        dice.ifopt = ifopt

        dice.solver_status = 'N'

        if pyo.value(dice.ifopt) == 0:

            print('Baseline Run:')

            # For base run, this subroutine calculates Hotelling rents
            # Carbon price is maximum of Hotelling rent or baseline price
            # The cprice equation is different from 2013R. Not sure what went wrong.

            dice.a2 = 0

            for i in np.arange(n_repeat):
                solution = dice.opt.solve(dice, tee=tee)

            if (solution.solver.status == pyo.SolverStatus.ok):
                print('Optimal solution. Objective Value: ', dice.OBJ())
                dice.solver_status = 'Y'
            else:
                raise Exception("Optimal solution not found!!")

            dice.a2 = pyo.value(dice.a20)

            for ttt in dice.t:
                dice.photel[ttt] = dice.CPRICE[ttt]
                if pyo.value(ttt) < pyo.value(dice.tnopol+1):
                    dice.CPRICE[ttt].setub(
                        max(pyo.value(dice.photel[ttt]), pyo.value(dice.cpricebase[ttt])))

        else:

            print('Optimization Run:')

            for ttt in dice.t:
                dice.CPRICE[ttt].setub(np.inf)
                dice.CPRICE[ttt].setlb(-np.inf)
            dice.MIU[dice.t.first()].fix(dice.miu0)

        for i in np.arange(n_repeat):
            solution = dice.opt.solve(dice, tee=tee)

        if (solution.solver.status == pyo.SolverStatus.ok):
            print('Optimal solution. Objective Value: ', dice.OBJ())
            dice.solver_status = 'Y'
        else:
            raise Exception("Optimal solution not found!!")

        for ttt in dice.t:
            dice.ppm[ttt] = dice.MAT[ttt] / 2.13
            dice.atfrac[ttt] = (dice.MAT[ttt] - dice.mateq) /                 (dice.CCATOT[ttt]+1E-6)
            dice.scc[ttt] = dice.dual[dice.EEQ[ttt]] /                 (dice.dual[dice.CC[ttt]]+1E-5)*-1000

        return solution


# In[ ]:





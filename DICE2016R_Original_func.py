#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pyomo.environ as pyo
from DICE2016R_Core_func import DICE2016R_Core

class DICE2016R_Original:

    def run_DICE(dice, ifopt=0, n_repeat=1, tee=False):

        solution = DICE2016R_Core.run_DICE(dice, ifopt, n_repeat, tee)

        return

    def set_bounds(dice):

        DICE2016R_Core.set_bounds(dice)

        return

    def define_constraints(dice):

        dice.EEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqEEQ, doc='Emissions equation')
        dice.MIUEQ = pyo.Constraint(dice.t, rule = DICE2016R_Core.eqMIUEQ)
        dice.EINDEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqEINDEQ, doc='Industrial emissions')
        dice.CCACCA = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqCCACCA, doc='Cumulative industrial carbon emissions')
        dice.CCATOTEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqCCATOTEQ, doc='Cumulative total carbon emissions')
        dice.FORCEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqFORCEQ, doc='Radiative forcing equation')
        dice.DAMFRACEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqDAMFRACEQ, doc='Equation for damage fraction')
        dice.DAMEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqDAMEQ, doc='Damage equation')
        dice.ABATEEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqABATEEQ, doc='Cost of emissions reductions equation')
        dice.MCABATEEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqMCABATEEQ, doc='Equation for MC abatement')
        dice.CARBPRICEEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqCARBPRICEEQ, doc='Carbon price equation from abatement')

        # Climate and carbon cycle =============

        dice.MMAT = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqMMAT, doc='Atmospheric concentration equation')
        dice.MML = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqMML, doc='Lower ocean concentration')
        dice.MMU = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqMMU, doc='Shallow ocean concentration')
        dice.TATMEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqTATMEQ, doc='Temperature-climate equation for atmosphere')
        dice.TOCEANEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqTOCEANEQ, doc='Temperature-climate equation for lower oceans')

        # Economic variables ===================

        dice.YGROSSEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqYGROSSEQ, doc='Output gross equation')
        dice.YNETEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqYNETEQ, doc='Output net of damages equation')
        dice.YY = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqYY, doc='Output net equation')
        dice.CC = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqCC, doc='Consumption equation')
        dice.CPCE = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqCPCE, doc='Per capita consumption definition')
        dice.SEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqSEQ, doc='Savings rate equation')
        dice.KK = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqKK, doc='Capital balance equation')
        dice.RIEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqRIEQ, doc='Interest rate equation')

        # Utility ================

        dice.CEMUTOTPEREQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqCEMUTOTPEREQ, doc='Period utility')
        dice.PERIODUEQ = pyo.Constraint(
            dice.t, rule=DICE2016R_Core.eqPERIODUEQ, doc='Instantaneous utility function equation')
        dice.UTIL = pyo.Constraint(rule=DICE2016R_Core.eqUTIL)
        dice.OBJ = pyo.Objective(rule=DICE2016R_Core.eqOBJ, sense=pyo.maximize)
#------------------------------------------------------        
        dice.CthawedPFEQ = pyo.Constraint(
            dice.t, rule = DICE2016R_Core.eqCthawedPFEQ, doc='thawed permafrost emissions')
        
        dice.CCumPFEQ = pyo.Constraint(
            dice.t, rule = DICE2016R_Core.eqCCumPFEQ,doc='CO2em')
        
        dice.CO2_PFEQ = pyo.Constraint(
            dice.t, rule = DICE2016R_Core.eqCO2_PFEQ, doc=' the amount of carbon in newly release permafrost(GtCO2)') 
    
        dice.CH4_PFEQ = pyo.Constraint(
            dice.t, rule = DICE2016R_Core.eqCH4_PFEQ, doc=' the amount of methane from permafrost(TgCH4)')      
        
        dice.E_PFEQ = pyo.Constraint(dice.t,rule = DICE2016R_Core.eqE_PFEQ, doc='CO2 emissions (GtCO2 5 year)') 

#------------------------------------------------------        
        dice.DDTpfEQ = pyo.Constraint(dice.t,rule=DICE2016R_Core.eqDDTpfEQ, doc='thawing index')
        dice.ALTpfEQ = pyo.Constraint(dice.t,rule=DICE2016R_Core.eqALTpfEQ, doc='active layer thickness')
        dice.C_PFEQ = pyo.Constraint(dice.t,rule=DICE2016R_Core.eqC_PFEQ, doc='soil organic carbon in the active layer') 
  
    
        return

    def initialize_variables(dice):

        DICE2016R_Core.initialize_variables(dice)

        return

    def precalculate_parameters(dice):

        DICE2016R_Core.precalculate_parameters(dice)

        return

    def initialize_parameters(dice):

        DICE2016R_Core.initialize_parameters(dice)

        return


# In[ ]:





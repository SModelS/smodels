{
    # channels: values (as python lists) are the pdg codes of the particles X,Y 
    # to be produced in the 2->2 processes p,p -> X,Y
    # names of keys are used for debugging only
    "channels" : {"1" : [1000024,-1000024]
    },
    # ------ MODE PARAMETER IS NOT USED IN THE CURRENT RESUMMINO INTERFACE VERSION ------ 
    # mode: if you want to calculate all channels at all orders, use 'all', 
    # if you want to check the cross section at LO before calculating NLO, use 'check'.
    ### Currently not used; if you want to always calculate the higher orders, put xsec_limit to 0.
    #"mode" : "check", 
    # -------------------------------------
    # The limit for the NLO calculation is determined by the 'xsec_limit' variable.
    # If the LO cross secion is below this value (in pb), no NLO (or NLO+NLL) cross-sections are
    # calculated.
    "xsec_limit" : 0.00001,
    # pdfs (case-sensitive): our default is cteq66, which gives results close 
    # to the PDFLHC2021_40 set. Any pdf can be used, but be sure to check
    # the name on the lhapdf website before putting any value here.
}

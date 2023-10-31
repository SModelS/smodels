{
    # channels: values (as python lists) are the pdg codes of the mother
    # particles to consider for productions
    # names of keys are used for debugging only
    #Defaut channels here are for one of the Bino-like scenarios, with M_1 < M_2, mu
    "channels" : {"1" : [1000023,1000024],"2" : [1000023, -1000024], "3" : [1000024,-1000024]
    },
    # ------------------------------------ MODE PARAMETERS IS NOT USED IN THE CURRENT RESUMMINO-INTERFACE VERSION
    # mode: If you want to calculate all the channels you choose at every order,
    # use 'all', if you want to check at LO before calculating NLO, use 'check'
    # instead. 
    ### For now, we choose to check by default, if you want to calculate all, put xsec_limit to 0.
    #"mode" : "check", 
    # -------------------------------------
    # The limit for the NLO calculation is determined by the 'xsec_limit' variable, 
    # below this value (in pb), no NLO cross-section are
    # calculated.
    "xsec_limit" : 0.00001,
    # pdfs (case-sensitive): use this to change the parton distribution
    # function used for the cross-section calculation.
    # PDFLHC2021 is default, but any pdf can be used. Be sure to check
    # the name on the lhapdf website before putting any value here.
    "pdfs" : {
        "pdf_lo" : "cteq66",
        "pdfset_lo" : 0,
        "pdf_nlo" : "cteq66",
        "pdfset_nlo" : 0
    }
}

{
    # channels: values (as python lists) are the pdg codes of the particles X,Y 
    # to be produced in the 2->2 processes p,p -> X,Y
    # names of keys are used for debugging only
    "channels" : {"1" : [1000023,1000024], "2" : [1000023, -1000024], 
                  "3" : [1000024,-1000024]
    },
    # -------------------------------------
    # The limit for the NLO calculation is determined by the 'xsec_limit' variable.
    # If the LO cross secion is below this value (in pb), no NLO (or NLO+NLL) cross-sections are
    # calculated.
    "xsec_limit" : 0.00001,
    # pdfs (case-sensitive): our default is cteq66, which gives results close 
    # to the PDFLHC2021_40 set. Any pdf can be used, but be sure to check
    # the name on the lhapdf website before putting any value here.
    "pdfs" : {
        "pdf_lo" : "cteq66",
        "pdfset_lo" : 0,
        "pdf_nlo" : "cteq66",
        "pdfset_nlo" : 0
    }
}

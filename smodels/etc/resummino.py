{
    # channels: values (as python lists) are the pdg codes of the mother
    # particles to consider for productions
    # names of keys are used for debugging only
    "channels" : {"1" : [1000023,1000025], "2" : [1000024,-1000024],
                "3" : [1000023,1000023]
    },
    # mode: If you want to calculate all the channel you chose at every order,
    # use 'all', if you want to check at LO before calculating NLO, use 'check'
    # instead. The limit for the NLO calculation is determined by the
    # 'mode_limit' variable, below this value (in pb), no NLO cross-section are
    # calculated.
    "mode" : "check",
    "mode_limit" : 0.00001,
    # pdfs: use this to change the parton distribution function used for the
    # cross-section calculation.
    # We used PDFLHC2021 by default, but any pdf can be used. Be sure to check
    # the name on the lhapdf website before putting any value here.
    "pdfs" : {
        "pdf_lo" : "cteq66",
        "pdfset_lo" : 0,
        "pdf_nlo" : "cteq66",
        "pdfset_nlo" : 0
    }
}

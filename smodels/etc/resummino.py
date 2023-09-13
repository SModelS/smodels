#channels: Value (in square bracket) represent the pdg code of the daugther particles to consider,
#key is not important here in general, you can choose 1,2,3,4,5, etc. by default

#mode: If you want to calculate all the channel you chose at every order, use 'all', if you want
#to check at LO before calculating NLO, use 'check' instead. The limit for the NLO calculation is determined 
#by the 'mode_limit' variable, below this value (in pb), no NLO cross-section are calculated.

#pdfs: You can use the pdfs part to change the parton distribution function used for the cross-section calculation.
#We used PDFLHC2021 by default, but any pdf can be used. Be sure to check the name on the lhapdf website before putting
#any value here.

{
    "channels" : {"1" : [1000023,1000025],
                "2" : [1000024,-1000024],
                "3" : [1000023,1000023]
    },
    "mode" : "check",
    "mode_limit" : 0.00001,
    "pdfs" : {
        "pdf_lo" : "cteq66",
        "pdfset_lo" : 0,
        "pdf_nlo" : "cteq66",
        "pdfset_nlo" : 0
    }
}




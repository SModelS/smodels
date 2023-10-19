#!/usr/bin/env python3

import numpy as np

def moment3(sigHi, sigLo):
    # Compute 3rd moment of bkg pdf from the upper, lower uncertainties
    from scipy import stats as st
    from scipy import integrate as integr
    def x3_times_bifurgauss(x, sigHi, sigLo):
        if x >= 0:
            return x*x*x * 2*sigHi/(sigLo + sigHi) * st.norm.pdf(x, 0, sigHi)
        else:
            return x*x*x * 2*sigLo/(sigLo + sigHi) * st.norm.pdf(x, 0, sigLo)
    m3, err = integr.quad(x3_times_bifurgauss, -np.inf, np.inf, args = (sigHi, sigLo))
    return m3

def parseCsvFile ( ):
    f=open("Results.csv")
    lines=f.readlines()
    f.close()
    numbers = []
    for line in lines:
        if line[0]=="#":
            continue
        if "fit" in line:
            numbers = []
        if line.startswith("Bin"):
            continue
        line = line.strip()
        if len(line)==0:
            continue
        tokens = list ( map ( float, line.split(",") ) )

        if len ( tokens ) == 4:
            idx = int ( tokens[0] )
            third = moment3 ( tokens[2], -tokens[3] )
            numbers.append ( third )
            # print ( f"{idx}, {third}" )
    print ( "numbers", numbers )

if __name__ == "__main__":
    parseCsvFile()
    print ( moment3( 13.561, 12.734 ) )
    # 149.74
    print ( moment3( 8.8574, 8.4841 ) )


#!/usr/bin/env python3

def convert ( frm, to ):
    f=open( frm )
    lines=f.readlines()
    f.close()
    f=open( to, "wt" )
    for line in lines:
        newline = line.strip()
        if newline.startswith("#") or newline.startswith("M"):
            f.write ( newline+"\n" )
        else:
            cols = newline.split(",")
            if cols == [""]:
                continue
            tokens = list(map(float,cols))
            mgl,mN=tokens[0],60
            x=tokens[1]
            mC=(mgl-mN)*x+mN
            if len(tokens)>2:
                ul=tokens[2]
                f.write ( "%.2f,%.2f,%.2f\n" % ( mgl,mC,ul ) )
            else:
                f.write ( "%.2f,%.2f\n" % ( mgl,mC ) )
    f.close()
    
convert("HEPData-ins1839446-v1-Upperlimit_GGridx.csv","ulT5WW60.csv")
convert("HEPData-ins1839446-v1-Exclusion_contour_2_(obs.).csv","obsT5WW60.csv")
convert("HEPData-ins1839446-v1-Exclusion_contour_2_(exp.).csv","expT5WW60.csv")
#convert("HEPData-ins1839446-v1-Upperlimit_GGx12.csv","ulT5WWp5.csv")
#convert("HEPData-ins1839446-v1-Exclusion_contour_1_(obs.).csv","obsT5WWp5.csv")
#convert("HEPData-ins1839446-v1-Exclusion_contour_1_(exp.).csv","expT5WWp5.csv")
convert("HEPData-ins1839446-v1-Upperlimit_SSgridx.csv","ulT6WW60.csv")
convert("HEPData-ins1839446-v1-Exclusion_contour_5_(obs.).csv","obsT6WW60.csv")
convert("HEPData-ins1839446-v1-Exclusion_contour_5_(exp.).csv","expT6WW60.csv")

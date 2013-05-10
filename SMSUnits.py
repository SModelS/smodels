from unum import Unum

# description see
# http://home.scarlet.be/be052320/Unum.html
# http://home.scarlet.be/be052320/docs.html

# install it with:
# sudo pip install ez_setup
# sudo easy_install Unum
# or sth similar

Unum.VALUE_FORMAT = "%0.2E"

fb=Unum.unit('fb')
pb=Unum.unit('pb', 1000 * fb)

eV=Unum.unit('eV')
keV=Unum.unit('keV',10**3*eV)
MeV=Unum.unit('MeV',10**6*eV)
GeV=Unum.unit('GeV',10**9*eV)
TeV=Unum.unit('TeV',10**12*eV)

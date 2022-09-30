#!/usr/bin/env python3

def find():
    try:
        import sysconfig
        ldlib = sysconfig.get_config_vars()["LDLIBRARY"]
        ldlib = ldlib.replace("lib","").replace(".so","").replace(".dll","")
        print ( f"-l{ldlib}" )
    except Exception as e:
        import sys
        print ( f"-lpython{sys.version_info.major}.{sys.version_info.minor}" )

find()

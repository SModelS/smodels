#!/usr/bin/env python3

def findLibFile():
    try:
        import sysconfig
        ldlib = sysconfig.get_config_vars()["LDLIBRARY"]
        ldlib = ldlib.replace("lib","").replace(".so","").replace(".dll","")
        ldlib = ldlib.replace(".dy","").replace(".a","")
        print ( f"-l{ldlib}" )
    except Exception as e:
        import sys
        print ( f"-lpython{sys.version_info.major}.{sys.version_info.minor}" )

if __name__ == "__main__":
    findLibFile()

def configure():
    import os
    import logging.config
    import set_path
    basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    logfile= '%s/etc/logging.conf' % basepath
    if not os.path.exists ( logfile ):
        # print "[tools/__init__.py] %s does not exist." % logfile
        import inspect, set_path
        base=os.path.realpath ( inspect.getabsfile(configure) )
        base=base.replace("tools/loggingConfiguration.py","")
        logfile=base+"etc/logging.conf"
        # print "base=",base
    logging.config.fileConfig( logfile, disable_existing_loggers=False)

def createLogger ( name ):
    configure()
    import logging
    return logging.getLogger(name)

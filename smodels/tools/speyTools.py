#!/usr/bin/env python3

"""
.. module:: speyTools
   :synopsis: a module that contains tools and convenience methods 
              that we use in connection with spey.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""


def getSpeyInitialisation ( dataset, allowNegativeSignals : bool = False,
                            initial_bracket : bool = True ):
    """ get decent initial bounds and initial values for a statModel
    :param allowNegativeSignals: if true, then bound the poi to positive values
    :param initial_bracket: also supply an initial bracket for CLs-alpha root
                            finding?
    """
    if dataset.getType()=="efficiencyMap":
        ini = getInitialisationForSingleRegions ( dataset, allowNegativeSignals )
        return filterInitialBracket ( ini, initial_bracket )
        
    if dataset.type=="simplified":
        ini = getInitialisationForSL ( dataset, allowNegativeSignals )
        return filterInitialBracket ( ini, initial_bracket )
    if dataset.type!="pyhf":
        raise Exception ( f"dont know that dataset type {dataset.type}" )
    ini = getInitialisationForPyhf ( dataset, allowNegativeSignals )
    return filterInitialBracket ( ini, initial_bracket )

def filterInitialBracket ( args : tuple, initial_bracket : bool ):
    """ possibly filter out the suggestion for an initial bracket """
    if initial_bracket:
        return dictionary
    args[2].pop ( "low_init", None )
    args[2].pop ( "hig_init", None )
    return args

def getInitialisationForSingleRegions ( dataset, allowNegativeSignals : bool = False ):
    statModel = dataset.statModel
    model = statModel.backend.model
    config = model.config()
    init = config.suggested_init
    bounds = config.suggested_bounds
    args = {}
    return init,bounds,args

def getInitialisationForPyhf ( dataset, allowNegativeSignals : bool = False ):
    """ get decent initial bounds and initial values for a statModel
    :param allowNegativeSignals: if true, then bound the poi to positive values
    """
    statModel = dataset.statModel
    model = statModel.backend.model
    config = model.config()
    init = config.suggested_init
    bounds = config.suggested_bounds
    """
    print ( "datasets", [ x.dataInfo.dataId for x in dataset._datasets ] )
    bg = [ x.dataInfo.expectedBG for x in dataset._datasets ]
    signal = model.signal[0]["value"]["data"]
    print ( "signals", model.signal[0]["value"]["data"] )
    print ( "nobs", model.background["channels"][0]["samples"][0]["data"] )
    print ( "suggested bounds", bounds )
    if allowNegativeSignals:
        bounds[config.poi_index] = (config.minimum_poi, 100)
    else:
        bounds[config.poi_index] = (0, 100)
    mui = []
    for i in range(len(nobs)):
        mui.append ( (nobs[i]-bg[i])/signal[i] )
    print ( "mui", mui )
    bounds = [(suggested[0]-800,suggested[1]+800) for suggested in config.suggested_bounds]
    bounds[config.poi_index] = (-200, 800) # for now!
    """
    # args = { "maxiter": 500, "method": "SLSQP" } ## extra args for the optimizers
    # method BFGS, SLSQP
    #args = { "maxiter": 500, "method": "BFGS", "ntrials": 3,
    #         "xrtol": 1e-6 "low_init": bounds[0][0], 
#                "hig_init": bounds[0][1] 
    #}
    # args = { "maxiter": 500, "ntrials": 1, "method": "SLSQP" }
    args = { "maxiter": 500, "ntrials": 1, "method": None }
    return init,bounds,args

def alternateMethod ( args : dict ):
    """ try a different method. i hope we wont need this in the long run """
    if not "method" in args or args["method"] is None:
        args["method"]="BFGS"
        return
    if args["method"]=="SLSQP":
        args["method"]="BFGS"
        return
    if args["method"]=="BFGS":
        args["method"]="SLSQP"
        return
    args["method"]="L-BFGS-B"

def getInitialisationForSL ( dataset, allowNegativeSignals : bool = False ):
    """ get decent initial bounds and initial values for an SL statModel
    :param allowNegativeSignals: if true, then bound the poi to positive values
    """
    import numpy as np
    statModel = dataset.statModel
    config = statModel.backend.model.config()
    init = config.suggested_init
    bounds = config.suggested_bounds
    args = {} ## extra args for the optimizers
    if False: ## these were the old values!
        bounds = [(suggested[0]-800,suggested[1]+800) for suggested in config.suggested_bounds]
        bounds[config.poi_index] = (-520, 800) # for now!
        return bounds,init,args
    assert config.poi_index == 0, f"Error: I assume the poi index to be zero, not {config.poi_index}"
    init[1:] = statModel.backend.model.observed - statModel.backend.model.background - statModel.backend.model.signal
    numerator, denominator = [], []
    for i in range ( len ( statModel.backend.model.observed ) ):
        cov = statModel.backend.model.covariance[i][i]
        x = np.sqrt ( cov )
        bounds[i+1]=(-5*x+init[i+1],5*x+init[i+1])
        sig = statModel.backend.model.signal[i]
        obs = statModel.backend.model.observed[i]
        bg = statModel.backend.model.background[i]
        if sig > 0. and cov > 0.:
            # for the given region, mui would be the best bet for muhat
            mui = ( obs - bg ) / sig
            # its variance is this
            cov_mui =  ( obs + cov ) / ( sig**2 )
            # the inverse of which will be our weights
            wi = 1. / cov_mui
            numerator.append ( wi * mui )
            denominator.append ( wi )
    ## ok so the bounds should be -5*x,5*x with x being np.sqrt(statModel.backend.model.covariance[i][i], the initial values just the diff between observation and expectation
    init_muhat = np.sum ( numerator ) / np.sum ( denominator)
    if not allowNegativeSignals and init_muhat<0.:
        init_muhat = 0.
    err_muhat = np.sqrt ( len(denominator) / np.sum ( denominator ) )
    # err_muhat = np.sqrt ( np.sum ( totweight )**(-2) * np.sum ( covest ) )

    init[config.poi_index] = init_muhat
    minmu, maxmu = -5*err_muhat + init_muhat, 5*err_muhat + init_muhat
    if not allowNegativeSignals and minmu < 0.:
        minmu = 0.
    bounds [ config.poi_index ] = ( minmu, maxmu )
    if False:
        print ( f"residual ({len(init)}) is", init_muhat, "+-", err_muhat )
        print ( "errors are", covest )
        print ( "residuals are", residuals )
        print ( "init pars in srCombinations are", init[:3] )
        print ( "signals   in srCombinations are", statModel.backend.model.signal[:] )
        print ( "deltas  in srCombinations are", statModel.backend.model.observed - statModel.backend.model.background )
        print ( "bounds    in srCombinations are", bounds[:3] )
    args = { "maxiter": 500, "method": None, "ntrials": 1,
             "low_init": bounds[0][0], 
                "hig_init": bounds[0][1] }
    # args["method"]="BFGS"
    args["tol"]=1e-3
    # args["xrtol"]=1e-6
    # print ( f"speyTools: initbracket is", args["low_init"], args["hig_init"] )
    return init,bounds,args


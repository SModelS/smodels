#!/usr/bin/env python3

"""
.. module:: expSMSMapObj
   :synopsis: Contains the ExpSMSMap used to collect the SMS in the database and their
              mapping to the corresponding experimental results, datasets and txnames.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

class ExpSMSMap(object):
    """ 
    Holds a list of experimental results, their SMS (ExpSMS) and facilities
    for mapping the SMS to their corresponding results, datasets and txnames.
    """

    def __init__(self, expResultList=[]):
        """
        :param exptResultList: List of ExptResult objects used to build the map
        """

        self.expResults = expResultList[:]
        self.smsDict = {}
        self.smsExpMap = {}
        self.canonNames = {}
                
        if len(self.expResults) > 0:
            self.getUniqueSMSDict()

        # Use smsDict to build a dictionary for
        # mapping canonical names to the list of SMS (indices)
        for isms, sms in self.smsDict.items():
            cName = sms.canonName
            if cName in self.canonNames:
                self.canonNames[cName].append(isms)
            else:
                self.canonNames[cName] = [isms]

        # Create empty dictionary to store matches
        # used for theory predictions
        self.smsMatches = {isms : [] for isms in self.smsDict}   

    def __len__(self):
        """
        Returns the number of unique SMS
        """

        return len(self.smsDict)

    def getUniqueSMSDict(self):
        """
        Iterates over all (active) experimental results and build two dictionaries:
        one mapping indices to the unique SMS and another one
        a with the unique SMS indices as keys and a dictionary
        {ExpResult_index : {DataSet_index : {TxName_index : smsLabel}}} as values.

        :return: Index->SMS dictionary and SMS_index -> Exp_index dictionary
        """

        self.smsDict = {}
        self.smsExpMap = {}
                
        # Loop over active experimental results:
        # (Note that expResultList returns only the active results)
        for iexp,exp in enumerate(self.expResults):
            # Loop over datasets:
            for ids,dataset in enumerate(exp.datasets):
                # Loop over txnames:
                for itx,tx in enumerate(dataset.txnameList):
                    # Sort TxName (if it is already tagged as sorted, do nothing)
                    tx.sortSMSMap()
                    # Loop over ExpSMS in the txname:
                    for sms,smsLabel in tx.smsMap.items():
                        ismsMatch = None
                        # Loop for an identical SMS in smsDict
                        # (there can only be one match, since the SMS within
                        # a given txname must be unique)
                        for isms,sms1 in self.smsDict.items():
                            if sms.identicalTo(sms1):
                                ismsMatch = isms
                                break
                                        
                        # Update dictionary
                        if ismsMatch is None:
                            isms = len(self.smsDict)
                            self.smsDict[isms] = sms
                            self.smsExpMap[isms] = {iexp : {ids : {itx : smsLabel}}}
                        else:
                            useDict = self.smsExpMap[ismsMatch]
                            entry = [iexp,ids]
                            # Loop over entries until one is not found
                            for key in entry:
                                if useDict.get(key) is None:
                                    useDict.update({key : {}})
                                useDict = useDict.get(key)
                            useDict[itx] = smsLabel

    def clearMatches(self):
        """
        Clear the matches dictionary
        """

        self.smsMatches = {isms : [] for isms in self.smsDict}   

    def getMatchesTo(self,sms):
        """
        Checks if any of the unique SMS in the map
        matches the sms.

        :param sms: A SMS object (TheorySMS)
        """

        # First filter the possibilities by canonical names
        # (inclusive canonical names are always selected)
        canonName = sms.canonName
        selectedSMSlist = []
        for cName,ismsList in self.canonNames.items():
            if cName == canonName:
                selectedSMSlist += ismsList
        
        # Loop over selected SMS and check for matches
        # Store the (correctly orderd) match in self.smsMatches
        for isms in selectedSMSlist:
            txsms = self.smsDict[isms]
            matchedSMS = txsms.matchesTo(sms)
            if matchedSMS is None:
                continue
            # Get all the experimental results containing the sms:
            expList = [self.expResults[iexp] for iexp in self.smsExpMap[isms]]
            expTypes = set([exp.globalInfo.type for exp in expList])
            # Tag the original SMS as covered by the exp types
            for t in expTypes:
                sms.setCoveredBy(t)
            effDict = {}
            for iexp in self.smsExpMap[isms]:
                exp = self.expResults[iexp]
                effDict[iexp] = {}
                for ids in self.smsExpMap[isms][iexp]:
                    dataset = exp.datasets[ids]
                    effDict[iexp][ids] = {}
                    for itx,smsLabel in self.smsExpMap[isms][iexp][ids].items():
                        tx = dataset.txnameList[itx]                        
                        eff = tx.getEfficiencyFor(matchedSMS)
                        if eff is None or abs(eff) < 1e-14:
                            continue
                        # Store the efficiency:
                        print('isms,iexp,ids,itx=',isms,iexp,ids,itx)
                        print('smsLabel=',smsLabel)
                        effDict[iexp][ids][itx] = {smsLabel : eff}
                        # Tag the original SMS as tested
                        sms.setTestedBy(exp.globalInfo.type)

            # Add the efficiency diction to the matchedSMS
            matchedSMS.effDict = effDict
            # Add the matches sms to the smsMatches dict
            self.smsMatches[isms].append(matchedSMS)


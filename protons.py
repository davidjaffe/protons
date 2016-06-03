#!/usr/bin/env python
'''
estimate proton content of materials
Data from NIST : http://www.nist.gov/pml/data/comp.cfm/
Official atomic weights, etc. :
Atomic weights of the elements 2013 (IUPAC Technical Report) Pure Appl. Chem. 2016, 88(3), 265-291
20160603
'''
import sys
import os

class protons():
    def __init__(self):# name     [Rel.atomic mass, unc]    [abundance, unc]    [lower, upper limit on atomic weight]
        self.isoInfo = {'1H'  : [ [1.007825032, 9.e-11],    [0.999885, 70.e-6], [1.00784, 1.00811] ],
                        '2H'  : [ [2.01410177812, 12.e-11], [0.000115, 70.e-6], ],
                        '12C' : [ [12.0,            0.]   , [0.9893,    8.e-4], [12.0096,12.0116]  ],
                        '13C' : [ [13.00335483507, 23.e-11],[0.0107,    8.e-4], ],
                        }
        return
    def massFraction(self,isotope,molecule):
        '''
        return mass fraction of isotope in molecule using central values and limits on atomic weight
        example: isotope = '2H', molecule = {'H':6, 'C':6} for benzene
        '''
        totalM,isoM = 0., 0.
        totMlo,totMhi = 0.,0.
        for element in molecule:
            n = molecule[element]
            # calculates total abundance of specified isotopes of element
            totalAbund = 0.
            for iso in self.isoInfo:
                if element in iso: totalAbund += self.isoInfo[iso][1][0]
            #print element,totalAbund
            for iso in self.isoInfo:
                if element in iso:
                    ram,dram = self.isoInfo[iso][0]   # relative atomic mass, unc
                    abund,dabund = self.isoInfo[iso][1]  # abundance, unc
                    if len(self.isoInfo[iso])>2:
                        awlo,awhi = self.isoInfo[iso][2]  # lower, upper limit on atomic weight
                        totMlo += float(n)*awlo
                        totMhi += float(n)*awhi
                    totalM += float(n)*abund*ram/totalAbund
                    if isotope==iso:
                        isoM = float(n)*abund*ram/totalAbund
        f,flo,fhi = -1.,-1.,-1.
        if totalM>0: f = isoM/totalM
        if totMlo>0:
            flo = isoM/totMlo
            fhi = isoM/totMhi
        return f,flo,fhi
    def moleculeName(self,molecule):
        '''
        string giving name of molecule
        '''
        name = ''
        for iso in molecule:
            name+= iso + '('+str(molecule[iso])+')'
        return name
    def report(self,words,cv,uv,lv):
        '''
        report words with  central value cv, upper and lower values uv,lv
        '''
        u = 1000.*(uv-cv) # per mil
        l = 1000.*(lv-cv)
        fu,fl = 0.,0.
        if cv!=0:
            fu = 100.*(uv-cv)/cv
            fl = 100.*(lv-cv)/cv
        print '{0} {1:f} absolute [{2:.3f},{3:.3f}] per mille, relative [{4:.3f},{5:.3f}] percent'.format(words,cv,l,u,fl,fu)
        return
if __name__ == '__main__' :
    P = protons()
    composition = {15:0.0, 16:0.0698, 17:0.306, 18:0.450, 19:0.174, 20:0.0}
    WF = {}
    for isotope in ['1H','2H']:
        hydrogen = 'H'
        print 'hydrogen is',hydrogen
        wf,wflo,wfhi = 0.,0.,0.
        for nC in sorted(composition):
            nH = 2*nC-6
            molecule = {hydrogen:nH, 'C':nC}
            f,flo,fhi = P.massFraction(isotope,molecule)
            name = P.moleculeName(molecule)
            #print name,f,'[',flo,',',fhi,']'
            P.report(name+' '+isotope+' fraction',f,flo,fhi)
            wf += f*composition[nC]
            wflo += flo*composition[nC]
            wfhi += fhi*composition[nC]
        #print 'weighted',isotope,'fraction',wf,'[',wflo,',',wfhi,']'
        P.report('weighted '+isotope+' fraction',wf,wflo,wfhi)
        WF[isotope] = wf
    wf = WF['1H']+WF['2H']
    print 'total 1H+2H weighted fraction',wf,'total 1H+2H fraction minus 1H fraction',wf-WF['1H']
    print '1H fraction relative to total 1H+2H',WF['1H']/wf,'1-rel',1.-WF['1H']/wf

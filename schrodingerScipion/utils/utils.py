# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import numpy as np
import os
from pyworkflow.utils.path import moveFile
import pyworkflow.object as pwobj

def putMol2Title(fn, title=""):
    i=0
    fhIn = open(fn)
    fhOut = open(fn+".aux",'w')
    for line in fhIn.readlines():
        if i!=1:
            fhOut.write(line)
        else:
            if title!="":
                fhOut.write(title+"\n")
            else:
                fhOut.write(os.path.splitext(os.path.split(fn)[1])[0]+"\n")
        i+=1
    fhIn.close()
    fhOut.close()
    moveFile(fn+".aux",fn)

def sortDockingResults(smallList):
    ds = []
    le = []
    leSA = []
    leLn = []
    for small in smallList:
        ds.append(small._energy.get())
        le.append(small.ligandEfficiency.get())
        leSA.append(small.ligandEfficiencySA.get())
        leLn.append(small.ligandEfficiencyLn.get())

    iN = 100.0 / len(ds)
    ds = np.asarray(ds)
    le = np.asarray(le)
    leSA = np.asarray(leSA)
    leLn = np.asarray(leLn)

    h = np.zeros(len(ds))
    i=0
    for small in smallList:
        hds = np.sum(ds >= small._energy.get()) * iN
        hle = np.sum(le >= small.ligandEfficiency.get()) * iN
        hleSA = np.sum(leSA >= small.ligandEfficiencySA.get()) * iN
        hleLn = np.sum(leLn >= small.ligandEfficiencyLn.get()) * iN
        hi = -0.25 * (hds + hle + hleSA + hleLn)
        small.Hrank = pwobj.Float(hi)
        h[i]=hi
        i+=1

    return np.argsort(h)

def parseLogPockets(fnLog):
    pocketsDic, pId = {}, 1
    with open(fnLog) as fh:
        for line in fh:
            if line.startswith("SiteScore"):
                # SiteScore size   Dscore  volume  exposure enclosure contact  phobic   philic   balance  don/acc
                pocketsDic[pId] = [float(x) for x in fh.readline().split()]
                pId += 1
    return pocketsDic

def parseLogProperties(fnLog, idx=None):
    pDic, pId = {}, 1
    with open(fnLog) as fh:
        for line in fh:
            if line.startswith("SiteScore"):
                # SiteScore size   Dscore  volume  exposure enclosure contact  phobic   philic   balance  don/acc
                if idx == pId:
                    keys = line.split()
                    values = [float(x) for x in fh.readline().split()]
                    pDic = dict(zip(keys, values))
                pId += 1
    return pDic

def getChargeFromMAE(maeFile):
    charge = 0
    with open(maeFile) as f:
        keys, values = False, False
        for line in f:
            if keys:
                kList.append(line.strip())
            if values and not line.strip().endswith('{') and not line.strip() == ':::':
                elements = maeLineSplit(line)
                newCharge = float(elements[chargeIdx])
                charge += newCharge

            if line.strip().endswith('{'):
                keys, values = True, False
                kList = []
            elif line.strip() == ':::':
                if keys:
                    keys = False
                    if 'i_m_formal_charge' in kList:
                        values = True
                        chargeIdx = kList.index('i_m_formal_charge')
                elif values:
                    values = False
    return charge

def maeLineSplit(maeLine):
    '''Some elements are strings surrounded by "" but with spaces in between,
        accounting for several elements of the list when they are actually just 1'''
    elements = []
    stri, ele = False, ''
    for a in maeLine.strip():
        if stri:
            ele += a
            if a == '"':
                stri = False
        else:
            if a == ' ':
                elements.append(ele)
                ele = ''
            else:
                ele += a
                if a == '"':
                    stri = True
    elements.append(ele)
    return elements
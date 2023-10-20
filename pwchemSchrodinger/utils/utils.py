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
# General imports
import numpy as np
import os, subprocess, time
from subprocess import check_call

# Scipion em imports
from pyworkflow.utils.path import moveFile
import pyworkflow.object as pwobj

# Scipion chem imports
from pwchem.objects import SmallMolecule

# Plugin imports
from ..constants import TIMESTEP, PRESSURE, BAROSTAT, BROWNIAN, TENSION, RESTRAINS, MSJ_SYSMD_SIM, MSJ_SYSMD_INIT

def putMol2Title(fn, title=""):
    with open(fn) as fhIn:
        with open(fn + ".aux", 'w') as fhOut:
            for i, line in enumerate(fhIn.readlines()):
                if i!=1:
                    fhOut.write(line)
                else:
                    if title!="":
                        fhOut.write(title+"\n")
                    else:
                        fhOut.write(os.path.splitext(os.path.split(fn)[1])[0]+"\n")
    moveFile(fn+".aux", fn)

def putSDFTitle(fn, title=''):
    with open(fn) as fhIn:
        with open(fn + ".aux", 'w') as fhOut:
            for i, line in enumerate(fhIn.readlines()):
                if i != 0:
                    fhOut.write(line)
                else:
                    if title != "":
                        fhOut.write(title + "\n")
                    else:
                        fhOut.write(os.path.splitext(os.path.split(fn)[1])[0] + "\n")
    moveFile(fn + ".aux", fn)

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

def unzipMaegz(inFile, outFile=None):
    if inFile.endswith('.mae'):
        return inFile
    elif inFile.endswith('.maegz'):
        if not outFile:
            outFile = inFile.replace('.maegz', '.mae')
        subprocess.check_call('zcat {} > {}'.format(inFile, outFile), shell=True)
        return outFile
    else:
        print('Format not recognized. File must be a .maegz')

def getChargeFromMAE(maeFile):
    maeFile = unzipMaegz(maeFile)

    charge = 0
    with open(maeFile) as f:
        keys, values = False, False
        for line in f:
            if keys:
                kList.append(line.strip())
            if values and not line.strip().endswith('{') and line.strip() != ':::':
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

def getNumberOfStructures(maeFile):
    maeFile = unzipMaegz(maeFile)
    mols = 0
    with open(maeFile) as f:
        for line in f:
            mols += line.count('f_m_ct')
    return mols

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

# ----------------------- Protocol utils -----------------------
def saveMolecule(protocols, molFn, molSet, oriMol):
    while protocols.saving:
        time.sleep(0.2)
    protocols.saving = True
    smallMolecule = SmallMolecule()
    smallMolecule.copy(oriMol, copyId=False)
    smallMolecule.setFileName(molFn)
    confId = getConfId(molFn, oriMol.getMolName())
    if confId:
        smallMolecule.setConfId(molFn.split('-')[-1].split('.')[0])

    molSet.append(smallMolecule.clone())
    protocols.saving = False

def getConfId(molFn, molName):
    try:
        return molFn.split(molName)[1].split('-')[1].split('.')[0]
    except Exception:
        return None

def createMSJDic(protocol):
    msjDic = {}
    for pName in protocol.getStageParamsDic(type='Normal').keys():
      if hasattr(protocol, pName):
        msjDic[pName] = getattr(protocol, pName).get()
      else:
        print('Something is wrong with parameter ', pName)

    for pName in protocol.getStageParamsDic(type='Enum').keys():
      if hasattr(protocol, pName):
        msjDic[pName] = protocol.getEnumText(pName)
      else:
        print('Something is wrong with parameter ', pName)
    return msjDic

def buildSimulateStr(protocol, msjDic):
    '''Checks the values stored in the msjDic and trnaslates them into msjStr.
        If a value is not found in the msjDic, the default is used'''
    msjDic = protocol.addDefaultForMissing(msjDic)

    glueArg = '[]'
    if msjDic['glueSolute']:
      glueArg = 'solute'

    # NearT and farT must be at least the boundT
    msjDic['nearT'] = max(msjDic['bondedT'], msjDic['nearT'])
    msjDic['farT'] = max(msjDic['bondedT'], msjDic['farT'])
    timeStepArg = TIMESTEP % (msjDic['bondedT'], msjDic['nearT'], msjDic['farT'])

    pressureArg, barostatArg = '', ''
    method = protocol._thermoDic[msjDic['thermostat']]
    ensemType = msjDic['ensemType']
    if ensemType not in ['NVE', 'NVT', 'Minimization (Brownian)']:
      pressureArg = PRESSURE % (msjDic['pressure'], msjDic['coupleStyle'].lower())
      barostatArg = BAROSTAT % (msjDic['presMDCons'])
      method = protocol._baroDic[msjDic['barostat']]

    tensionArg, brownianArg = '', ''
    if ensemType == 'Minimization (Brownian)':
      ensemType = 'NVT'
      method = 'Brownie'
      brownianArg = BROWNIAN % (msjDic['deltaMax'])
    elif ensemType == 'NPgT':
      tensionArg = TENSION % msjDic['surfTension']

    restrainArg = ''
    if msjDic['restrains'] != 'None':
      restrainArg = RESTRAINS % (msjDic['restrains'].lower(), msjDic['restrainForce'])

    if not msjDic['annealing']:
      annealArg = 'off'
      tempArg = msjDic['temperature']
    else:
      annealArg = 'on'
      tempArg = protocol.parseAnnealing(msjDic['annealTemps'])

    msjStr = MSJ_SYSMD_SIM % (annealArg, os.path.abspath(protocol._getTmpPath()),
                               glueArg, msjDic['simTime'], timeStepArg, tempArg, pressureArg,
                               tensionArg, ensemType, method, msjDic['tempMDCons'], barostatArg, brownianArg,
                               restrainArg, msjDic['velResamp'], msjDic['trajInterval'])
    return msjStr

def getJobName(protocol):
    files = os.listdir(protocol._getTmpPath())
    for f in files:
        if f.endswith('.msj'):
            return f.replace('.msj', '')

def getSchJobId(protocol, jobControlProg):
    jobId = None
    jobListFile = os.path.abspath(protocol._getTmpPath('jobList.txt'))
    if getJobName(protocol):
        check_call(jobControlProg + ' -list {} | grep {} > {}'.
                    format(getJobName(protocol), getJobName(protocol), jobListFile), shell=True)
        with open(jobListFile) as f:
            jobId = f.read().split('\n')[0].split()[0]
    return jobId

def setAborted(protocol, jobControlProg):
    jobId = getSchJobId(protocol, jobControlProg)
    if jobId:
        print('Killing job: {} with jobName {}'.format(jobId, getJobName(protocol)))
        check_call(jobControlProg + ' -kill {}'.format(jobId), shell=True)

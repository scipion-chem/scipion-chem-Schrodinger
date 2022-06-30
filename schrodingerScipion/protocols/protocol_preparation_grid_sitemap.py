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
import os, shutil, glob
from subprocess import CalledProcessError

from pyworkflow.protocol.params import MultiPointerParam, STEPS_PARALLEL, PointerParam, BooleanParam, \
    FloatParam, IntParam, EnumParam
from pyworkflow.object import String, Float
from pyworkflow.utils.path import createLink, makePath

from pwem.protocols import EMProtocol
from schrodingerScipion.objects import SchrodingerGrid, SetOfSchrodingerGrids

from schrodingerScipion import Plugin as schrodinger_plugin

structConvertProg = schrodinger_plugin.getHome('utilities/structconvert')

class ProtSchrodingerGridSiteMap(EMProtocol):
    """Calls glide to prepare a grid with an input that is a set of binding sites"""
    _label = 'grid definition binding site (glide)'
    _program = ""

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructROIs', PointerParam, pointerClass="SetOfStructROIs",
                      label='Sets of Structural ROIs:', allowsNull=False,
                      help='Sets of known or predicted protein structural ROIs to center the grid on')

        group = form.addGroup('Inner box')
        group.addParam('innerAction', EnumParam, default=0, label='Determine inner box: ',
                      choices=['Manually', 'PocketDiameter'], display=EnumParam.DISPLAY_HLIST,
                      help='How to set the inner box.'
                           'Manually: you will manually set the same x,y,z for every ROI'
                           'PocketDiameter: the diameter * n of each ROI will be used. You can set n')

        line = group.addLine('Inner box (Angstroms)', condition='innerAction==0',
                             help='The docked ligand mass center must be inside the inner box radius')
        line.addParam('innerX', IntParam, default=10, label='X')
        line.addParam('innerY', IntParam, default=10, label='Y')
        line.addParam('innerZ', IntParam, default=10, label='Z')

        group.addParam('diameterNin', FloatParam, default=0.8, condition='innerAction==1',
                       label='Size of inner box vs diameter: ',
                       help='The diameter * n of each ROI will be used as inner box side')

        group = form.addGroup('Outer box')
        group.addParam('outerAction', EnumParam, default=0, label='Determine outer box: ',
                       choices=['Manually', 'PocketDiameter'], display=EnumParam.DISPLAY_HLIST,
                       help='How to set the outer box.'
                            'Manually: you will manually set the same x,y,z for every structural ROI'
                            'PocketDiameter: the diameter * n of each pocket will be used. You can set n')

        line = group.addLine('Outer box (Angstroms)', condition='outerAction==0',
                            help='The docked ligand atoms must be inside the outer box radius.')
        line.addParam('outerX', IntParam, default=30, label='X')
        line.addParam('outerY', IntParam, default=30, label='Y')
        line.addParam('outerZ', IntParam, default=30, label='Z')

        group.addParam('diameterNout', FloatParam, default=1.2, condition='outerAction==1',
                       label='Size of outer box vs diameter: ',
                       help='The diameter * n of each structural ROI will be used as outer box side')

        group = form.addGroup('Hydrogen bonds')
        group.addParam('HbondDonorAromH', BooleanParam, default=False, label='Aromatic H as H-bond donors:',
                      help='Accept aromatic hydrogens as potential H-bond donors.')
        group.addParam('HbondDonorAromHCharge', FloatParam, default=0.0, label='Aromatic H as H-bond donors Charge:',
                      condition='HbondDonorAromH',
                      help='Partial charge cutoff for accepting aromatic hydrogens as potential H-bond donors. '
                           'The cutoff is applied to the actual (signed) charge, not the absolute value.')
        group.addParam('HbondAcceptHalo', BooleanParam, default=False, label='Halogens as H-bond acceptors:',
                      help='Accept halogens (neutral or charged, F, Cl, Br, or I) as H-bond acceptors.')
        group.addParam('HbondDonorHalo', BooleanParam, default=False, label='Halogens as H-bond donors:',
                      help='Accept the halogens (Cl, Br, I, but not F) as potential H-bond '
                           '(noncovalent interaction) donors')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        cStep = self._insertFunctionStep('convertStep', prerequisites=[])
        prepSteps = []
        for pocket in self.inputStructROIs.get():
            pStep = self._insertFunctionStep('preparationStep', pocket.clone(), prerequisites=[cStep])
            prepSteps.append(pStep)

        self._insertFunctionStep('createOutputStep', prerequisites=prepSteps)

    def convertStep(self):
        maeFile = self.inputStructROIs.get().getMAEFile()
        if not maeFile:
            pdbFile = self.inputStructROIs.get().getProteinFile()
            maeFile = self._getExtraPath('inputReceptor.maegz')
            prog = schrodinger_plugin.getHome('utilities/prepwizard')
            args = ' -WAIT -noprotassign -noimpref -noepik {} {}'. \
                format(os.path.abspath(pdbFile), os.path.abspath(maeFile))
            self.runJob(prog, args, cwd=self._getExtraPath())
        else:
            ext = os.path.splitext(maeFile)[1]
            shutil.copy(maeFile, self._getExtraPath('inputReceptor{}'.format(ext)))

    def preparationStep(self, pocket):
        x, y, z = pocket.calculateMassCenter()

        fnGridDir = self.getGridDir(pocket.getObjId())
        gridName = self.getGridName(pocket.getObjId())
        makePath(fnGridDir)

        fnJob = os.path.abspath(os.path.join(fnGridDir, gridName)) + '.inp'
        fh = open(fnJob, 'w')
        fh.write("GRIDFILE %s.zip\n" % gridName)
        fh.write("OUTPUTDIR %s\n" % fnGridDir)
        fh.write("RECEP_FILE %s\n" % os.path.abspath(self.getInputMaeFile()))
        fh.write("REC_MAECHARGES True\n")
        fh.write("HBOND_DONOR_AROMH %s\n" % self.HbondDonorAromH.get())
        if self.HbondDonorAromH.get():
            fh.write("HBOND_DONOR_AROMH_CHARGE %f\n" % self.HbondDonorAromHCharge.get())
        fh.write("HBOND_ACCEP_HALO %s\n" % self.HbondAcceptHalo.get())
        fh.write("HBOND_DONOR_HALO %s\n" % self.HbondDonorHalo.get())
        fh.write("INNERBOX %d,%d,%d\n" % (self.getInnerBox(pocket)))
        fh.write("ACTXRANGE %d\n" % self.getOuterBox(pocket)[0])
        fh.write("ACTYRANGE %d\n" % self.getOuterBox(pocket)[1])
        fh.write("ACTZRANGE %d\n" % self.getOuterBox(pocket)[2])
        fh.write("OUTERBOX %d,%d,%d\n" % (self.getOuterBox(pocket)))
        fh.write("GRID_CENTER %s,%s,%s\n" % (x, y, z))
        fh.close()

        args = "-WAIT -LOCAL %s.inp" % (gridName)
        self.runJob(schrodinger_plugin.getHome('glide'), args, cwd=fnGridDir)

    def createOutputStep(self):
        outGrids = SetOfSchrodingerGrids(filename=self._getPath('SchGrids.sqlite'))
        for pocket in self.inputStructROIs.get():
            gridId = pocket.getObjId()
            fnDir = self.getGridDir(gridId)
            if os.path.exists(fnDir):
                fnBase = os.path.split(fnDir)[1]
                fnGrid = os.path.join(fnDir, '%s.zip' % fnBase)
                if os.path.exists(fnGrid):
                    SchGrid = SchrodingerGrid(filename=fnGrid, **self.getGridArgs(pocket))
                    SchGrid.structureFile = String(self.getInputMaeFile())
                    if str(pocket.getScore()) != 'None':
                        SchGrid.pocketScore = Float(pocket.getScore())
                    SchGrid.setObjId(gridId)
                    # gridFile.bindingSiteDScore = Float(dscore)
                    outGrids.append(SchGrid)

        outGrids.buildBBoxesPML()
        self._defineOutputs(outputGrids=outGrids)


    ############################ Utils functions ###########################################

    def pdb2maestro(self, pdbIn, maeOut):
        prog = schrodinger_plugin.getHome('utilities/pdbconvert')
        args = '-ipdb {} -omae {}'.format(os.path.abspath(pdbIn), os.path.abspath(maeOut))
        try:
            self.runJob(prog, args, cwd=self._getExtraPath())
        except CalledProcessError as exception:
            # ask to Schrodinger why it returns a code 2 if it worked properly
            if exception.returncode != 2:
                raise exception
        return maeOut

    def getInputMaeFile(self):
        maeFile = glob.glob(self._getExtraPath('inputReceptor.mae*'))[0]
        return maeFile

    def getGridDir(self, pocketId):
        return self._getExtraPath(self.getGridName(pocketId))

    def getGridName(self, pocketId):
        return 'grid_{}'.format(pocketId)

    def getInnerBox(self, pocket):
        if self.innerAction.get() == 0:
            return self.innerX.get(), self.innerY.get(), self.innerZ.get()
        else:
            diam = pocket.getDiameter() * self.diameterNin.get()
            return diam, diam, diam

    def getOuterBox(self, pocket):
        if self.outerAction.get() == 0:
            return self.outerX.get(), self.outerY.get(), self.outerZ.get()
        else:
            diam = pocket.getDiameter() * self.diameterNout.get()
            return diam, diam, diam

    def getGridArgs(self, pocket):
        oBox, iBox = self.getOuterBox(pocket), self.getInnerBox(pocket)
        cMass = pocket.calculateMassCenter()
        return {'centerX': cMass[0], 'centerY': cMass[1], 'centerZ': cMass[2],
                'innerX': iBox[0], 'innerY': iBox[1], 'innerZ': iBox[2],
                'outerX': oBox[0], 'outerY': oBox[1], 'outerZ': oBox[2],
                'proteinFile': pocket.getProteinFile()}


    def _validate(self):
        errors = []
        if self.HbondAcceptHalo.get() and self.HbondDonorHalo.get():
            errors.append('Halogens cannot be simultaneously acceptors and donors of H-bonds')
        return errors

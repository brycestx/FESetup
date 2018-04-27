#  Copyright (C) 2012-2016  Hannes H Loeffler, Julien Michel
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#  For full details of the license please see the COPYING file
#  that should have come with this distribution.

r"""
A class to build a complex.  Derives from Common.

The complex Setup class composes a complex object from a protein and a ligand
object.  The class can create an AMBER topology file for the complex.
"""


__revision__ = "$Id: complex.py 580 2016-05-13 08:24:28Z halx $"



import FESetup
from FESetup import const, errors, logger
import utils

from ligand import Ligand
from protein import Protein
from common import *



class Complex(Common):
    """The complex setup class."""

    from FESetup.prepare.ligutil import flex as lig_flex

    SSBONDS_OFFSET = 1


    def __init__(self, protein, ligands):
        """
        :param protein: the protein for complex composition
        :type protein: Protein or string
        :param ligand: the ligand for complex composition
        :type ligand: Ligand or string
        :raises: SetupError
        """

        if type(ligands) is not list:
            ligands = [ligands]

        self.leap_added = False
        self.additional_bond_lines = []

        # FIXME: remove when ModelConfig is done
        #        this is still used for the Morph class
        if type(protein) == str and type(ligands[0]) == str:
            super(Complex, self).__init__(protein + const.PROT_LIG_SEP +
                                          ligands[0])

            self.protein_file = protein
            self.ligand_files = ligands
            self.ligands = [[] for i in range(len(ligands))]

            # FIXME: quick fix to allow dGprep to redo complex morph
            for i in range(len(ligands)):
                self.ligands[i] = Ligand(ligands[i], '')

            return

        assert type(protein) == Protein
        assert type(ligands) == list
        assert type(ligands[0]) == Ligand

        self.complex_name = protein.mol_name + const.PROT_LIG_SEP
        for i in range(len(ligands)):
            self.complex_name += ligands[i].mol_name
            if i < len(ligands)-1:
                self.complex_name += const.PROT_LIG_SEP

        super(Complex, self).__init__(self.complex_name)

        self.ligand_files = [[] for i in range(len(ligands))]
        self.protein_file = protein.orig_file
        self.charge = protein.charge

        if abs(self.charge) > const.TINY_CHARGE:
            logger.write('Warning: non-zero complex charge (%f)' % self.charge)

        self.protein = protein
        self.ligands = ligands
        self.frcmods = [[] for i in range(len(ligands))]
        for i in range(len(ligands)):
            self.frcmods[i] = os.getcwd() + "/" + const.LIGAND_WORKDIR + "/" + ligands[i].mol_name + "/" + ligands[i].frcmod
            self.ligand_files[i] = ligands[i].orig_file
            self.charge += ligands[i].charge



    def getResidueIndex(self, protein, residue):
        res_list = []
        name_set = set()
        for atom in protein.GetAtoms():
            name = oe.OEAtomGetResidue(atom).GetResidueNumber()
            if name not in name_set:
                name_set.add(name)
                res_list.append(name)
        res_list.sort()
        res_map = {}
        for i in range(len(res_list)):
            res_map[res_list[i]] = i+1
        return res_map[residue.GetResidueNumber()]


    @report
    def prepare_top(self, gaff='gaff', pert=None, add_frcmods=[], generate_zinc_bonds = False):
        """
        Prepare for parmtop creation i.e. add molecules to Leap structure.
        This needs to be run before create_top() to ensure that the molecule
        has been added but also to not trigger parmtop generation to early.
        Pmemd needs to have a second molecule added in alchemical free
        energy setups.
        :param generate_zinc_bonds: determine whether or not we should automatically create bond between Histidine-coordinated Zincs and the ligands
        """

            # ensure ligand is in MOL2/GAFF format
            lig = self.ligands[i]
            ac = os.getcwd() + "/../../" + const.LIGAND_WORKDIR + "/" + lig.mol_name + "/" + const.LIGAND_AC_FILE
            if os.access(ac, os.F_OK):

                mol_file = os.getcwd() + "/../../" + const.LIGAND_WORKDIR + "/" + lig.mol_name + "/" + const.GAFF_MOL2_FILE
                #frc = const.LIGAND_WORKDIR/lig.mol_name/lig.frcmod

                antechamber = utils.check_amber('antechamber')
                utils.run_amber(antechamber,
                                '-i %s -fi ac -o %s -fo mol2 -j 1 -at %s -pf y' %
                                (ac, mol_file, gaff) )
                self.ligand_fmt = 'mol2'
            else:
                # antechamber has trouble with dummy atoms
                mol_file = self.ligand_files[i]

                if self.ligand_fmt != 'mol2' and self.ligand_fmt != 'pdb':
                    raise errors.SetupError('unsupported leap input format: %s ('
                                            'only mol2 and pdb)' % self.ligand_fmt)

            if i == 0:
                frcmods = self.frcmods
                if add_frcmods:
                    frcmods.extend(add_frcmods)

            if not self.leap_added:
                print("  adding ligand '" + lig.mol_name + "' to complex")
                if i == 0:
                    self.leap.add_force_field(gaff)
                self.leap.add_mol(mol_file, self.ligand_fmt, frcmods, pert=pert)


    @report
    def create_top(self, boxtype='', boxlength=10.0, align=False,
                   neutralize=False, addcmd='', addcmd2='',
                   remove_first=False, conc = 0.0, dens = 1.0):
        """Generate an AMBER topology file via leap.

        :param boxtype: rectangular, octahedron or set (set dimensions explicitly)
        :param boxlength: side length of the box
        :param align: align solute along the principal axes
        :param neutralize: neutralise the system
        :param addcmd: inject additional leap commands
        :param remove_first: remove first unit/residue
        :param conc: ion concentration
        :type conc: float
        :param dens: expected target density
        :type dens: float
        :type boxtype: string
        :type boxlength: float
        :type align: bool
        :type neutralize: int
        :type addcmd: string
        :type remove_first: bool
        """

        if not self.leap_added:
            self.leap.add_mol(self.protein_file, 'pdb')
            self.leap_added = True

        if os.access('ffgen_leap.in', os.F_OK) and read_ffgen:
            self.amber_top = 'solvated.parm7'
            self.amber_crd = 'solvated.rst7'
            self.sander_crd = 'solvated.rst7'
            utils.run_leap(self.amber_top, self.amber_crd, program = 'tleap',
                           script = 'ffgen_leap.in')
        else:

            # FIXME: there can be problems with the ordering of commands, e.g.
            #        when tip4pew is used the frcmod files are only loaded after
            #        reading PDB and MOL2
            leapin = self._amber_top_common(boxtype, boxlength,
                                            neutralize, addcmd,
                                            addcmd2, align=align,
                                            remove_first=remove_first,
                                            conc=conc, dens=dens, additional_bond_lines=self.additional_bond_lines)

            utils.run_leap(self.amber_top, self.amber_crd, 'tleap', leapin)


    @report
    def prot_flex(self, cut_sidechain = 15.0, cut_backbone = 15.0):
        """
        Create a flexibility file for the protein describing how the input
        molecule can be moved by Sire.

        :param cut_sidechain: side chain cutoff
        :type cut_sidechain: float
        :param cut_backbone: backbone cutoff
        :type cut_backbone: float
        :raises: SetupError
        """


        if cut_sidechain < 0.0 or cut_backbone < 0.0:
            raise errors.SetupError('Cutoffs must be positive')


        import Sire.IO

        amber = Sire.IO.Amber()
        molecules, space = amber.readCrdTop(self.sander_crd, self.amber_top)

        moleculeNumbers = molecules.molNums()
        moleculeNumbers.sort()
        moleculeList = []

        for moleculeNumber in moleculeNumbers:
            molecule = molecules.molecule(moleculeNumber).molecule()
            moleculeList.append(molecule)

        id_first_nonligand = len(self.ligands)
        ligands = moleculeList[0:id_first_nonligand]
        not_ligand = moleculeList[id_first_nonligand:]

        sc_bb_residues = []

        logger.write('Computing flexible protein residues from %s' %
                     self.sander_crd)

        for cut in cut_sidechain, cut_backbone:
            cut_residues = []
            cut2 = cut**2

            for molecule in not_ligand:

                #for residue in molecule.residues():
                nmolresidues = molecule.nResidues()
                for z in range(0,nmolresidues):
                    residue = molecule.residues()[z]
                    # FIXME: a better way to skip unwanted residues would be to
                    # examine amber.zmatrices directly
                    # .value() returns a QtString!
                    if (str(residue.name().value()) not in
                        const.AMBER_PROTEIN_RESIDUES):
                        continue

                    shortest_dist2 = float('inf')

                    #for resat in residue.atoms():
                    nresatoms = residue.nAtoms()
                    for x in range(0,nresatoms):
                        resat = residue.atoms()[x]

                        if (resat.property('mass').value() <
                            const.MAX_HYDROGEN_MASS):
                            continue

                        rescoords = resat.property('coordinates')

                        for ligand in ligands:
                            #for ligat in ligand.atoms():
                            nligatoms = ligand.nAtoms()
                            for y in range(0,nligatoms):
                                ligat = ligand.atoms()[y]

                                if (ligat.property('mass').value() <
                                    const.MAX_HYDROGEN_MASS):
                                    continue

                                ligcoords = ligat.property('coordinates')
                                dist2 = space.calcDist2(rescoords, ligcoords)

                                if dist2 < shortest_dist2:
                                    shortest_dist2 = dist2

                    if shortest_dist2 < cut2:
                        cut_residues.append(residue)

            sc_bb_residues.append(cut_residues)

        lines = ['''# Flexible residues were only selected from the following list of residue names
# %s
# Cut-off for selection of flexible side chains %s Angstroms
# Number of residues with flexible side chains: %s
# Cut-off for selection of flexible backbone: %s Angstroms
# Number of residues with flexible backbone: %s
''' % (', '.join(const.AMBER_PROTEIN_RESIDUES), cut_sidechain,
       len(sc_bb_residues[0]), cut_backbone, len(sc_bb_residues[1])) ]

        htext = ['flexible sidechain', 'flexible backbone']

        for i in 0, 1:
            lines.append("%s\n" % htext[i])
            nums = []

            for residue in sc_bb_residues[i]:
                nums.append(residue.number().value() )

            nums.sort()

            line = ''
            for num in nums:
                if len(line) > 75:
                    lines.append('%s\n' % line)
                    line = ''
                line += ' %4d' % num

            lines.append('%s\n' % line)

        with open(const.PROTEIN_FLEX_FILE, 'w') as output:
            output.write(''.join(lines))



if __name__ == '__main__':
    pass

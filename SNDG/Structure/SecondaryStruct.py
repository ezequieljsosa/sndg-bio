import tempfile
from collections import defaultdict

from SNDG import execute, Struct
from Bio.SeqUtils import  seq3

class SecondaryStruct:
    '''
    http://swift.cmbi.ru.nl/gv/dssp/HTML/descrip.html
    '''

    def __init__(self, pdb_path):
        self.pdb_path = pdb_path
        self.dssp = []

    def run_dssp(self):
        out = tempfile.mkstemp(suffix=".dssp")[1]
        execute("dssp -i {pdb_path} -o {out}",
                pdb_path=self.pdb_path, out=out)

        with open(out) as h:
            start = False
            for l in h:
                if start:
                    res = int(l[5:10])
                    aa = l[10:14].strip()
                    ss = l[14:17].strip()
                    bbl1 = l[23:24]
                    bbl2 = l[24:25]
                    bp1 = int(l[25:29])
                    bp2 = int(l[29:33])
                    bslabel = l[33:34]
                    self.dssp.append(Struct(res=res, aa=aa, ss=ss, bp1=bp1, bp2=bp2,
                                            bbl1=bbl1, bbl2=bbl2, bslabel=bslabel))
                else:
                    if l.startswith("  #  RESIDUE AA"):
                        start = True

    def new_with_secondary(self, output_pdb_path):
        """
        Creates a new pdb file with secondary structure
        :param output_pdb_path: output pdb file
        """
        header_helix = self.header_helix()
        header_beta = self.header_beta()
        with open(output_pdb_path,"w") as w:
            with open(self.pdb_path) as r:
                header_added = False
                for l in r:
                    if l.startswith("ATOM") and not header_added:
                        header_added = True
                        for h in header_helix:
                            w.write(h + "\n")
                        for h in header_beta:
                            w.write(h + "\n")
                    w.write(l)

    def header_beta(self):
        """
        http://www.wwpdb.org/documentation/file-format-content/format33/sect5.html
        http://www.csb.yale.edu/userguides/databases/dssp/dssp_man.html
        """

        SHEET = "SHEET "                # 1 -  6        Record name   "SHEET "
        s1 = " "
        strand =  lambda x: x           # 8 - 10        Integer       strand         Strand  number which starts at 1 for each
                                        #     strand within a sheet and increases by one.
        s2 = " "
        sheetID = lambda x: x           # 12 - 14        LString(3)    sheetID        Sheet  identifier.
        numStrands =  lambda x: x       # 15 - 16        Integer       numStrands     Number  of strands in sheet.
        s3 = " "
        initResName = lambda x: x       # 18 - 20        Residue name  initResName    Residue  name of initial residue.
        s4 = " "
        initChainID = lambda x: x       # 22             Character     initChainID    Chain identifier of initial residue
                                        # in strand.
        initSeqNum =  lambda x: x       # 23 - 26        Integer       initSeqNum     Sequence number of initial residue
                                        # in strand.
        initICode =   lambda x: x       # 27             AChar         initICode      Insertion code of initial residue
                                        # in  strand.
        s5 = " "
        endResName =  lambda x: x       # 29 - 31        Residue name  endResName     Residue name of terminal residue.
        s6 = " "
        endChainID =  lambda x: x       # 33             Character     endChainID     Chain identifier of terminal residue.
        endSeqNum =   lambda x: x       # 34 - 37        Integer       endSeqNum      Sequence number of terminal residue.
        endICode =    lambda x: x       # 38             AChar         endICode       Insertion code of terminal residue.
        sense =       lambda x: x       # 39 - 40        Integer       sense          Sense of strand with respect to previous
                                        # strand in the sheet. 0 if first strand,
                                        # 1 if  parallel,and -1 if anti-parallel.
        s7 = " "
        curAtom =     lambda x: x       # 42 - 45        Atom          curAtom        Registration.  Atom name in current strand.
        curResName  = lambda x: x       # 46 - 48        Residue name  curResName     Registration.  Residue name in current strand
        s8 = " "
        curChainId =  lambda x: x       # 50             Character     curChainId     Registration. Chain identifier in
                                        # current strand.
        curResSeq =   lambda x: x       # 51 - 54        Integer       curResSeq      Registration.  Residue sequence number
                                        # in current strand.
        curICode =    lambda x: x       # 55             AChar         curICode       Registration. Insertion code in
                                        # current strand.
        s9 = " "
        prevAtom =     lambda x: x      # 57 - 60        Atom          prevAtom       Registration.  Atom name in previous strand.
        prevResName =  lambda x: x      # 61 - 63        Residue name  prevResName    Registration.  Residue name in
                                        # previous strand.
        s10 = " "
        prevChainId =  lambda x: x      # 65             Character     prevChainId    Registration.  Chain identifier in
                                        # previous  strand.
        prevResSeq =   lambda x: x      # 66 - 69        Integer       prevResSeq     Registration. Residue sequence number
                                        # in previous strand.
        prevICode =    lambda x: x      # 70             AChar         prevICode      Registration.  Insertion code in
                                        # previous strand.


        current = "-"

        #TODO init num strands
        betas = defaultdict(lambda:{})
        for ss_entry in self.dssp:
            if "first" not in betas[ss_entry.bslabel]:
                betas[ss_entry.bslabel]["first"] = ss_entry
            if "strands" not  in betas[ss_entry.bslabel]:
                betas[ss_entry.bslabel]["strands"] = []
            betas[ss_entry.bslabel]["strands"].append(  (ss_entry.bbl1 , ss_entry.bbl2)  )
        header = []

        size = None
        for ss_entry in self.dssp:
            ss = ss_entry.ss
            sheet = ss_entry.bslabel
            if current != ss_entry.ss:
                if current == "E":
                    newline += endResName(ss_entry.aa) + s6 + endChainID(" ")
                    newline += endSeqNum(ss_entry.res) + endICode(" ")

                    if betas[sheet] == 1:
                        newline += " 0"
                    else:
                        newline += sense(ss_entry) + s7 + curAtom("N") + curResName(ss_entry.aa) + s8 + curChainId(" ")
                        newline += curResSeq() + curICode() + s9 + prevAtom("O") + prevResName() + s10
                        newline += prevChainId() + prevResSeq() + prevICode()
                    header.append(newline)
                current = ss_entry.ss
                if ss == "E":
                    strand = 0
                    beta += 0
                    newline = SHEET + s1 + strand(betas[sheet]) + s2 + sheetID(sheet)
                    newline += numStrands(sheet) + s3 + initResName(ss_entry.aa) + s4
                    newline += initChainID(" ") + initSeqNum(ss_entry.res) + initICode(" ") + s5

            elif current == ss:
                if ss in ["H", "G", "I"]:
                    size += 1

    def header_helix(self):
        """
        http://www.wwpdb.org/documentation/file-format-content/format33/sect5.html
                                             CLASS NUMBER
        TYPE OF  HELIX                     (COLUMNS 39 - 40)
        --------------------------------------------------------------
        Right-handed alpha (default)                1  --> DSSP: H	Alpha helix (4-12)
        Right-handed pi                             3  --> DSSP: I	Pi helix
        Right-handed 3 - 10                         5  --> DSSP: G	3-10 helix
        """

        HELIX = "HELIX "  # 1 -  6        Record name    "HELIX "
        s1 = " "
        serNum = lambda helix: str(helix).rjust(3)  # 8 - 10  serNum Serial number of the helix. This starts
        # at 1  and increases incrementally.
        s2 = " "
        helixID = serNum  # 12 - 14        LString(3)     helixID       Helix  identifier. In addition to a serial
        #   number, each helix is given an
        #   alphanumeric character helix identifier.
        s3 = " "
        initResName = lambda aa: seq3(aa).upper().rjust(3)  # 16 - 18  Residue name   initResName   Name of the initial residue.
        s4 = " "
        initChainID = lambda chain: chain  # 20  initChainID   Chain identifier for the chain containing this  helix.

        s5 = " "
        initSeqNum = lambda res: str(res).rjust(4)  # 22 - 25  initSeqNum    Sequence number of the initial residue.
        initICode = lambda icode: icode  # 26  AChar     initICode     Insertion code of the initial residue.
        s6 = " "
        endResName = initResName  # 28 - 30 Residue  name  endResName    Name of the terminal residue of the helix.
        s7 = " "
        endChainID = initChainID  # 32  endChainID    Chain identifier for the chain containing this  helix.
        s8 = " "
        endSeqNum = initSeqNum  # 34 - 37         endSeqNum     Sequence number of the terminal residue.
        endICode = initICode  # 38             AChar          endICode      Insertion code of the terminal residue.
        helixClass = lambda ss: " 1" if ss == "H" else (
            " 3" if ss == "I" else " 5")  # 39 - 40        helixClass    Helix class (see below).
        comment = " ".rjust(30)  # 41 - 70        comment       Comment about this helix.
        s9 = " "
        length = lambda l: str(l).rjust(5)  # 72 - 76    length        Length of this helix.

        current = "-"
        helix = 0
        header = []
        size = None
        for ss_entry in self.dssp:
            ss = ss_entry.ss
            if current != ss_entry.ss:
                if current in ["H", "G", "I"]:
                    newline += endResName(ss_entry.aa) + s7 + endChainID(" ") + s8
                    newline += endSeqNum(ss_entry.res) + endICode(" ")
                    newline += helixClass(ss) + comment + s9 + length(size)
                    header.append(newline)
                current = ss_entry.ss
                if ss in ["H", "G", "I"]:
                    size = 1
                    helix += 1
                    newline = HELIX + s1 + serNum(helix) + s2 + helixID(helix) + s3
                    newline += initResName(ss_entry.aa) + s4 + initChainID(" ") + s5
                    newline += initSeqNum(ss_entry.res) + initICode(" ") + s6

            elif current == ss:
                if ss in ["H", "G", "I"]:
                    size += 1
        return header




if __name__ == '__main__':
    from SNDG import init_log

    init_log()
    secs = SecondaryStruct("test/P9WKI1_1gr0A.pdb")
    secs.run_dssp()
    secs.new_with_secondary("test/P9WKI1_1gr0A_V2.pdb")

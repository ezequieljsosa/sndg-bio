'''
Created on Apr 21, 2014

@author: eze
'''

import os
import sys
import logging

from modeller import alignment, model, environ, profile
from modeller.automodel import assess, refine
from modeller.automodel.automodel import automodel
from modeller.parallel import local_slave, job

from Bio.PDB.PDBParser import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from SNDG.Structure.ChainSplitter import ChainSplitter

_log = logging.getLogger(__name__)


class ModellerOverflowError(Exception):
    '''
    https://salilab.org/archives/modeller_usage/2013/msg00084.html
    
On 06/03/2013 10:56 PM, Arjun Ray wrote:
I am trying to mutate to ARG26 and getting the following error:

randomi_498_> Atoms,selected atoms,random_seed,amplitude:     2031
11        1        4.0000
randomi_496_> Amplitude is > 0; randomization is done.
check_inf__E> Atom 206 has out-of-range coordinates (usually infinity).
               The objective function can thus not be calculated.

This typically happens when huge forces on your system cause it to become unstable. This can sometimes happen if you have a very bad starting geometry. I can't tell anything further without seeing your input files.

    Ben Webb, Modeller Caretaker
    
    -------------------------
    
    https://salilab.org/archives/modeller_usage/2013/msg00160.html
    
On 10/10/13 8:21 AM, Anna Marabotti wrote:
I created several mutants of my protein (which is a homodimeric protein)
using the script mutate_model.py.and all has gone perfectly, but I have
a problem in introducing a mutation.
...
OverflowError:check_inf__E> Atom 5347 has out-of-range coordinates
(usually infinity).
               The objective function can thus not be calculated.

Usually you can avoid this problem by editing the script and changing the value of rand_seed.

    Ben Webb, Modeller Caretaker
    
    
    
    '''

    def __init__(self, outputs, base_models):
        self.outputs = outputs
        self.base_models = base_models


class Modeller(object):
    '''
    classdocs
    '''

    @staticmethod
    def process_modpipe_profiles():
        # https://salilab.org/modbase-download/projects/genomes/M_tuberculosis/2018/m_tuberculosis_2018.tar
        # https://salilab.org/modbase-download/projects/genomes/L_major/2016/l_major_2016.tar
        # https://salilab.org/modbase-download/modbase_models_academic-20120928.txt.gz
        # https://salilab.org/modeller/supplemental.html
        # https://salilab.org/modeller/downloads/20200513_pdb95.pir.gz
        # https://salilab.org/modeller/downloads/20200513_pdb95_profiles.tar.bz2 -> data/i8/1i8oA

        env = environ()
        # /mnt/data/data/databases/pdb/modpipe/data/i8/1i8oA
        prf = profile(env, file="./1i8oA-uniprot90.prf", profile_format="TEXT")
        aln = prf.to_alignment()
        aln.write(file='pepe.fasta', alignment_format='FASTA')

    """
    data = "./pdball.pir".read() 
    from Bio.SeqRecord import SeqRecord 
    from Bio.Seq import Seq 
    with open("modeller.fasta","w") as h: 
        for record in data.split("C; Produced by MODELLER"): 
             record = [line.strip() for line in record.split("\n") if line.strip()] 
             if record:             
                 seqrecord = SeqRecord(id=record[0][1:].strip(), 
                     name="",description=record[1].strip(), 
                     seq=Seq( "".join(record[2:])  )) 
                 bpio.write(seqrecord,h,"fasta") 
    """

    def __init__(self, base_path, pdbs_dir):
        '''
        Constructor
        '''
        self.model_count = 4
        self._assess_methods = [assess.DOPE, assess.GA341]
        self._refinement = refine.slow
        self.parallel_jobs = 4
        self.pdbs_dir = pdbs_dir
        self.base_path = base_path

    def log_all(self):
        from modeller.util.logger import log
        log.level(1, 1, 1, 1, 0)

    def _create_aligment(self, env, base_models):

        _log.debug("creating alignments for %s with %s pdbs" % (self.seqrecord.id, len(base_models)))
        aligned_models = []

        env.io.atom_files_directory = [self.out_folder + '/']

        aln = alignment(env)

        aln.append_sequence(str(self.seqrecord.seq))
        aln[0].code = str(self.seqrecord.id)

        for i, pdb_chain_file_path in enumerate(base_models, 1):
            # TODO sacar parseo feo
            code = pdb_chain_file_path.split("/")[-1].replace(".ent", "").replace("pdb", "")
            m = model(env, file=code)
            aln.append_model(m, align_codes=code)
            aln[i].code = code
            aligned_models.append(code)

        aln.malign()
        aln.id_table(matrix_file=self.seqrecord.id + '_family.mat')

        aln.write(file=self.model_directory() + "/" + self.seqrecord.id + '.ali', alignment_format='PIR')
        assert os.path.exists(
            self.model_directory() + "/" + self.seqrecord.id + '.ali'), "NOOOOOOOOOOOO!!!!:  " + os.getcwd() + "/" + self.seqrecord.id + '.ali'
        aln.write(file=self.model_directory() + self.seqrecord.id + '.pap', alignment_format='PAP')
        return aligned_models

    def aln2pir(self, model_id, alignment, query, template):

        p = PDBParser()

        # En el alienamiento que se toma como entrada, la primer
        # secuencia tiene que ser el target y la segunda el template
        pdb_chain = template.id.split("_")[1]
        pdb = str("_".join(template.id.split("_")[0:2]))
        residues = list(p.get_structure(pdb, os.path.join(self.pdbs_dir, pdb + ".pdb")).get_residues())
        aligned_residues = residues[alignment.aln_hit.start:]
        pdb_start = str(aligned_residues[0].id[1])
        #         pdb_start = str( int(template.id.split("_")[2]) + alignment.aln_hit.start)
        #         pdb_end = str(int(template.id.split("_")[3]) -  alignment.aln_hit.start)
        pdb_end = "+" + str(len(str(template.seq).replace("-", "").replace(".", "")))
        #         pdb_end = "+" + str( len(aligned_residues)  )
        start = int(alignment.aln_query.start) + 1
        length = str(int(alignment.aln_query.end) - int(alignment.aln_query.start))
        targetheader = f">P1;{query.id}\nsequence:{query.id}:{start}:A:+{length}:A::::\n"
        templateheader = ">P1;%s\nstructureX:%s:%s:%s:%s:%s::::"

        # targetrec = alignment["query"]
        # templaterec = alignment["hit"]

        outfile = self._aln_file(model_id, query.id)

        with open(outfile, "w") as piralignment:
            piralignment.write(targetheader)
            piralignment.write(str(query.seq) + "*" + "\n")
            piralignment.write(templateheader % (pdb, pdb, pdb_start, pdb_chain, pdb_end, pdb_chain) + "\n")
            piralignment.write(str(template.seq) + "*" + "\n")


        return [pdb]

    def _create_models(self, env, base_models, model_id, query_id):

        class MyModel(automodel):
            def special_patches(self, aln):
                # Rename both chains and renumber the residues in each
                start = int(aln[1].range[0].split(":")[0])
                self.rename_segments(segment_ids=['A'], renumber_residues=start)
                # self.chains[0].name = 'A'

        a = MyModel(env, alnfile=self._aln_file(model_id, query_id),
                    knowns=base_models, sequence=str(query_id),
                    assess_methods=self._assess_methods)
        a.starting_model = 1
        a.ending_model = self.model_count
        a.md_level = self._refinement

        # TODO assessing with {assess_methods}
        _log.debug("Generating {model_count} model/s for {orf} with {base_models} ".format(
            orf=query_id, base_models=base_models, model_count=self.model_count))

        if self.parallel_jobs > 1:
            j = job()
            for _ in range(0, self.parallel_jobs):
                j.append(local_slave())
            a.use_parallel_job(j)

        model = a.make()

        if "OverflowError" in str(a.outputs[0]["failure"].__class__):
            raise ModellerOverflowError(a.outputs[0]["failure"], base_models)
        _log.debug("{model_count} models for {orf} with {base_models} where generated ".format(
            orf=query_id, base_models=base_models, model_count=self.model_count))

    def create_model(self, model_id, alignment):
        # base_models_tuple, target_id,pdb_id,pdb_chain,aln_target,aln_pdb
        current_dir = os.getcwd()
        org_stdout = sys.stdout  # store original stdout object for later

        # query,template = (alignment.query_record(),alignment.hit_record())
        query, template = (SeqRecord(Seq(alignment.aln_query.seq), id=alignment.aln_query.name),
                           SeqRecord(Seq(alignment.aln_hit.seq), id=alignment.aln_hit.name))

        self._create_directory_structure(model_id, query.id)

        sys.stdout = open(self.log_path(model_id, query.id), "w")
        try:
            env = environ()
            env.io.atom_files_directory = self.pdbs_dir

            #             #align_models = self._create_aligment( env, base_models)
            align_models = self.aln2pir(model_id, alignment, query, template)

            os.chdir(self.model_directory(model_id, query.id))
            self._create_models(env, align_models, model_id, query.id)
            #

            return [self.pdb_path(model_id, query.id, i) for i in range(1, self.model_count + 1)]
        #
        #
        #
        finally:
            sys.stdout.close()
            sys.stdout = org_stdout
            os.chdir(current_dir)

    def file_list(self, model_id, query_id):
        for i in range(1, self.model_count + 1):
            yield self.pdb_path(model_id, query_id, i)

    def _aln_file(self, model_id, query_id):
        return str(os.path.join(self.model_directory(model_id, query_id), 'seq.ali'))

    def pdb_path(self, model_id, query_id, i=1):
        return os.path.join(self.model_directory(model_id, query_id),
                            query_id + ".B9999" + str(i).zfill(4) + ".pdb")

    def log_path(self, model_id, query_id):
        return os.path.join(self.model_directory(model_id, query_id), "process._log")

    def model_directory(self, model_id, query_id):
        return os.path.join(self.base_path, query_id, model_id)

    def _create_directory_structure(self, model_id, query_id):
        if not os.path.exists(self.model_directory(model_id, query_id)):
            os.makedirs(self.model_directory(model_id, query_id))


if __name__ == '__main__':
    from SNDG import init_log, Struct

    init_log()

    workdir = "/media/eze/Data/data/organismos/Pext14-3B/analysis/struct/good"
    modeler = Modeller(workdir, "/tmp")
    model_id = "PE143B_RS25640_3u52_B_6_498"
    alignment = Struct(
        aln_query=Struct(name="PE143B_RS25640",
                         seq="KKLNAKDKYRLLTRDLAWEPSYRTEEEIFPYIAYEGLKIHDWNKWEDPFRLTMDAYWKYQAEKERKFYAIIDAHAQNNGHLNITDARYLSALKIFLQAISPGEYAAHKGFARAGREFRGVGTQVACQMQAIDELRHAQTQIHALSNYNKFYNGFHAFADQRDRIWYTSVARSFFDDAMSAGPFEFMIAIGFSFEYVLTNLLFVPFMSGAAYNGDMATVTFGFSAQSDEARHMTLGLECIKFMLEQDPANLPIVQGWIDKWFWRGFRVLGLVSTMMDYMLPKRVMSWREAWEIYGAENGGALFKDLARYGIRPPKSWDDAEASIDHMSHQFMLGLYQWSFGTAFHAWIPSDDDMQWLSAKYPTTFDKYYRPRWEHIKKMEAAGTPFKNYGLAKLCQCCQLPTVFTEPDDPTLICHRQVQYKGDKYHFCSDHCMGIFNNEPEKYIQAWLPMPALFQAPTN-GDLGAWMD-WVSLKDGQDNGDFADSQDRRN",
                         start=7, end=494),  # .ungap("-")
        aln_hit=Struct(name="3u52_B_6_498",
                       seq="KKLNLKDKYQYLTRDMAWEPTYQDKKDIFPEEDFEGIKITDWSQWEDPFRLTMDAYWKYQAEKEKKLYAIFDAFAQNNGHQNISDARYVNALKLFISGISPLEHAAFQGYSKVGRQFSGAGARVACQMQAIDELRHSQTQQHAMSHYNKHFNGLHDGPHMHDRVWYLSVPKSFFDDARSAGPFEFLTAISFSFEYVLTNLLFVPFMSGAAYNGDMATVTFGFSAQSDEARHMTLGLEVIKFILEQHEDNVPIVQRWIDKWFWRGFRLLSLVSMMMDYMLPNKVMSWSEAWEVYYEQNGGALFKDLERYGIRPPKYQDVANDAKHHLSHQLWTTFYQYCQATNFHTWIPEKEEMDWMSEKYPDTFDKYYRPRYEYLAKEAAAGRRFYNNTLPQLCQVCQIPTIFTEKDAPTMLSHRQIEHEGERYHFCSDGCCDIFKHEPEKYIQAWLPVHQIYQGNCEGGDLETVVQKYYHINIGEDNFDYVGSPDQKH",
                       start=0, end=489))
    ChainSplitter("/tmp/").make_pdb(pdb_path="/data/databases/pdb/divided/u5/pdb3u52.ent", pdb_id="3u52", chain="B",
                                    overwrite=True)
    models = modeler.create_model(model_id, alignment)

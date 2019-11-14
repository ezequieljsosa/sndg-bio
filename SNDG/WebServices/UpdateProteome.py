'''
Created on Jun 15, 2017

@author: eze
'''

import subprocess as sp
import os
import sys
import logging


from Bio import Entrez
import Bio.SeqIO as bpio 
from SNDG import init_log
from SNDG.WebServices.NCBI import NCBI

from SNDGInt.Submitter import ExternalAssembly
from BIADeploy.BioMongoDB import BioMongoDB
from BIADeploy.BiaSql import BiaSql
from Bia.External.uniprot import Uniprot
from BIADeploy.UniprotAnnotator import UniprotAnnotator
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bia.Model.sequence import Protein
from VARDB.Aln2Sql import Aln2Sql
from SNDGComp.Taxon import Tax, tax_db
from peewee import MySQLDatabase, DoesNotExist

from Bia.Programs.Search.hmmer import Hmmer
from BiaOntology import so
from Bia.Model.collections import SeqCollection
from BIADeploy.JBrowse import JBrowse
import traceback
from mysql.connector.errors import OperationalError
import shutil
from BCBio import GFF



_log = logging.getLogger(__name__)

class Uniprot2(object):
    '''

    '''
    

    def __init__(self):
        '''
        Constructor
        '''
        self.mdb = "saureus"
        self.sqldb = "bioseq"
        
        Entrez.email = 'A.N.Other@example.com' 
        self.ncbi = NCBI()
        self.uniprot = Uniprot()
        self.uniprot_annotator = UniprotAnnotator()
    
    def init(self):
        self.mdb = BioMongoDB(self.mdb)

    
    def save_sqldb(self,seq_col_name,annotation_tax=None,tmp_dir="/tmp/" ):

        if not self.mdb.seq_col_exists(seq_col_name):
               
            if  (not self.sqldb.seq_col_exists(seq_col_name + "_unip_tr")             
                and not os.path.exists(tmp_dir + "/tr.xml")):    
                _log.debug("Downloading uniprot annotation for tax %s" % str(annotation_tax))
                sp_path,tr_path =  self.uniprot.download_proteome_from_tax(annotation_tax, tmp_dir)
            else:
                sp_path,tr_path = (tmp_dir + "/sp.xml",tmp_dir + "/tr.xml")
            
            if not self.sqldb.seq_col_exists(seq_col_name + "_unip_sp"):    
                if os.path.getsize( sp_path ) > 0 :    
                    _log.debug("Loading swissprot to sqldb" )                    
                    self.sqldb.save_collection(seq_col_name+ "_unip_sp", bpio.parse(sp_path,"uniprot-xml"))
                    os.remove(sp_path)
            if not self.sqldb.seq_col_exists(seq_col_name + "_unip_tr"):
                if os.path.getsize( tr_path ) > 0 :
                    _log.debug("Loading trembl to sqldb" )
                    self.sqldb.save_collection(seq_col_name+ "_unip_tr", bpio.parse(tr_path,"uniprot-xml"))
                    os.remove(tr_path)
        
        
            _log.debug("Creanting mongodb collection..." )
            self.mdb.load_from_db(seq_col_name, seq_col_sql)
        
        
        
        protein_fasta = tmp_dir + "/proteins.fasta"
        if not os.path.exists(protein_fasta) or (not os.path.getsize(protein_fasta)):
            with open(protein_fasta,"w" ) as h:
                    for p in Protein.objects(organism=seq_col_name).no_cache():
                        bpio.write(SeqRecord(id=p.gene[0],seq=Seq(p.seq)), h, "fasta")

        
        genome = SeqCollection.objects(name=seq_col_name).get()
        genome.ncbi_assembly = seq_col_name
        if not genome.statistics:

            
            self.mdb.index_seq_collection(seq_col_name, pathways=False)
            self.mdb.build_statistics(seq_col_name)
            
            _log.info("Sequence collection %s created correctly " % seq_col_name)

if __name__ == '__main__':
    init_log("/tmp/sndg2.log")
    logger = logging.getLogger('peewee')
    logger.setLevel(logging.ERROR)
    dep = Deployer()
    connect_to_db(password="mito")
    import mysql.connector
    dep.mdb = "saureus"
    #dep.annotation_tax = "158879"
    dep.init()
    tax_db.initialize(MySQLDatabase('bioseq', user='root',passwd="mito"))
    
    





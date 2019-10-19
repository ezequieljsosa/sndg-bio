from tqdm import tqdm
import traceback
from SNDG.Comparative.Pangenome import Pangenome, Strain, sqldb
from SNDG.WebServices.NCBI import NCBI
from Bio import Entrez, SeqIO
import multiprocessing

from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(driver="MySQLdb", user="root",
                                      passwd="mito", host="localhost", db="bioseqdb")


def save_sequences(strain):
    db = server[strain.acc]
    sequencesList = Entrez.read(Entrez.esearch(
        db="nuccore", term=strain.acc + "[Assembly]", idtype="acc", retmax=100000))["IdList"]

    n = 150
    list_of_lists = [sequencesList[i:i+n] for i in range(0, len(sequencesList), n)]
    def seq_iterator():
        for seq_ids in list_of_lists:
            handle = Entrez.efetch(db="nuccore", id=",".join(seq_ids), rettype="gbwithparts", retmode="text")
            for seq in SeqIO.parse(handle, "genbank"):
                yield seq


    pbar2.set_description(genome.name)

    with tqdm(seq_iterator(), total=len(sequencesList)) as pbar3:
        pbar3.set_description("loading sequences")
        db.load(pbar3)
    server.commit()


if __name__ == '__main__':


    from peewee import MySQLDatabase

    mysql_db = MySQLDatabase('bioseqdb', user="root", password="mito")
    sqldb.initialize(mysql_db)
    tables = [Pangenome,Strain]

    # for x in tables:
    #     x.create_table()

    Entrez.email = "Your.Name.Here@example.org"
    genomes = list(Pangenome.select())
    with tqdm(genomes) as pbar1:
        for genome in pbar1:
            pbar1.set_description(genome.name)
            strains = list(Strain.select().where(Strain.pangenome == genome))

            with tqdm(strains) as pbar2:
                for strain in pbar2:
                    try:
                        if strain.acc not in server:

                            db = server.new_database(strain.acc, description="")
                            server.commit()

                            def save_strain():
                                save_sequences(strain)
                            p = multiprocessing.Process(target=save_strain)
                            p.start()
                            p.join(80)
                            if p.is_alive():
                                p.terminate()
                                p.join()
                                raise multiprocessing.TimeoutError()

                            strain.loaded = True
                            strain.save()

                    except Exception as ex:

                                traceback.print_exc()
                                server.rollback()
                                if strain.acc in server:
                                    server.remove_database(strain.acc)
                                    server.commit()
                                mysql_db = MySQLDatabase('bioseqdb', user="root", password="mito")
                                server = BioSeqDatabase.open_database(driver="MySQLdb", user="root",
                                                                      passwd="mito", host="localhost", db="bioseqdb")


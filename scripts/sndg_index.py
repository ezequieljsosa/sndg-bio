#!/usr/bin/python
# encoding: utf-8
'''
scripts.load_bia_proteome -- shortdesc

scripts.load_bia_proteome is a description

It defines classes_and_methods

@author:     Ezequiel Sosa
@copyright:  2015 BIA. All rights reserved.
@license:    license
@contact:    user_email

'''
import sys
import logging
from tqdm import tqdm
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from peewee import fn
import pymongo
from SNDG import init_log
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.Taxon import Tax, TaxName
from SNDG.BioMongo.Model.Taxonomy import Taxonomy

from mongoengine.errors import DoesNotExist

init_log()
_log = logging.getLogger(__name__)

__all__ = []
__version__ = 0.1
__date__ = '2015-09-18'
__updated__ = '2015-09-18'


def get_or_load_by_name(tax_name, tax_cache):
    return get_or_load(tax_name, tax_cache,
                       lambda x: Tax.select().join(TaxName).where(fn.Lower(TaxName.name) == fn.Lower(x)).get()
                       )


def get_or_load_by_id(tax_id, tax_cache):
    def by_id(x):
        return Tax.select().where(Tax.ncbi_taxon_id == x).get()
    return get_or_load(tax_id, tax_cache,by_id
                       )


def get_or_load(tax_name, tax_cache, tax_fn):
    if tax_name in tax_cache:
        return tax_cache[tax_name]
    else:
        try:
            tax = tax_fn(tax_name)

            t = Taxonomy(name=[x.name for x in tax.names if x.name_class == "scientific name"][0],
                         names=[x.name for x in tax.names], ncbi_taxon_id=tax.ncbi_taxon_id, rank=tax.node_rank,
                         genetic_code=tax.genetic_code, mito_genetic_code=tax.mito_genetic_code)
            t.keywords = reduce(lambda x, y: x + y, [x.name.lower().split() for x in tax.names])
            t.keywords.append(str(t.ncbi_taxon_id))
            for parent in Tax.parents(tax, return_self=False):
                names = {x.name_class: x.name for x in parent.names}
                t.parents = [names["scientific name"]] + t.parents
                t.keywords += reduce(lambda x, y: x + y, [x.lower().split() for x in names.values()])
                t.keywords.append(str(parent.ncbi_taxon_id))
            t.keywords = list(set(t.keywords))
            tax_cache[tax_name] = t
            t.save()
            return t
        except Tax.DoesNotExist:
            _log.warn(str(tax_name) + " does not exists")
            return None


def update_element(val, collection, elem, idx_name, tax_cache,org):
    if val:
        collection.update({"_id": elem["_id"]}, {"$set": {idx_name + ".tax": list(val.keywords)}})
        collection.update({"_id": elem["_id"]}, {"$addToSet": {"keywords": {"$each":
                                                                                list(val.keywords)}}})
    else:
        tax_cache[str(org).lower()] = None
        _log.warn(str(org )+ " not found")


def main(argv=None):  # IGNORE:C0111

    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by user_name on %s.
  Copyright 2015 BIA. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                        help="set verbosity level [default: %(default)s]")

    parser.add_argument("-db_structure", "--db_structure", help="Mongo structure db", default='pdb')
    parser.add_argument("-db_genome", "--db_genome", help="Mongo proteins db", default='saureus')
    parser.add_argument('-o', '--overwrite', default=True, action='store_true')
    parser.add_argument("-host", "--db_host", default='127.0.0.1')
    parser.add_argument('-V', '--version', action='version', version=program_version_message)

    args = parser.parse_args()

    db = pymongo.MongoClient(args.db_host)[args.db_structure]
    BioMongoDB(args.db_genome)
    logging.getLogger("peewee").setLevel(logging.WARN)
    from peewee import MySQLDatabase
    from SNDG.BioMongo.Process.Taxon import tax_db
    tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))

    tax_cache = {}
    for t in Taxonomy.objects().no_cache():
        for n in t.names:
            tax_cache[n.lower()] = t
            tax_cache[t.ncbi_taxon_id] = t
    query = {}
    idx_name = "sndg_index"
    if not args.overwrite:
        query = {idx_name: {"$exists": 0}}
    # total = db.structures.count(query)
    # with tqdm(db.structures.find(query, {"organism": 1}), total=total) as pbar:
    #     for struct in pbar:
    #         if "organism" in struct:
    #             for org in [x for x in set(struct["organism"].lower().split(";") + struct["organism"].lower().split(",") +
    #                   [struct["organism"].lower().split("(")[0]]) if ";" not in x and "," not in x and "(" not in x]:
    #                 org = org.strip()
    #                 val = get_or_load_by_name(org, tax_cache)
    #                 if val:
    #                     db.structures.update({"_id": struct["_id"]}, {"$set": {idx_name + ".tax": list(val.keywords)}})
    #                 else:
    #                     tax_cache[org.lower()] = None
    #                     _log.warn(org + " not found")


    # db.structures.update({"ligands.0":{"$exists",1}},  {"$set": {idx_name + ".ligand": 1}},multi=True);

    db = pymongo.MongoClient(args.db_host)[args.db_genome]

    # total = db.barcodes.count(query)
    # with tqdm(db.barcodes.find(query, {"tax": 1}), total=total) as pbar:
    #     for barcode in pbar:
    #         val = get_or_load_by_id(barcode["tax"], tax_cache)
    #         update_element(val, db.barcodes, barcode, idx_name, tax_cache,barcode["tax"])

    total = db.sequence_collection.count(query)
    with tqdm(db.sequence_collection.find(query, {"name": 1, "tax": 1,"assemblyStatus":1},no_cursor_timeout=True), total=total) as pbar:
        for genome in pbar:
            if "tax" in genome:
                val = get_or_load_by_id(int(genome["tax"]["tid"]), tax_cache)
                update_element(val, db.sequence_collection, genome, idx_name, tax_cache,genome["tax"]["tid"])
                if val:
                    select = {"organism": genome["name"]}
                    kws = list(val.keywords)
                    db.proteins.update(select, {"$set": {idx_name + ".tax": kws}},multi=True)
                    db.proteins.update(select, {"$addToSet": {"keywords": {"$each":
                                                                               kws}}},multi=True)
                    db.contig_collection.update(select, {"$set": {idx_name + ".tax": kws,
                                                                  idx_name + ".assemblyStatus" :genome["assemblyStatus"]

                                                                  }},multi=True)
                    db.contig_collection.update(select, {"$addToSet": {"keywords": {"$each":
                                                                                    kws}}},multi=True)


    print ("Ok")


if __name__ == "__main__":
    sys.exit(main())

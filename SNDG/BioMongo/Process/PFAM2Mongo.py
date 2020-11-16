"""

"""
from SNDG.BioMongo.Model.Ontology import Ontology
from SNDG.BioMongo.Process.KeywordIndexer import KeywordIndexer


class PFAM2Mongo(object):
    """

    """

    def __init__(self, pfam_file_path, ontology_name='pfam'):

        self.pfam_file_path = pfam_file_path
        self.ontology_name = ontology_name
        self.ki = KeywordIndexer()

    def load(self):
        ont_doc = None
        term = None
        with open(self.pfam_file_path) as pfam_handle:
            for line in pfam_handle:
                if "DESC" in line:
                    assert ont_doc, line
                    assert term, line
                    ont_doc.description = line.split("DESC")[1].strip()
                    ont_doc.keywords = self.ki.extract_keywords([ont_doc.description, ont_doc.name, ont_doc.term]) + [
                        term]
                    ont_doc.save()
                    ont_doc = None
                    term = None
                elif "ACC" in line:
                    term = line.split("ACC")[1].strip().lower()
                    if ont_doc:
                        ont_doc.term = term
                    else:
                        name = term
                        ont_doc = Ontology(name=name, ontology=self.ontology_name, term=term)
                elif "NAME" in line:
                    name = line.split("NAME")[1].strip()

                    if ont_doc:
                        ont_doc.name = name
                    else:
                        term = line.split("NAME")[1].strip()
                        ont_doc = Ontology(name=name, ontology=self.ontology_name, term=term)

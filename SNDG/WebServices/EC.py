

class EC:
    enzclass = "ftp://ftp.expasy.org/databases/enzyme/enzclass.txt"
    ec2go = "http://current.geneontology.org/ontology/external2go/ec2go"

    @staticmethod
    def read_class(enzclass_path):
        term_desc = {}
        with open(enzclass_path) as h:
            for l in h:
                l = l.strip()
                if l and  l[1] == ".":
                    term = l[:9].replace(" ","")
                    desc = l[9:].strip()
                    term_desc[term] = desc
        return term_desc

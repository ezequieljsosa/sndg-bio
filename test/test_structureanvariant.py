import unittest

# from SNDG.Structure.StructureVariant import StructureVariant
from SNDG import execute
from io import StringIO
import tempfile
import json


class TestVcfSnpeffIO(unittest.TestCase):

    def test_residues_mapping(self):
        stdout = tempfile.NamedTemporaryFile()
        execute("python3 -m SNDG.Structure.StructureVariant residues -i ./test/prot2.fasta --pdb_data /tmp/data -s 2VYI"
                , stdout=stdout)

        with open(stdout.name) as output:
            contents = output.read()
        stdout.close()
        expected = ' '.join("""pdb     chain   resid   alt     ref     pos     pdb_pos
2vyi    A       81      P       G       23      1
2vyi    A       160     K       A       160     77
2vyi    B       160     K       A       160     74""".split())
        self.assertEqual(
            expected, ' '.join(contents.strip().split())
        )

    def test_residues_ann(self):
        stdout = tempfile.NamedTemporaryFile()
        execute("cat test/test.tbl | python -m SNDG.Structure.StructureVariant ann --pdb_data /tmp/data "
                , stdout=stdout)
        with open(stdout.name) as output:
            contents = output.read()
        stdout.close()
        anns = [json.loads(x) for x in contents.split("###") if x.strip()]
        self.assertEqual(4, len(anns))
        r2vyi_B_160_K = [x["ann"] for x in anns if x["residue"] == "2vyi_B_160_K"][0]

        self.assertEqual(0.705, r2vyi_B_160_K["pockets"][0]["druggabilitty"])

        r1azm_A_91_P = [x["ann"] for x in anns if x["residue"] == "1azm_A_91_P"][0]
        self.assertTrue("BINDING SITE FOR RESIDUE AZM A 262" in
                        [x["details"] for x in r1azm_A_91_P["binding"]])


if __name__ == '__main__':
    unittest.main()

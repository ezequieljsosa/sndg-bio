import logging
import requests
from collections import defaultdict

_log = logging.getLogger(__name__)


class PDBsWS:

    @staticmethod
    def is_obsolete(pdb: str):
        req = requests.get(f'https://data.rcsb.org/rest/v1/core/entry/{pdb}')
        if req.status_code == 404:
            return True
        elif req.status_code != 200:
            raise ConnectionError(f"error fetching https://data.rcsb.org/rest/v1/core/entry/{pdb}")
        return False

    @staticmethod
    def ligands(pdb: str):
        url = f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb}'
        req = requests.get(url)
        if req.status_code == 200:
            data = req.json()
            data_dict = defaultdict(list)
            for ligand in data[pdb]:
                data_dict[ligand["chem_comp_id"]].append(ligand)
            data_dict = dict(data_dict)
            return data_dict
        elif req.status_code == 404:
            return {}
        else:
            raise RuntimeError(f"error fetching {pdb} at {url}")

    @staticmethod
    def binding(pdb: str):
        """

        :param pdb:
        :return:
        {"1azm":[
        {
      "evidence_code": "Software",
      "ligand_residues": [
        {
          "entity_id": null,
          "residue_number": null,
          "author_insertion_code": null,
          "chain_id": null,
          "author_residue_number": null,
          "chem_comp_id": null,
          "struct_asym_id": null
        }
      ],
      "site_id": "AC2",
      "details": "BINDING SITE FOR RESIDUE AZM A 262",
      "site_residues": [
        {
          "entity_id": 1,
          "residue_number": 91,
          "author_insertion_code": null,
          "symmetry_symbol": "1_555",
          "chain_id": "A",
          "author_residue_number": 91,
          "chem_comp_id": "PHE",
          "struct_asym_id": "A"
        }, ...
        {
          "entity_id": 2,
          "residue_number": 1,
          "author_insertion_code": null,
          "symmetry_symbol": "1_555",
          "chain_id": "A",
          "author_residue_number": 261,
          "chem_comp_id": "ZN",
          "struct_asym_id": "B"
        }...
      ]
    },
        """
        req = requests.get(f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/{pdb}')
        if req.status_code == 200:
            data = req.json()
            return data[pdb]

        else:
            raise ConnectionError(f"error fetching {pdb}: {req.text}")

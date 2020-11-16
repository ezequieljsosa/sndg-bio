'''
Created on Aug 11, 2014

@author: eze
'''

SO_ROOT_TERMS = ["SO:0000400", "SO:0001260", "SO:0000110", "SO:0001060"]

EXCLUDED_SO = ["SO:0001080", "SO:0001128", "SO:0000408", "SO:0001068",
               "SO:0001114", "SO:0001111", "SO:0001089", "SO:0000001", "SO:0000419"];
SO_TERM_NAMES = {
    "SO:0000236": "ORF",
    "SO:0000704": "Gene",
    "SO:0000996": "predicted_gene",
    "SO:0001079": "polypeptide_structural_motif",
    "SO:0000857": "homologous",
    "SO:0000417": "polypeptide_domain",
    "SO:0000253": "tRNA",
    "SO:0000419": "mature_protein_region",
    "SO:0000419": "chain",
    "SO:0000839": "polypeptide_region",
    "SO:1000002": "substitution",
    "SO:0000409": "binding_site",
    "SO:0001104": "catalytic_residue",
    "SO:0001128": "polypeptide_turn_motif",
    "SO:0001077": "transmembrane_polypeptide_region",
    "SO:0001630": "splice_region_variant",
    "SO:0001656": "metal_binding_site",
    "SO:0001811": "phosphorylation_site",
    "SO:0001971": "zinc_finger_binding_site",
    "SO:0001429": "DNA_binding_site",
    "SO:0001114": "peptide_helix",
    "SO:0001068": "polypeptide_repeat",
    "SO:0001111": "beta_strand",
    "SO:0001089": "post_translationally_modified_region",
    "SO:0001088": "disulfide_bond",
    "SO:0001062": "propeptide",
    "SO:0001060": "sequence_variant",
    "SO:0001067": "polypeptide_motif",
    "SO:0000001": "region",
    "SO:0001066": "compositionally_biased_region_of_peptide",
    "SO:0000418": "signal_peptide",
    "SO:0000419": "mature_protein_region",
    "SO:0000691": "cleaved_initiator_methionine"
};
SO_TERMS = {v: k for k, v in SO_TERM_NAMES.items()}

SO_RESIDUES = {
    "asparagine": "SO:0001449",
    "Asn": "SO:0001449",
    "N": "SO:0001449",

    "alanine": "SO:0001435",
    "Ala": "SO:0001435",
    "A": "SO:0001435",

    "tryptophan": "SO:0001440",
    "Trp": "SO:0001440",
    "W": "SO:0001440",

    "selenocysteine": "SO:0001455",
    "Cys": "SO:0001455",
    "C": "SO:0001455",

    "tyrosine": "SO:0001446",
    "Tyr": "SO:0001446",
    "Y": "SO:0001446",

    "modified_amino_acid_feature": "SO:0001385",

    "proline": "SO:0001439",
    "Pro": "SO:0001439",
    "P": "SO:0001439",

    "histidine": "SO:0001452",
    "His": "SO:0001452",
    "H": "SO:0001452",

    "phenylalanine": "SO:0001441",
    "Phe": "SO:0001441",
    "F": "SO:0001441",

    "serine": "SO:0001444",
    "Ser": "SO:0001444",
    "S": "SO:0001444",

    "leucine": "SO:0001437",
    "Leu": "SO:0001437",
    "L": "SO:0001437",

    "pyrrolysine": "SO:0001456",

    "cysteine": "SO:0001447",
    "Cys": "SO:0001447",
    "C": "SO:0001447",

    "lysine": "SO:0001450",
    "Lys": "SO:0001450",
    "K": "SO:0001450",

    "arginine": "SO:0001451",
    "Arg": "SO:0001451",
    "R": "SO:0001451",

    "glutamine": "SO:0001448",
    "Gln": "SO:0001448",
    "Q": "SO:0001448",

    "glycine": "SO:0001443",
    "Gly": "SO:0001443",
    "G": "SO:0001443",

    "valine": "SO:0001436",
    "Val": "SO:0001436",
    "V": "SO:0001436",

    "glutamic_acid": "SO:0001454",
    "Glu": "SO:0001454",
    "E": "SO:0001454",

    "aspartic_acid": "SO:0001453",
    "Asp": "SO:0001453",
    "D": "SO:0001453",

    "threonine": "SO:0001445",
    "Thr": "SO:0001445",
    "T": "SO:0001445",

    "catalytic_residue": "SO:0001104",

    "methionine": "SO:0001442",
    "Met": "SO:0001442",
    "M": "SO:0001442",

    "isoleucine": "SO:0001438",
    "Ile": "SO:0001438",
    "I": "SO:0001438"
}

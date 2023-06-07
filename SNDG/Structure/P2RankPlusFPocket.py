if __name__ == "__main__":
    import os
    import argparse
    import subprocess as sp
    from tqdm import tqdm
    import pandas as pd
    from SNDG.Structure.FPocket import FpocketOutput
    from SNDG import mkdir
    import json

    fpocket_template = "docker run --user $UID:$GID --rm -i -v {workdir}:{workdir} ezequieljsosa/fpocket fpocket -f {pdbinput} -m 2 -D 6"
    p2rank_output = "/home/eze/workspace/pockets/p2rank_2.4/8IWM/8IWM.pdb_predictions.csv"
    pdb = "/home/eze/workspace/pockets/p2rank_2.4/8IWM.pdb"
    pocket_tmp = "/tmp/pocket/"
    output_path = "/home/eze/workspace/pockets/p2rank_2.4/8IWM/pockets.json"

    mkdir(pocket_tmp)

    df = pd.read_csv(p2rank_output)
    df.columns = [x.strip() for x in df.columns]
    with open(pdb) as h:
        pdb_lines = h.readlines()
    outpockets =[]
    for _, record in tqdm(list(df.iterrows())):

        residues = [x.strip() for x in record.residue_ids.split()]
        atoms = [x.strip() for x in record.surf_atom_ids.split()]
        pocket_name = record["name"].strip()
        pdb_pocket_file = pocket_tmp + pocket_name + ".pdb"

        with open(pdb_pocket_file, "w") as h:
            for x in pdb_lines:
                if x.startswith("ATOM"):
                    rid = x.split()[4].strip() + "_" + x.split()[5].strip()
                    if rid in residues:
                        h.write(x)
                else:
                    h.write(x)

        fpo = FpocketOutput(directory=pocket_tmp + pocket_name + "_out")

        if not os.path.exists(fpo._info_file_path()):
            cmd = fpocket_template.format(pdbinput=pdb_pocket_file, workdir=pocket_tmp)
            sp.call(cmd, shell=True)

        if os.path.exists(fpo._info_file_path()):
            fpo.parse()
            if fpo.pockets:
                pprops = fpo.pockets[0].properties
                pprops["name"] = pocket_name
                pprops["residues"] = residues
                pprops["atoms"] = atoms
                pprops["score"] = record["score"]
                pprops["probability"] = record["probability"]
                pprops["sas_points"] = record["sas_points"]
                pprops["alpha_spheres"] = fpo.pockets[0].alpha_spheres
                outpockets.append(pprops)
    with open(output_path,"w") as h:
        json.dump(outpockets,h)
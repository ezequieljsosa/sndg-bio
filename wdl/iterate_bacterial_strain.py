#!/usr/bin/env python
# encoding: utf-8

import sys
import re
import os
from tqdm import tqdm
from glob import glob
import json
from collections import defaultdict
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import subprocess as sp


def process_sample(sample_name, fastq1, fastq2, cpus, reference_filename, reference_directory,
                   organism_name, cromwell, wdl, output_dir):
    sample_dir = os.sep.join([output_dir, sample_name])
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)

    # Create input
    input_json = {
        "processBacterialStrain.strain": sample_name,
        "processBacterialStrain.fastq1": fastq1,
        "processBacterialStrain.fastq2": fastq2,
        "processBacterialStrain.reference_dir": reference_directory,
        "processBacterialStrain.reference_filename": reference_filename,
        "processBacterialStrain.species": organism_name,
        "processBacterialStrain.cpus": str(cpus)
    }
    input_json_path = os.sep.join([output_dir, sample_name + "_in.json"])
    with open(input_json_path, "w") as h:
        json.dump(input_json, h, indent=4, sort_keys=True)

    # call workflow
    stderr_file = os.sep.join([output_dir, "%s.err" % sample_name])
    stdout_file = os.sep.join([output_dir, "%s.out" % sample_name])
    with open(stderr_file,"w") as he:
        with open(stdout_file,"w") as ho:
            sp.run(["java", "-jar", cromwell, "run", wdl, "-i", input_json_path], stdout=ho,
                   stderr=he)

    # collect results
    with open(stdout_file) as h:
        outstr = h.read()
    output = outstr.split("Final Outputs:")[1].split("\n}\n")[0] + "}"

    output = json.loads(output)

    for k, outfile in output.items():
        task = k.split(".")[1]
        task_dir = os.sep.join([sample_dir, task])
        if not os.path.exists(task_dir):
            os.makedirs(task_dir)
        sp.run(["ln", outfile, task_dir],stdout=sys.stderr)


if __name__ == "__main__":
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument("-s", "--samples_glob", help="Samples glob expression", default="./*.fastq")
    parser.add_argument("-n", "--name_extraction", help="regexp to extract the sample name from the fastq filename",
                        default="(.+)_.*\.fastq")

    parser.add_argument("-rf", "--reference_filename", help="Reference fasta filename. NOT the PATH", required=True)
    parser.add_argument("-rd", "--reference_directory", help="Directory containing reference_filename and its indexes",
                        required=True)
    parser.add_argument("-on", "--organism_name", default='org', help="Name of the organism")
    parser.add_argument('--cpus', default=1, help="cpus to use")

    parser.add_argument('-c', '--cromwell', default=1, help="cromwell jar path")
    parser.add_argument('--wdl', default="./process_bacterial_strain.wdl", help="wdl path")

    parser.add_argument('-o', "--output_dir", default="./results", help="output directory")

    args = parser.parse_args()

    assert os.path.exists(args.reference_directory), "'%s' does not exists" % args.reference_directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    samples = defaultdict(lambda: [])

    sys.stderr.write("searching files...\n")
    for filepath in tqdm(glob(args.samples_glob), file=sys.stderr):
        filename = os.path.basename(filepath)
        sample_name = re.search(args.name_extraction, filename).group(1)
        samples[sample_name].append(filename)

    sys.stderr.write("processing samples...\n")
    pbar = tqdm(samples.items(), file=sys.stderr)
    for sample_name, fastqs in pbar:
        pbar.set_description(sample_name)
        if len(fastqs) == 2:
            process_sample(sample_name, fastqs[0], fastqs[1], args.cpus, args.reference_filename,
                           args.reference_directory,
                           args.organism_name, args.cromwell, args.wdl, args.output_dir)
        else:
            sys.stderr.write("%s has only one fastq\n" % sample_name)
    sys.stderr.write("finished!\n")

"""

https://github.com/tseemann/prokka
"""
from SNDG import docker_wrap_command, DOCKER_MAPPINGS, execute


class Prokka:
    """

    """
    DEFAULT_DOCKER_IMG = "staphb/prokka:latest"

    def __init__(self):
        pass

    def add_spades_data(self, fasta_path, gbk_path):
        pass

    def annotate(self, fasta_path, output, training=None,
                 locustag="PROKKA", gram=None, genus="", species="", strain="", kingdom="Bacteria",
                 gcode=0, rfam=False, increment=5, prefix="ann", cpus=1,centre=""):

        if not os.path.exists(os.path.abspath(output + "/../")):
            raise FileNotFoundError(f'{os.path.abspath(output + "/../")} not found, cant create output dir')
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f'{os.path.abspath(fasta_path)} not found')

        if training and not os.path.exists(training):
            raise FileNotFoundError(f'{os.path.abspath(training)} not found')

        db = f"--proteins {training}" if training else ""
        rfam = f"--rfam" if rfam else ""
        species = f"--species '{species}'" if species else ""
        strain = f"--strain '{strain}'" if strain else locustag
        kingdom = f"--kingdom '{kingdom}'" if kingdom else ""
        gram = f"--gram {gram}" if gram else ""
        centre = f"--centre '{centre}'" if centre else ""
        genus = f"--genus '{genus}'" if genus else ""


        cmd = f'''prokka --compliant --cpus {cpus} {gram} --addgenes {rfam} --locustag {locustag} --outdir {output} \
                  --prefix {prefix} --force {db} {fasta_path} {kingdom} {strain}  \
                  --increment {increment} --gcode {gcode}  {centre} {genus} {species}'''

        cmd2 = docker_wrap_command(cmd)
        execute(cmd2)


DOCKER_MAPPINGS["prokka"] = Prokka.DEFAULT_DOCKER_IMG

if __name__ == '__main__':
    import argparse
    import os

    parser = argparse.ArgumentParser(description='Mapping to variant calls pipeline.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input_fasta', required=True)
    required.add_argument('-o', '--output_dir', required=True)
    required.add_argument('-t', '--training_fasta', required=False)
    required.add_argument('-lt', '--locus_tag', required=True)
    required.add_argument('--rfam', action="store_true")
    required.add_argument('--centre', required=False)
    required.add_argument('--genus', required=False)
    required.add_argument('--species', required=False)
    required.add_argument('-strain', '--strain', required=False)

    required.add_argument("--cpus", default=1)

    args = parser.parse_args()

    prokka = Prokka()
    prokka.annotate(args.input_fasta, args.output_dir, rfam=args.rfam, training=args.training_fasta,cpus=args.cpus,
                    locustag=args.locus_tag,centre=args.centre,genus=args.genus,species=args.species)

import os
import sys

def reduce_seqs(first,second,idx):
    """
    TA T --> error
    TA A --> TA
    T T --> TT
    A TT --> ATT
    AA * --> AA
    """
    if not first:
        return second
    if len(first) == 1:
        return first + second
    if second == "*":
        return first
    if (len(first) > (idx + 1 )):
        assert (first[idx:] == second[idx:len(first)-1])
    return first[:idx+1] + second[len(first)- idx:]
    # return first[:idx] + second[:len(first)-1]



def fix_lines(lines, samples):

    if len(lines) == 1:
        return lines[0]
    final_ref = ""

    sample_gts1 = {s: "" for s in samples}
    sample_gts2 = {s: "" for s in samples}
    for idx_line, l in enumerate(lines):
        vec = l.split()
        ref = vec[3]
        # pos = int(vec[1])
        alts = {i: x for i, x in enumerate(vec[4].split(","))}

        gt_options = {k + 1: v for k, v in alts.items()}
        gt_options[0] = ref

        final_ref=reduce_seqs(final_ref,ref,idx_line)
        # ref = "".join([ref_tmp[i] for i in range(pos - i, pos)]) + ref

        format_f = vec[8]
        for idx, sample in enumerate(samples):
            gt_index = [i for i, x in enumerate(format_f.split(":")) if x == "GT"][0]
            gt = vec[9 + idx].split(":")[gt_index].replace("|","/")
            gt1,gt2 = gt.split("/")
            gt1,gt2 = int(gt1),int(gt2)
            sample_gts1[sample] = reduce_seqs(sample_gts1[sample],gt_options[gt1],idx_line)
            sample_gts2[sample] = reduce_seqs(sample_gts2[sample],gt_options[gt2],idx_line)
    alts = list(set(list(sample_gts1.values()) + list(sample_gts2.values())) - set([final_ref]) )
    ref_and_alts = [final_ref] + alts
    gts = []
    for sample in samples:
        alt1 = sample_gts1[sample]
        alt2 = sample_gts2[sample]
        alt = str(ref_and_alts.index(alt1)) + "/" + str(ref_and_alts.index(alt2))
        gts.append(alt)

    new_lines = lines[0].split("\t")
    new_lines[3] = final_ref
    new_lines[4] = ",".join(alts)

    for idx,s in enumerate(samples):
        new_lines[9 + idx] = gts[idx] + ":" + ":".join(new_lines[9 + idx].split(":")[1:])



    return "\t".join(new_lines)





if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='fix vcf errors')
    parser.add_argument('vcf_in')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-o', '--output',default=None)

    args = parser.parse_args()

    if args.verbose:
        os.environ["verbose"] = "y"
    if args.output:
        hw = open(args.output,"w")
    else:
        hw = sys.stdout

    try:
        if args.vcf_in != "-":
            h = open(args.vcf_in)
        else:
            h = sys.stdin
        prevlines = []
        prevpos = None
        for l in h:
            if l.startswith("#CHROM"):
                samples = l.split()[9:]
            if l.startswith("#"):
                hw.write(l)
            else:
                vec = l.split()
                pos = int(vec[1])
                if not prevlines:
                    prevlines.append(l)
                    prevpos = pos
                else:
                    if prevpos == (pos - 1):
                        prevlines.append(l)
                        prevpos = pos
                    else:
                        try:
                            fixed_lines = fix_lines(prevlines, samples)
                            hw.write(fixed_lines)
                        except:
                            for l in  prevlines:
                                hw.write(l)


                        prevlines = [l]
                        prevpos = pos
        hw.write(fix_lines(prevlines, samples))
    finally:
        if args.vcf_in != "-":
            h.close()
        hw.close()

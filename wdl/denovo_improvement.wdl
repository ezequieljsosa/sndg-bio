workflow denovoImprovement {
    String cpus
    File contigs
    File fastq1
    File fastq2
    File insert_size_metrics
    File bin_SSPACE
    File bin_GapFiller

    call scaffolding {
        input:
            fastq1=fastq1,
            fastq2=fastq2,
            insert_size_metrics=insert_size_metrics,
            incontigs=contigs,
            bin_SSPACE=bin_SSPACE  ,
	    cpus=cpus
    }
    call gap_filler {
            input:
                fastq1=fastq1,
                fastq2=fastq2,
                insert_size_metrics=insert_size_metrics,
                incontigs=scaffolding.contigs,
                bin_GapFiller=bin_GapFiller  ,
                cpus=cpus
        }
    call scaffolding2 {
            input:
                fastq1=fastq1,
                fastq2=fastq2,
                insert_size_metrics=insert_size_metrics,
                incontigs=gap_filler.contigs,
                bin_SSPACE=bin_SSPACE  ,
    	    cpus=cpus
        }

}

task scaffolding {
    File fastq1
    File fastq2
    File insert_size_metrics
    File incontigs
    File bin_SSPACE
    String cpus

    command {
        export PERL5LIB=$PERL5LIB:${bin_SSPACE}:${bin_SSPACE}/dotlib
        isize=`grep -v "#" ${insert_size_metrics} | head -n 3 | tail -n 1 | cut -f6`
        isize="$( cut -d '.' -f 1 <<< "$isize" )";
        orientation=`grep -v "#" ${insert_size_metrics} | head -n 3 | tail -n 1 | cut -f9`
        echo "lib bwa ${fastq1} ${fastq2} $isize 0.25 $orientation" > libraries.txt
        perl ${bin_SSPACE}/SSPACE_Standard_v3.0.pl -l libraries.txt -s ${incontigs} -x 1 -m 24 -k 3  -n 15 -g 1 -T ${cpus} -b sspace1
    }

    output {
            File contigs = "sspace1/sspace1.final.scaffolds.fasta"
    }

}

task scaffolding2 {
    File fastq1
    File fastq2
    File insert_size_metrics
    File incontigs
    File bin_SSPACE
    String cpus

    command {
        export PERL5LIB=$PERL5LIB:${bin_SSPACE}:${bin_SSPACE}/dotlib
        isize=`grep -v "#" ${insert_size_metrics} | head -n 3 | tail -n 1 | cut -f6`
        isize="$( cut -d '.' -f 1 <<< "$isize" )";
        orientation=`grep -v "#" ${insert_size_metrics} | head -n 3 | tail -n 1 | cut -f9`
        echo "lib bwa ${fastq1} ${fastq2} $isize 0.25 $orientation" > libraries.txt
        perl ${bin_SSPACE}/SSPACE_Standard_v3.0.pl -l libraries.txt -s ${incontigs} -x 1 -m 24 -k 3  -n 15 -g 1 -T ${cpus} -b sspace2
    }

    output {
            File contigs = "sspace2/sspace2.final.scaffolds.fasta"
    }

}

task gap_filler {

    File fastq1
    File fastq2
    File insert_size_metrics
    File incontigs
    File bin_GapFiller
    String cpus

    command {
        export PERL5LIB=$PERL5LIB:${bin_GapFiller}
        isize=`grep -v "#" ${insert_size_metrics} | head -n 3 | tail -n 1 | cut -f6`
        orientation=`grep -v "#" ${insert_size_metrics} | head -n 3 | tail -n 1 | cut -f9`
        echo "lib bwa ${fastq1} ${fastq2} $isize 0.25 $orientation" > libraries.txt
        perl ${bin_GapFiller}/GapFiller.pl -l libraries.txt -s ${incontigs} -m 20 -o 2 -g 1 -T ${cpus} -b gapfiller
    }
    output {
        File contigs = "gapfiller/gapfiller.gapfilled.final.fa"
    }

}




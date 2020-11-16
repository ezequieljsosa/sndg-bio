workflow processBacterialStrain {
    String cpus
    File reference_dir
    String reference_filename
    Array[File] vcfs

    Int? headcrop

    call trimming {
        input:
            fastq1=fastq1,
            fastq2=fastq2,
            cpus=cpus,
            headcrop=headcrop,
    }



    call alignment {
        input:
            reference_dir=reference_dir,
            reference_filename=reference_filename,
            strain=strain,
            species=species,
            cpus=cpus  ,
	    fastq1=trimming.trimmed1,
	    fastq2=trimming.trimmed2,
    }
    call variant_call {
      input:
        alignment_bam=alignment.mapped_bam,
        mapped_bam_idx=alignment.mapped_bam_idx,
    	reference_dir=reference_dir,
	reference_filename=reference_filename,
     }

     call genotype_gvcf {
      input:
         gvcf=variant_call.vcf,
         gvcf_idx=variant_call.vcfidx,
         reference_dir=reference_dir,
         reference_filename=reference_filename,
         ploidy=1,
     }

}

task denovo {
    File fastq1
    File fastq2
    String cpus

    command {
        spades.py --cov-cutoff 5 -t ${cpus} --pe1-1 ${fastq1} --pe1-2 ${fastq2} -o ./
    }
    runtime {
            docker: "ezequieljsosa/reads_processing"
    }
    output {
            File scaffolds = "scaffolds.fasta"
            File contigs = "contigs.fasta"
    }

}

task trimming {
    File fastq1
    File fastq2
    String cpus
    Int? headcrop

    command {
        java -jar $TRIMMOMATIC PE -phred33 -threads ${cpus} ${fastq1} ${fastq2} 1.fastq.gz  1U.fastq.gz 2.fastq.gz 2U.fastq.gz \
	ILLUMINACLIP:/app/Trimmomatic-0.38/adapters/TruSeq2-PE.fa:2:30:10  ILLUMINACLIP:/app/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:10 \
	ILLUMINACLIP:/app/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 ILLUMINACLIP:/app/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10  \
	HEADCROP:${default=18 headcrop}  TRAILING:3 SLIDINGWINDOW:10:5 MINLEN:36
	zcat 1U.fastq.gz  2U.fastq.gz > U.fastq
	gzip U.fastq
	rm 1U.fastq.gz  2U.fastq.gz
    }
    output {
        File trimmed1 = "1.fastq.gz"
	File trimmed2 = "2.fastq.gz"
	File singles = "U.fastq.gz"
    }
    runtime {
        docker: "ezequieljsosa/reads_processing"
    }
}

task alignment {
    String strain
    String species
    String cpus
    File reference_dir
    String reference_filename
    File fastq1
    File fastq2

    command {

        bwa mem -t ${cpus} -M -R '@RG\tID:group1\tSM:${strain}\tPL:illumina\\tLB:${species}' \
        ${reference_dir}/${reference_filename} ${fastq1} ${fastq2} > aligned_reads.sam

        samtools view -@ ${cpus} -S -b -h aligned_reads.sam > aligned_reads.bam
        samtools flagstat -@ ${cpus} aligned_reads.bam > flagstat.txt

        samtools view -@ ${cpus} -F 4 -S -b -h aligned_reads.bam > mapped_reads.bam
        samtools view -@ ${cpus} -f 4 -S -b -h aligned_reads.bam > unmapped_reads.bam
        bedtools bamtofastq -i unmapped_reads.bam -fq unmapped_1.fastq -fq2 unmapped_2.fastq


        java -jar $PICARD SortSam INPUT=mapped_reads.bam OUTPUT=sorted_reads.bam SORT_ORDER=coordinate
        rm mapped_reads.bam
        java -jar $PICARD MarkDuplicates INPUT=sorted_reads.bam OUTPUT=mapped_reads.bam METRICS_FILE=dedup.txt

        rm aligned_reads.bam unmapped_reads.bam aligned_reads.sam sorted_reads.bam
        java -jar $PICARD BuildBamIndex INPUT=mapped_reads.bam


    }
    output {
        File mapped_bam = "mapped_reads.bam"
        File mapped_bam_idx = "mapped_reads.bai"
        File unmapped_1_fastq = "unmapped_1.fastq"
	File unmapped_2_fastq = "unmapped_2.fastq"
	File dedup = "dedup.txt"
	File flagstat = "flagstat.txt"
    }
    runtime {
        docker: "ezequieljsosa/reads_processing"
    }
}

task variant_call {
    File alignment_bam
    File mapped_bam_idx
    File reference_dir
    String reference_filename
    String? ploidy

    command {
        java -jar /gatk/gatk-package-4.1.0.0-local.jar  HaplotypeCaller -ERC GVCF \
	    -R ${reference_dir}/${reference_filename} -ploidy ${default=1 ploidy} \
	    -I ${alignment_bam} \
	    -O raw.g.vcf.gz
    }
    output {
        File vcf = "raw.g.vcf.gz"
        File vcfidx = "raw.g.vcf.gz.tbi"
    }
    runtime {
        docker: "broadinstitute/gatk:4.1.0.0"
    }
}

task genotype_gvcf {

    File gvcf
    File gvcf_idx
    File reference_dir
    String reference_filename
    String? ploidy

    command {
        java -jar /gatk/gatk-package-4.1.0.0-local.jar GenotypeGVCFs  \
	    -R ${reference_dir}/${reference_filename} -ploidy ${default=1 ploidy} \
	    -V raw.g.vcf.gz  \
	    -O output.vcf.gz
    }
    output {
        File vcf = "output.vcf.gz"
        File vcfidx = "output.vcf.gz.tbi"
    }
    runtime {
        docker: "broadinstitute/gatk:4.1.0.0"
    }
}
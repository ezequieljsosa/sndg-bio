workflow processBacterialStrain {

    File gvcf
    File gvcf_idx
    File reference_dir
    String reference_filename
    Int? ploidy

    #TODO
    #https://software.broadinstitute.org/gatk/documentation/article?id=11065
    #https://gatkforums.broadinstitute.org/gatk/discussion/1247/what-should-i-use-as-known-variants-sites-for-running-tool-x
    #more Hard filtering https://software.broadinstitute.org/gatk/documentation/article.php?id=1255

#      call genotype_gvcf {
#              input:
#                 gvcf=gvcf,
#                 gvcf_idx=gvcf_idx,
#                 reference_dir=reference_dir,
#                 reference_filename=reference_filename,
#                 ploidy=1,
#             }


    # http://evodify.com/gatk-in-non-model-organism/
    call separate_snp_indel {
        input:
            gvcf=gvcf,#genotype_gvcf.vcf,
            gvcf_idx=gvcf_idx,#genotype_gvcf.vcfidx,
            reference_filename=reference_filename,
            reference_dir=reference_dir,

    }
    #https://software.broadinstitute.org/gatk/documentation/article?id=11097
    call variant_filtration {
            input:
                snp_gvcf=separate_snp_indel.snp_vcf,
                snp_gvcf_idx=separate_snp_indel.snp_vcfidx,
                indel_gvcf=separate_snp_indel.indel_vcf,
                indel_gvcf_idx=separate_snp_indel.indel_vcfidx,
                reference_filename=reference_filename,
                reference_dir=reference_dir,
        }
     call merge_vcfs {
            input:
                snp_gvcf=variant_filtration.snp_vcf,
                indel_gvcf=variant_filtration.indel_vcf,

        }



}

task separate_snp_indel {

    File gvcf
    File gvcf_idx
    File reference_dir
    String reference_filename
    String? ploidy

    command {
        java -jar /gatk/gatk-package-4.1.0.0-local.jar \
	     SelectVariants \
         -R ${reference_dir}/${reference_filename} \
         -V ${gvcf} \
         --select-type-to-include SNP \
         -O snp.vcf.gz

        java -jar /gatk/gatk-package-4.1.0.0-local.jar \
         SelectVariants \
         -R ${reference_dir}/${reference_filename} \
         -V ${gvcf} \
         --select-type-to-include INDEL \
         -O indel.vcf.gz

    }
    output {
        File snp_vcf = "snp.vcf.gz"
        File snp_vcfidx = "snp.vcf.gz.tbi"
        File indel_vcf = "indel.vcf.gz"
        File indel_vcfidx = "indel.vcf.gz.tbi"
    }
    runtime {
        docker: "broadinstitute/gatk:4.1.0.0"
    }
}


task variant_filtration {

    File snp_gvcf
    File snp_gvcf_idx
    File indel_gvcf
    File indel_gvcf_idx


    File reference_dir
    String reference_filename
    String? ploidy

    command {
        java -jar /gatk/gatk-package-4.1.0.0-local.jar \
	      VariantFiltration \
          -R ${reference_dir}/${reference_filename} \
          -V ${snp_gvcf} \
          --filter-expression "QUAL < 0 || MQ < 40.00 || SOR > 3.000 || QD < 2.00 || FS > 60.000 || MQRankSum < -12.5 || ReadPosRankSum < -8.000"  \
          --filter-name "sndg_filter1" \
          -O snp.vcf

        grep -E '^#|PASS' snp.vcf > snp2.vcf


        java -jar /gatk/gatk-package-4.1.0.0-local.jar \
              VariantFiltration \
              -R ${reference_dir}/${reference_filename} \
              -V ${indel_gvcf} \
              --filter-expression "QUAL < 0 || MQ < 30.00 || SOR > 10.000 || QD < 2.00 || FS > 200.000 || ReadPosRankSum < -20.000"  \
              --filter-name "sndg_filter2" \
              -O indel.vcf

        grep -E '^#|PASS' indel.vcf > indel2.vcf



    }
    output {
        File snp_vcf = "snp2.vcf"
        File indel_vcf = "indel2.vcf"

    }
    runtime {
        docker: "broadinstitute/gatk:4.1.0.0"
    }
}

task merge_vcfs {
    File snp_gvcf
    File indel_gvcf


    command {
    java -jar /app/picard.jar MergeVcfs \
            I=${snp_gvcf} \
            I=${indel_gvcf} \
            O=merged.vcf.gz
    }

    output {
            File vcf = "merged.vcf.gz"
            File vcfidx = "merged.vcf.gz.tbi"
        }
        runtime {
            docker: "ezequieljsosa/reads_processing"
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
	    -V ${gvcf}  \
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

#task snpEff {
#
#    String db
#    File vcf
#
#    command {
#        java -jar /app/snpEff/snpEff.jar  ann ${db} ${vcf} > output.ann.vcf
#    }
#    output {
#        File vcf = "output.ann.vcf"
#    }
#    runtime {
#        docker: "ezequieljsosa/reads_processing"
#    }
#}
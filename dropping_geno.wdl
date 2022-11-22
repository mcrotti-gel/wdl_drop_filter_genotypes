version 1.0

workflow remove_genotypes {
    input {
        File input_vcf
        Array[String] vcf_list=read_lines(input_vcf)
        }
    
    scatter(vcf in vcf_list) {
        String vcf_name = basename(vcf)
        call drop_genotypes {
            input: 
                vcf=vcf, 
                vcf_name=vcf_name
        }

        call medianGQ_filter {
            input:
                drop_genotypes_output=drop_genotypes.drop_genotypes_output,
                vcf_name=vcf_name
        }
    }

    output {
        Array[File] medianGQ_file = medianGQ_filter.medianGQ_filter_output
    }
}

task drop_genotypes {
    
    input {
        File vcf
        String vcf_name
    }
    
    runtime {
        cpu: 1
        memory: "2 GB"
        docker: "082963657711.dkr.ecr.eu-west-2.amazonaws.com/bcftools:1.13"
    }

    command {

        bcftools view -G -O z -o ${vcf_name}_sites.vcf.gz ${vcf} 
        bcftools index ${vcf_name}_sites.vcf.gz
    }

    output {
        File drop_genotypes_output = "${vcf_name}_sites.vcf.gz"
        File drop_genotypes_output_index = "${vcf_name}_sites.vcf.gz.csi"
    }
}

task medianGQ_filter {

    input {
        File drop_genotypes_output
        String vcf_name
    }

    runtime {
        cpu: 1
        memory: "2 GB"
        docker: "082963657711.dkr.ecr.eu-west-2.amazonaws.com/bcftools:1.13"
    }

    command {
        bcftools query -i 'medianGQ > 30' -f '%CHROM\t%POS\t%REF\t%ALT\t%medianGQ\n' ${drop_genotypes_output} > ${vcf_name}_medianGQ.tsv
    }

    output {
        File medianGQ_filter_output = "${vcf_name}_medianGQ.tsv"
    }

}
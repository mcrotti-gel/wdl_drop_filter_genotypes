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

        call minGQ_filter {
            input:
                drop_genotypes_output=drop_genotypes.drop_genotypes_output,
                vcf_name=vcf_name
        }
    }

    output {
        Array[File] minGQ_file = minGQ_filter.minGQ_filter_output
    }
}

task drop_genotypes {
    
    input {
        File vcf
        String vcf_name
    }
    
    runtime {
        cpus: 1
        memory: "2 GB"
        docker: "quay.io/lifebitai/bcftools"
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

task minGQ_filter {

    input {
        File drop_genotypes_output
        String vcf_name
    }

    runtime {
        cpus: 1
        memory: "2 GB"
        docker: "quay.io/lifebitai/bcftools"
    }

    command {
        bcftools query -i 'medianGQ < 80' -f '[%CHROM\t%POS\t%REF\t%ALT\t%medianGQ\n]' ${drop_genotypes_output} > ${vcf_name}_minGQ.tsv
    }

    output {
        File minGQ_filter_output = "${vcf_name}_minGQ.tsv"
    }

}
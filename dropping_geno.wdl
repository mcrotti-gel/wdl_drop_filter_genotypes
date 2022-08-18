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
    }

    output {
        Array[File] drop_geno_vcf = drop_genotypes.drop_genotypes_output
        Array[File] drop_geno_vcf_index = drop_genotypes.drop_genotypes_output_index
    }
}

task drop_genotypes {
    
    input {
        File vcf
        String vcf_name
    }
    

    runtime {
        cpus: 1
        memory: "10 GB"
        docker: "quay.io/lifebitai/bcftools"
    }

    command {

        bcftools view -G -O z -o ${vcf_name}_sites.vcf.gz ${vcf} 
        bcftools index -t ${vcf_name}_sites.vcf.gz
    }

    output {
        File drop_genotypes_output = "${vcf_name}_sites.vcf.gz"
        File drop_genotypes_output_index = "${vcf_name}_sites.vcf.gz.tbi"
    }
}


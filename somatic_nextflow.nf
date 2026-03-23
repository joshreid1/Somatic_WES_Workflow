#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Input parameters
params.output_prefix       = './somatic_nextflow'
params.sample_info         = "${projectDir}/pipeline_files/manifests/test_manifest.tsv"
params.ref_fasta           = '/stornext/Bioinf/data/lab_bahlo/projects/epilepsy/hg38/reference/fasta/Homo_sapiens_assembly38.fasta'
params.gene_bed            = '/vast/projects/reidj-project/gene_lists/Somatic_Gene_List_v2025-06.bed'
params.spliceai_distance   =  500

// vep files
params.vep_cache_dir		= '/vast/projects/bahlo_cache/vep_cache/'
params.vep_alphamissense 	= '/vast/projects/bahlo_cache/annotation/alphamissense/AlphaMissense_hg38.tsv.gz'
params.vep_revel			= '/vast/projects/bahlo_cache/annotation/REVEL/revel_1.3.hg38.vep.tsv.gz'

// vcfanno toml files
params.clinvar_toml         		= "${projectDir}/pipeline_files/vcfanno_files/clinvar_20250330.toml"
params.cosmic_toml          		= "${projectDir}/pipeline_files/vcfanno_files/cosmic_20251203.toml"
params.gnomad_toml        			= "${projectDir}/pipeline_files/vcfanno_files/gnomad_v4.0.0.toml"
params.gnomad_postprocess_toml 		= "${projectDir}/pipeline_files/vcfanno_files/gnomad_v4.0.0_postprocess.toml"
params.spliceai_lua					= "${projectDir}/pipeline_files/vcfanno_files/spliceai.lua"
params.spliceai_postprocess_toml 	= "${projectDir}/pipeline_files/vcfanno_files/spliceai_postprocess.toml"

process ProcessSample {
    cpus 1
    memory { 5 * task.attempt + ' GB' }
    time { 1 * task.attempt + ' h'}

    input:
    val(meta)
    
    output:
    tuple val(meta), path("merged.vcf"), emit: MergedVCF
    tuple val(meta), env(BAM_FILE), emit: BamFILE
    tuple val(meta), 
          path('mutect_v4.vcf.gz*', arity: '2'), 
          path('strelka_v5.vcf.gz*', arity: '2'), 
          path('freebayes_v4.vcf.gz*', arity: '2'),
          env(HAS_STRELKA),
          emit: IndividualTools

    script:
    def SUB_DIR = "${params.output_prefix}/${meta.id}"
    
    // Determine BAM file path
    def bam_path = "${meta.bam_path}/${meta.tumor}/${meta.tumor}.recal.bam"
    def cram_path = "${meta.bam_path}/${meta.tumor}/${meta.tumor}.recal.cram"
    def BAM_FILE = ""
   
    if (file(bam_path).exists()) {
        BAM_FILE = bam_path
    } else if (file(cram_path).exists()) {
        BAM_FILE = cram_path
    } else {
        error "No .bam or .cram file found for ${meta.tumor}"
    }

    // Determine VCF paths based on sample type
    def MUTECT_VCF = ""
    def STRELKA_VCF = ""
    def STRELKA_SNV_VCF = ""
    def STRELKA_INDEL_VCF = ""
    def FREEBAYES_ORIGINAL_VCF = ""
    def IS_TUMOR_ONLY = false

    if (meta.normal == "NA" && meta.tumor) {
        log.info "Processing tumor-only sample: ID=${meta.id}, TUMOR=${meta.tumor}"
        IS_TUMOR_ONLY = true
        MUTECT_VCF = "${meta.vcf_path}/mutect2/${meta.tumor}/${meta.tumor}.mutect2.filtered.vcf.gz"
        STRELKA_VCF = "NA"
        STRELKA_SNV_VCF = "NA"
        STRELKA_INDEL_VCF = "NA"
        FREEBAYES_ORIGINAL_VCF = "${meta.vcf_path}/freebayes/${meta.tumor}/${meta.tumor}.freebayes.vcf.gz"
    } else if (meta.normal && meta.normal != "NA" && meta.tumor && meta.tumor != "NA") {
        log.info "Processing paired tumor/normal sample: ID=${meta.id}, NORMAL=${meta.normal}, TUMOR=${meta.tumor}"
        IS_TUMOR_ONLY = false
        MUTECT_VCF = "${meta.vcf_path}/mutect2/${meta.id}/${meta.id}.mutect2.filtered.vcf.gz"
        STRELKA_SNV_VCF = "${meta.vcf_path}/strelka/${meta.tumor}_vs_${meta.normal}/${meta.tumor}_vs_${meta.normal}.strelka.somatic_snvs.vcf.gz"
        STRELKA_INDEL_VCF = "${meta.vcf_path}/strelka/${meta.tumor}_vs_${meta.normal}/${meta.tumor}_vs_${meta.normal}.strelka.somatic_indels.vcf.gz"
        STRELKA_VCF = "NA"
        FREEBAYES_ORIGINAL_VCF = "${meta.vcf_path}/freebayes/${meta.tumor}_vs_${meta.normal}/${meta.tumor}_vs_${meta.normal}.freebayes.vcf.gz"
    }

    """
    # Export BAM_FILE for output
    export BAM_FILE="${BAM_FILE}"
    
    # Set HAS_STRELKA flag for downstream processes
    if [ "${IS_TUMOR_ONLY}" == "true" ]; then
        export HAS_STRELKA="false"
    else
        export HAS_STRELKA="true"
    fi
    
    # Start of bash script
    echo "Running variant processing for sample ID: ${meta.id}"
    echo "  NORMAL: ${meta.normal}"
    echo "  TUMOR: ${meta.tumor}"
    echo "  bam_file: ${BAM_FILE}"
    echo "  IS_TUMOR_ONLY: ${IS_TUMOR_ONLY}"
    echo "  HAS_STRELKA: \${HAS_STRELKA}"

    # Print input paths
    echo "  MUTECT_VCF: ${MUTECT_VCF}"
    echo "  STRELKA_VCF: ${STRELKA_VCF}"
    echo "  STRELKA_SNV_VCF: ${STRELKA_SNV_VCF}"
    echo "  STRELKA_INDEL_VCF: ${STRELKA_INDEL_VCF}"
    echo "  FREEBAYES_ORIGINAL_VCF: ${FREEBAYES_ORIGINAL_VCF}"

    # Handle Strelka processing based on sample type
    if [ "${IS_TUMOR_ONLY}" == "true" ]; then
        echo "Tumor-only sample detected - creating empty Strelka placeholder VCF"
        
        # Create a minimal empty VCF with proper header for Strelka
        cat <<'EOF' > strelka_v2.vcf
##fileformat=VCFv4.2
##source=strelka_placeholder_tumor_only
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=Strelka,Number=1,Type=Float,Description="Strelka caller placeholder">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
EOF
        bgzip -c strelka_v2.vcf > strelka_v2.vcf.gz
        tabix -f strelka_v2.vcf.gz
        
    elif [ "${meta.normal}" != "NA" ]; then
        # Paired tumor/normal processing
        python /vast/projects/reidj-project/software/VCFpytools/add_vaf_strelka2.py --input ${STRELKA_SNV_VCF} --output strelka_snv_v1.vcf.gz --variant snv
        python /vast/projects/reidj-project/software/VCFpytools/add_vaf_strelka2.py --input ${STRELKA_INDEL_VCF} --output strelka_indel_v1.vcf.gz --variant indel

        # Split multiallelic variants for each caller
        bcftools norm -m -any -f ${params.ref_fasta}  --output-type z -o strelka_snv_v2.vcf.gz strelka_snv_v1.vcf.gz
        bcftools norm -m -any -f ${params.ref_fasta}  --output-type z -o strelka_indel_v2.vcf.gz strelka_indel_v1.vcf.gz

        tabix -f strelka_snv_v2.vcf.gz
        tabix -f strelka_indel_v2.vcf.gz
        
        # Join Strelka calls
        bcftools concat -a -O z -o strelka_v2.vcf.gz strelka_snv_v2.vcf.gz strelka_indel_v2.vcf.gz
        tabix -f strelka_v2.vcf.gz
        STRELKA_VCF="strelka_v2.vcf.gz"

        NORMAL_INDEX_MUTECT=\$(bcftools query -l ${MUTECT_VCF} | grep -n ${meta.normal} | cut -d: -f1 | awk '{print \$1-1}')
        NORMAL_INDEX_FREEBAYES=\$(bcftools query -l ${FREEBAYES_ORIGINAL_VCF} | grep -n ${meta.normal} | cut -d: -f1 | awk '{print \$1-1}')
        NORMAL_INDEX_STRELKA=\$(bcftools query -l strelka_v2.vcf.gz | grep -n NORMAL | cut -d: -f1 | awk '{print \$1-1}')
    fi

    TUMOR_INDEX_MUTECT=\$(bcftools query -l ${MUTECT_VCF} | grep -n ${meta.tumor} | cut -d: -f1 | awk '{print \$1-1}')
    TUMOR_INDEX_FREEBAYES=\$(bcftools query -l ${FREEBAYES_ORIGINAL_VCF} | grep -n ${meta.tumor} | cut -d: -f1 | awk '{print \$1-1}')
    
    # Only query Strelka index if not tumor-only
    if [ "${IS_TUMOR_ONLY}" != "true" ]; then
        TUMOR_INDEX_STRELKA=\$(bcftools query -l strelka_v2.vcf.gz | grep -n TUMOR | cut -d: -f1 | awk '{print \$1-1}')
    fi

    ## Fix freebayes "FORMAT=GQ" line ##
    bgzip -d -c ${FREEBAYES_ORIGINAL_VCF} | sed 's/##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">/##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">/' | bgzip -f -c > freebayes_${meta.id}.vcf.gz
    tabix -f freebayes_${meta.id}.vcf.gz

    # Split multiallelic variants for each caller
    bcftools norm -m -any -f ${params.ref_fasta} --threads 4 --output-type z -o "mutect_v1.vcf.gz" ${MUTECT_VCF}
    bcftools norm -m -any -f ${params.ref_fasta}  --threads 4 --output-type z -o "freebayes_v1.vcf.gz" "freebayes_${meta.id}.vcf.gz"

    tabix -f mutect_v1.vcf.gz
    tabix -f freebayes_v1.vcf.gz

    bcftools annotate --remove INFO,FILTER,^FORMAT/AD "strelka_v2.vcf.gz" -Oz -o "strelka_v3.vcf.gz"
    bcftools annotate --remove INFO,FILTER,^FORMAT/AD "mutect_v1.vcf.gz" -Oz -o "mutect_v2.vcf.gz"
    bcftools annotate --remove INFO,FILTER,^FORMAT/AD "freebayes_v1.vcf.gz" -Oz -o "freebayes_v2.vcf.gz"

    if [ "${meta.normal}" == "NA" ] && [ -n "${meta.tumor}" ]; then
        # Tumor-only: Filter for variants present in TUMOR for mutect and freebayes only
        bcftools view --threads 4 -i "AD[0:1]>=1" --output-type z -o mutect_v3.vcf.gz mutect_v2.vcf.gz
        bcftools view --threads 4 -i "AD[0:1]>=1" --output-type z -o freebayes_v3.vcf.gz freebayes_v2.vcf.gz
        
        # For Strelka, just copy the empty placeholder through
        cp strelka_v3.vcf.gz strelka_v4.vcf.gz

    elif [[ -n "${meta.normal}" && "${meta.normal}" != "NA" && -n "${meta.tumor}" ]]; then
        # Paired tumor/normal processing
        FULL_NORMAL="${meta.id}_${meta.normal}"
        FULL_TUMOR="${meta.id}_${meta.tumor}"

        # build sample_order once
        printf "%s\n%s\n" "\$FULL_NORMAL" "\$FULL_TUMOR" > sample_order.txt

        # callers: name raw-vcf tumor-index-var normal-index-var out-vcf-suffix
        callers=(
            "mutect:mutect_v2.vcf.gz:TUMOR_INDEX_MUTECT:NORMAL_INDEX_MUTECT:mutect_v3.vcf.gz"
            "freebayes:freebayes_v2.vcf.gz:TUMOR_INDEX_FREEBAYES:NORMAL_INDEX_FREEBAYES:freebayes_v3.vcf.gz"
            "strelka:strelka_v3.vcf.gz:TUMOR_INDEX_STRELKA:NORMAL_INDEX_STRELKA:strelka_v4.vcf.gz"
        )

        for entry in "\${callers[@]}"; do
            IFS=":" read -r caller rawvcf BIDX LIDX outvcf <<< "\$entry"

            # skip missing
            if [[ ! -s "\$rawvcf" ]]; then
                echo "WARNING: \$rawvcf not found—skipping \$caller" >&2
                continue
            fi

            # detect old sample names in header
            samples=\$(bcftools query -l "\$rawvcf")
            # for Strelka we also check NORMAL/TUMOR literal
            if [[ "\$caller" == "strelka" && \$(grep -xc "NORMAL" <<<"\$samples") -gt 0 ]]; then
                OLD_B="NORMAL"; OLD_R="TUMOR"
            else
                # generic normal
                if grep -qx "${meta.normal}" <<<"\$samples"; then
                    OLD_B="${meta.normal}"
                else
                    OLD_B="\$FULL_NORMAL"
                fi
                # generic tumor
                if grep -qx "${meta.tumor}" <<<"\$samples"; then
                    OLD_R="${meta.tumor}"
                else
                    OLD_R="\$FULL_TUMOR"
                fi
            fi

            # only reheader if names differ
            if [[ "\$OLD_B" != "\$FULL_NORMAL" || "\$OLD_R" != "\$FULL_TUMOR" ]]; then
                printf "%s\t%s\n%s\t%s\n" \
                    "\$OLD_B" "\$FULL_NORMAL" \
                    "\$OLD_R" "\$FULL_TUMOR" \
                > reheader.txt

                bcftools reheader -s reheader.txt "\$rawvcf" -o tmp_reheader.vcf.gz
                tabix tmp_reheader.vcf.gz

                invcf="tmp_reheader.vcf.gz"
            else
                invcf="\$rawvcf"
            fi

            # now filter tumor>1 & normal<1
            bcftools view --threads 4 \
                -i "AD[\${!BIDX}:1]>1 && AD[\${!LIDX}:1]<1" \
                "\$invcf" \
                | bcftools view -S sample_order.txt -Oz -o "\$outvcf"

            # clean up
            [[ -f tmp_reheader.vcf.gz ]] && rm tmp_reheader.vcf.gz tmp_reheader.vcf.gz.tbi
        done
    fi

    tabix -f mutect_v3.vcf.gz
    tabix -f freebayes_v3.vcf.gz
    tabix -f strelka_v4.vcf.gz

    # Add variant caller INFO field to each VCF
    bcftools filter -e '1==1' -s "1.0" mutect_v3.vcf.gz | bcftools annotate -c INFO/Mutect:=FILTER -o mutect_v4.vcf
    sed -i '/##INFO=<ID=Mutect/ s/Type=String/Type=Float/' mutect_v4.vcf
    bgzip -f mutect_v4.vcf
    tabix -f mutect_v4.vcf.gz

    bcftools filter -e '1==1' -s "1.0" strelka_v4.vcf.gz | bcftools annotate -c INFO/Strelka:=FILTER -o strelka_v5.vcf
    sed -i '/##INFO=<ID=Strelka/ s/Type=String/Type=Float/' strelka_v5.vcf
    bgzip -f strelka_v5.vcf
    tabix -f strelka_v5.vcf.gz

    bcftools filter -e '1==1' -s "1.0" freebayes_v3.vcf.gz | bcftools annotate -c INFO/Freebayes:=FILTER -o freebayes_v4.vcf
    sed -i '/##INFO=<ID=Freebayes/ s/Type=String/Type=Float/' freebayes_v4.vcf
    bgzip -f freebayes_v4.vcf
    tabix -f freebayes_v4.vcf.gz

    # Conditional merge: only include Strelka if not tumor-only
    if [ "${IS_TUMOR_ONLY}" == "true" ]; then
        # For tumor-only, merge without Strelka (it's empty anyway)
        bcftools concat -a -D mutect_v4.vcf.gz freebayes_v4.vcf.gz -Oz -o merged.vcf.gz
    else
        # For paired samples, include all three callers
        bcftools concat -a -D mutect_v4.vcf.gz strelka_v5.vcf.gz freebayes_v4.vcf.gz -Oz -o merged.vcf.gz
    fi
    
    tabix -f merged.vcf.gz

    bcftools annotate --remove FILTER,INFO merged.vcf.gz > merged.vcf

    echo "Finished processing sample ID: ${meta.id}"
    """
}

process Split_Vcf {

    cpus = 1
    memory = { 5 * task.attempt + ' GB' }
    time   = { 1 * task.attempt + ' h' }

    container 'community.wave.seqera.io/library/snpsift_gzip:e2aa3af43b84c914'

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("merged.*.vcf")

    shell:
    '''
    variant_number=$(zgrep -v "^#" !{vcf} | wc -l)

    if (( variant_number > 1000 )); then
        task_number=500
    elif (( variant_number > 10 )); then
        task_number=10
    else
        task_number=2
    fi

    split_number=$(( variant_number / task_number ))
    (( split_number < 5 )) && split_number=5

    echo "Variant number = $variant_number"            >> error_check.txt
    echo "Number of tasks = $task_number"              >> error_check.txt
    echo "Number of variants per split = $split_number" >> error_check.txt

    /opt/conda/share/snpsift-5.4.0c-0/scripts/snpSift split -l $split_number !{vcf}
    '''
}



process mpileup_check {

    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'

    tag "${meta.id}"

    cpus = 1
    memory = { 5 * task.attempt + ' GB' }
    time = { 1 * task.attempt + ' h' }

    input:
        tuple val(meta), path(split_vcf), path(BamFILE)

    output:
        tuple val(meta), path("mpileup_version.vcf.gz")

    script:
        """
        bcftools mpileup -x --indels-2.0 --max-depth 1000000 --max-idepth 100000 \
            --gap-frac 0.0001 --min-BQ 20 --min-MQ 1 \
            --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR \
            --fasta-ref ${params.ref_fasta} \
            --targets-file ${split_vcf} ${BamFILE} |
        bcftools norm --multiallelics -any -f ${params.ref_fasta} |
        bcftools call -mA |
        bcftools view -i "FORMAT/ADF[0:1]>=1 && FORMAT/ADR[0:1]>=1 && FORMAT/AD[0:1]>=3" \
            -Oz -o mpileup_version.vcf.gz
        """
}


process ensembl_vep {

    container = 'quay.io/biocontainers/ensembl-vep:115.2--pl5321h2a3209d_1'

    tag "${meta.id}"

    cpus = 1
    memory = { 5 * task.attempt + ' GB' }
    time = { 1 * task.attempt + ' h' }

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("merged_vep.vcf")

    script:
        """
        echo "Running ensembl-vep..."
        vep --cache --dir ${params.vep_cache_dir} --cache_version 115 --assembly GRCh38 \
            -i ${vcf} -o merged_vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs \
            --fasta ${params.ref_fasta} --offline --sift b --polyphen b --ccds --hgvs --hgvsg \
            --symbol --numbers --protein --af --af_1kg --max_af --variant_class \
            --pick_allele_gene --force_overwrite
        echo "ensembl-vep complete."
        """
}

process gnomad {

    tag "${meta.id}"

    cpus = 1
    memory = { 5 * task.attempt + ' GB' }
    time = { 1 * task.attempt + ' h' }

    container 'quay.io/biocontainers/vcfanno:0.2.6--0'

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("somatic_variants_annotated.vcf")

    script:
        """
        vcfanno -p 4 ${params.gnomad_toml} ${vcf} > somatic_variants_vep_gnomad.vcf
        vcfanno -p 4 ${params.gnomad_postprocess_toml} somatic_variants_vep_gnomad.vcf > somatic_variants_annotated.vcf
        rm somatic_variants_vep_gnomad.vcf
        """
}

process CADD_Run_Container {

    tag "${meta.id}"

	cpus = 1
	memory = { 10 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	container 'oras://docker.io/joshreid1/cadd-scoring:v1.6_edit'

	containerOptions '-B /vast/projects/bahlo_epilepsy/somatic_annotation_data/CADD-scripts/data/annotations:/CADD-scripts/data/annotations --writable-tmpfs'

	input: 
		tuple val(meta), path(vcf) 
		
	output:
		val(meta)
		path(vcf)
		path("tmp.tsv.gz")

	shell:
	'''
	#Link local software
	ln -s /stornext/System/data/software/rhel/9/base/tools/snakemake/8.11.3/bin/snakemake /usr/local/bin/snakemake
	ln -s /stornext/System/data/software/rhel/9/base/bioinf/bcftools/1.20/bin/bcftools /usr/local/bin/bcftools

	if [[ $(bcftools query -f '%ALT\n' !{vcf} | uniq) == "*" ]]; then
		cp !{vcf} tmp.vcf
		continue
	else
		grep -v "^#" !{vcf} | cut -f 1-5 | sed 's/^chr//' > tmp.vcf
		/CADD-scripts/CADD.sh tmp.vcf
	fi
	'''
}

process Process_CADD {
	cpus = 1
	memory = { 10 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

    container 'community.wave.seqera.io/library/htslib_vcfanno:8044b99f5458cd69'
		
	input:
		val(meta)
		path(vcf)
		path(cadd_tsv) 
		
	output:
		val(meta)
		path("*.cadd_run.vcf")
	
	shell:
	'''
	output=$(basename $PWD)
	
	cp "$(readlink -f !{cadd_tsv})" ./cadd_tmp

	tabix -f -b 2 -e 2 -s 1 ./cadd_tmp

	cadd_path=$(realpath ./cadd_tmp)
	echo '[[annotation]]' > cadd.toml
	echo "file= '${cadd_path}'" >> cadd.toml
	echo 'names=["CADD_Score"]' >> cadd.toml
	echo 'ops=["mean"]' >> cadd.toml
	echo 'columns=[6]' >> cadd.toml

	vcfanno cadd.toml !{vcf} > $output.cadd_run.vcf 
	'''
}

process ClinVar {
	label 'C1M1T1'

    container 'quay.io/biocontainers/vcfanno:0.2.6--0'

	input:
		val(meta)
        path(vcf)

    output:
    tuple val(meta), path("clinvar.vcf"), emit: clinvar
				
	shell:
	'''
  	vcfanno !{params.clinvar_toml} !{vcf} > clinvar.vcf
	'''
}

process SpliceAI_Run {
    tag "${meta.id}"

	cpus = 1
	memory = { 16 * task.attempt + ' GB' }
	time = { 2 * task.attempt + ' h'}

	container 'community.wave.seqera.io/library/python_pip_keras_setuptools_pruned:1c71801b2a7b49db'

	input:
		val(meta)
		path(vcf) 

	output:
		tuple val(meta), path("spliceai.vcf"), emit: vcf

	script:
	""" 
	export PYTHON_EGG_CACHE=./
	spliceai -I ${vcf} \
	-O spliceai.vcf \
	-R ${params.ref_fasta} \
	-D ${params.spliceai_distance} \
	-A /vast/projects/bahlo_epilepsy/somatic_annotation_data/SpliceAI/gencode.v43.canonical.annotation.txt
	"""
}
 
process Join_VCF {

    tag "${meta.id}"

	cpus = 2
	memory = { 10 * task.attempt + ' GB' }
	time = { 2 * task.attempt + ' h'}

	container 'community.wave.seqera.io/library/bcftools_snpsift:6e1309d39c8edffd'

	input:
		tuple val(meta), path('vcf', arity: '1..*')

	output:
		tuple val(meta), path("*.gz"), path("*.tbi"), emit: JoinedVCF

	shell:
	'''
	if [ $(ls vcf* 2>/dev/null | wc -l) -gt 1 ]; then
		/opt/conda/share/snpsift-5.4.0c-0/scripts/snpSift split -j vcf* | bcftools sort --temp-dir ./ -Oz -o final.vcf.gz --write-index=tbi
	else
		bcftools sort --temp-dir ./ -Oz -o final.vcf.gz --write-index=tbi vcf*
	fi
	'''	
}


process Check_Tools {

    tag "${meta.id}"    
    
    cpus = 1
    memory = { 1 * task.attempt + ' GB' }
    time = { 1 * task.attempt + ' h'}

    container 'quay.io/biocontainers/vcfanno:0.2.6--0'

    publishDir "results", mode: "copy"

    input:
        tuple val(meta), path(vcf), path(vcf_index), path(mutect), path(freebayes), path(strelka), val(has_strelka)
        
    output:
        tuple val(meta), path("*_somatic_variants_merged.vcf")
    
    shell:
    '''
    # Determine if Strelka should be included based on has_strelka flag
    if [ "!{has_strelka}" == "true" ]; then
        echo "Including Strelka annotations for paired tumor/normal sample"
        
        # Step 1: Generate the TOML configuration file for VCFAnno with Strelka
        cat <<EOL > "annotations_config.toml"
[[annotation]]
file="mutect_v4.vcf.gz"
fields=["Mutect"]
ops=["self"] 

[[annotation]]
file="strelka_v5.vcf.gz"
fields=["Strelka"]
ops=["self"]

[[annotation]]
file="freebayes_v4.vcf.gz"
fields=["Freebayes"]
ops=["self"]

[[postannotation]]
fields=["Mutect", "Strelka", "Freebayes"]
op="sum"
name="Tool_Count"
type="Float"
EOL

    else
        echo "Excluding Strelka annotations for tumor-only sample"
        
        # Step 1: Generate the TOML configuration file for VCFAnno without Strelka
        cat <<EOL > "annotations_config.toml"
[[annotation]]
file="mutect_v4.vcf.gz"
fields=["Mutect"]
ops=["self"] 

[[annotation]]
file="freebayes_v4.vcf.gz"
fields=["Freebayes"]
ops=["self"]

[[postannotation]]
fields=["Mutect", "Freebayes"]
op="sum"
name="Tool_Count"
type="Float"
EOL

    fi

    # Step 2: Apply annotations using VCFAnno
    vcfanno -p 4 annotations_config.toml !{vcf} > !{meta.id}_somatic_variants_merged.vcf
    '''
}

process Compress_Index {

    tag "${meta.id}"	

	cpus = 1
	memory = { 1 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	input:
        tuple val(meta), path(candidate_variants)

	output:
		tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi")

	script:
	"""
	bgzip ${candidate_variants}
	tabix ${candidate_variants}.gz
	"""
}

process Filter_Variants {
	
	cpus = 1
	memory = { 1 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	publishDir "results", mode: "copy"

	input:
		tuple val(meta), path(vcf), path(vcf_index)
		
	output:
		tuple val(meta), path("*_somatic_variants_filtered_gene_list.vcf.gz"), path("*_somatic_variants_filtered_gene_list.vcf.gz.tbi")
	
	script:
	"""
	TUMOR_INDEX=`bcftools query -l ${vcf} | grep -n ${meta.tumor} | cut -d: -f1 | awk '{print \$1-1}'`
	NORMAL_INDEX=`bcftools query -l ${vcf} | grep -n ${meta.normal} | cut -d: -f1 | awk '{print \$1-1}'`

	bcftools view --regions-file ${params.gene_bed} ${vcf} |

	#tee >(grep -v ^# | wc -l | awk '{print "Variants in gene list regions: " \$1}' >&2) |

	bcftools view -e "AC_gnomad_total_4.0>1000" |

	#tee >(grep -v ^# | wc -l | awk '{print "..and gnomAD v4 AC < 100: " \$1}' >&2) |

	bcftools view -i "AD["\$TUMOR_INDEX":1]>=4" |

	#tee >(grep -v ^# | wc -l | awk '{print "..and ALT AD >=4: " \$1}' >&2) |

	bcftools view -i "AD["\$TUMOR_INDEX":0] + AD["\$TUMOR_INDEX":1]>=10" |

	#tee >(grep -v ^# | wc -l | awk '{print "..and total AD >=10: " \$1}' >&2) |

	bcftools view -i "ADF["\$TUMOR_INDEX":1]>=2 && ADR["\$TUMOR_INDEX":1]>=2" |

	#tee >(grep -v ^# | wc -l | awk '{print "..and >= 2 ALT reads in each direction: " \$1}' >&2) |

	bcftools view -i "VDB>=0.05" |

	#tee >(grep -v ^# | wc -l | awk '{print "..and Variant Distance Bias (VDB) >=0.05: " \$1}' >&2) |

	#grep -E "^#|MODERATE|HIGH" |

	#tee >(grep -v ^# | wc -l | awk '{print "..and MODERATE or HIGH impact: " \$1}' >&2) |

	#bcftools view -i "INFO/CADD_Score > 10"|

	#tee >(grep -v ^# | wc -l | awk '{print "..and CADD_Score > 20: " \$1}' >&2) |

	bcftools view --output-type z --output-file ${meta.id}_somatic_variants_filtered_gene_list.vcf.gz

	tabix ${meta.id}_somatic_variants_filtered_gene_list.vcf.gz
	"""
}

process Create_Report {

    module 'curl/8.6.0'
    module 'R/4.5.1'
    module 'gcc/14.2'

	cpus = 1
	memory = { 8 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	publishDir "results", mode: "copy"

	input:
		tuple val(meta), path(vcf), path(vcf_index)
		
	output:
		tuple val(meta), path("*_somatic_variants_filtered_summary.xlsx")
	
    script:
    """
    Rscript /vast/projects/reidj-project/filter_somatic_variants.R ${vcf} ${meta.id}_somatic_variants_filtered_summary.xlsx
    
    """
}

process Create_BAMlet {

    module 'samtools'

    cpus = 1
    memory = { 8 * task.attempt + ' GB' }
    time = { 1 * task.attempt + ' h' }

    publishDir "results", mode: "copy"

    input:
    tuple val(meta), path(variant_summary)

    output:
    tuple val(meta), path("*.cram"), path("*.cram.crai")

    script:
    // Determine BAM or CRAM file path
    def bam_path = "${meta.bam_path}/${meta.tumor}/${meta.tumor}.recal.bam"
    def cram_path = "${meta.bam_path}/${meta.tumor}/${meta.tumor}.recal.cram"
    def BAM_FILE = ""

    if (file(bam_path).exists()) {
        BAM_FILE = bam_path
    } else if (file(cram_path).exists()) {
        BAM_FILE = cram_path
    } else {
        error "No .bam or .cram file found for ${meta.tumor}"
    }

    // Output files
    def out_cram = "${meta.tumor}.subset.cram"
    def out_crai = "${meta.tumor}.subset.cram.crai"

    """
    samtools view -b -L ${params.gene_bed} -o subset.bam ${BAM_FILE}

    samtools sort -o ${out_cram} -T temp_subset -O cram -T ${meta.tumor}_tmp.bam --reference ${params.ref_fasta} subset.bam

    samtools index ${out_cram}
    """
}


'''
bcftools +fill-tags --threads 4 somatic_variants_annotated_cadd.vcf.gz -O z -o ${meta.tumor}_somatic_variants_annotated.vcf.gz -- -t AF,AN,AC	
tabix -f ${meta.tumor}_somatic_variants_annotated.vcf.gz				
'''

workflow {

    ch_samples = Channel
        .fromPath(params.sample_info)
        .splitCsv(header: true, sep: "\t")
        .map { row ->
            def meta = [
                id:       row.ID,
                normal:    row.NORMAL,
                tumor:    row.TUMOR,
                bam_path: row.BAM_PATH,
                vcf_path: row.VCF_PATH
            ]
            meta
        }

    ch_processed = ProcessSample(ch_samples)

    ch_split_vcf = Split_Vcf(ch_processed.MergedVCF)
        .transpose()
        .combine(ch_processed.BamFILE, by: 0)

    ch_annotated = ch_split_vcf
        | mpileup_check
        | ensembl_vep
        | gnomad
        | CADD_Run_Container
        | Process_CADD

    ch_grouped = ClinVar(ch_annotated)
        .groupTuple(by: 0)

    ch_joined = Join_VCF(ch_grouped)

    ch_joined.JoinedVCF
        .join(ch_processed.IndividualTools, by: 0)
        | Check_Tools
        | Compress_Index
        | Filter_Variants
        | Create_Report
        | Create_BAMlet
}


'''

workflow {

    // Create a channel from the sample info file with proper meta map structure
    ch_samples = Channel
        .fromPath(params.sample_info)
        .splitCsv(header: true, sep: "\t")
        .map { row ->
            def meta = [
                id: row.ID,
                normal: row.NORMAL,
                tumor: row.TUMOR,
                bam_path: row.BAM_PATH,
                vcf_path: row.VCF_PATH
            ]
            return meta
        }

    // Process the sample channel - ensure meta is carried through
    ch_processed = ProcessSample(ch_samples)

    // Split VCF while maintaining sample identity
    ch_split_vcf = Split_Vcf(ch_processed.MergedVCF)
        .transpose()

    // Combine split VCF with corresponding BAM file using the meta map
    ch_combined = ch_split_vcf
        .combine(ch_processed.BamFILE, by: 0)

    // Annotation pipeline - split into three processes
    ch_mpileup    = mpileup_check(ch_combined)
    ch_vep        = ensembl_vep(ch_mpileup)
    ch_annotated  = gnomad(ch_vep)


    ch_cadd = CADD_Run_Container(ch_annotated)

    ch_cadd_processed = Process_CADD(ch_cadd)
   
    ch_clinvar = ClinVar(ch_cadd_processed)

    // Group by sample ID before joining VCFs
    ch_grouped = ch_clinvar
        .groupTuple(by: 0)

    ch_joined = Join_VCF(ch_grouped)

    // Merge channels
    ch_merged_for_check = ch_joined.JoinedVCF
        .join(ch_processed.IndividualTools, by: 0)

    ch_check_tools = Check_Tools(ch_merged_for_check)

    ch_compessed = Compress_Index(ch_check_tools)

    ch_filtered = Filter_Variants(ch_compressed)

    ch_report = Create_Report(ch_filtered)

    ch_bamlet = Create_BAMlet(ch_report)

   
} 

'''



process Annotate {

	module 'ensembl-vep/112'

	cpus = 1
	memory = { 5 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	input: tuple val(meta), path(split_vcf), path(BamFILE)

	output: tuple val(meta), path ("somatic_variants_annotated.vcf")

    script:
		"""
		#Check all positions
		bcftools mpileup  -x --indels-2.0 --max-depth 1000000 --max-idepth 100000 --gap-frac 0.0001 --min-BQ 20 --min-MQ 1 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR --fasta-ref ${params.ref_fasta} --targets-file ${split_vcf} ${BamFILE} |
		bcftools norm  --multiallelics -any -f ${params.ref_fasta} |
		bcftools call -mA |
		bcftools view -i "FORMAT/ADF[0:1]>=1 && FORMAT/ADR[0:1]>=1 && FORMAT/AD[0:1]>=3" -Oz -o mpileup_version.vcf.gz

		#Run ensembl-vep
		echo "Running ensembl-vep..."
		vep --cache --dir /stornext/Bioinf/data/lab_bahlo/ref_db/vep-cache/ --cache_version 104 --assembly GRCh38 \
				-i mpileup_version.vcf.gz -o merged_vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs \
				--fasta ${params.ref_fasta} --offline --sift b --polyphen b --ccds --hgvs --hgvsg --symbol \
				--numbers --protein --af --af_1kg --max_af --variant_class --pick_allele_gene --force_overwrite
		echo "ensembl-vep complete."

		#gnomAD annotations
		vcfanno -p 4 /stornext/Bioinf/data/lab_bahlo/users/reid.j/vcfanno/gnomad_v4.0.0.toml merged_vep.vcf > somatic_variants_vep_gnomad.vcf
		vcfanno -p 4 /stornext/Bioinf/data/lab_bahlo/users/reid.j/vcfanno/gnomad_v4_postprocess.toml somatic_variants_vep_gnomad.vcf > somatic_variants_annotated.vcf
		"""
}


process Check_Tools_OG {
	
	cpus = 1
	memory = { 1 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	publishDir "results", mode: "copy"

	input:
		tuple val(meta), path(vcf), path(vcf_index), path(mutect), path(freebayes), path(strelka)
		
	output:
		tuple val(meta), path("*_somatic_variants_merged.vcf.gz"), path("*_somatic_variants_merged.vcf.gz.tbi")
	
	shell:
	'''
	# Step 1: Generate the TOML configuration file for VCFAnno
	cat <<EOL > "annotations_config.toml"
	[[annotation]]
	file="mutect_v4.vcf.gz"
	fields=["Mutect"]
	ops=["self"] 

	[[annotation]]
	file="strelka_v5.vcf.gz"
	fields=["Strelka"]
	ops=["self"]

	[[annotation]]
	file="freebayes_v4.vcf.gz"
	fields=["Freebayes"]
	ops=["self"]

	[[postannotation]]
	fields=["Mutect", "Strelka", "Freebayes"]
	op="sum"
	name="Tool_Count"
	type="Float"
	EOL

	# Step 2: Apply annotations using VCFAnno
	vcfanno -p 4 annotations_config.toml !{vcf} > !{meta.id}_somatic_variants_merged.vcf
	bgzip -f !{meta.id}_somatic_variants_merged.vcf
	tabix -f !{meta.id}_somatic_variants_merged.vcf.gz
	'''
}
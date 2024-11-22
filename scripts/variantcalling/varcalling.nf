#!/usr/bin/env nextflow

nextflow.enable.dsl=1

ref_file_37 = params.ref_file_37
ref_file_38 = params.ref_file_38
intervals = "/path/to/intervals"
gvcfs_dir = "/path/to/gvcfs"
gvcfs_list = "/path/to/gvcflist"
gvcfsdb_new = "/path/to/gvcfsdb"
in_file = Channel.fromFilePairs(params.in_fastq, type: 'file')
in_file2 = Channel.fromFilePairs(params.in_fastq2, type: 'file')
in_bam = Channel.fromFilePairs(params.in_bam, type: 'file')
// in_file = Channel.fromPath (params.in_fastq)

// in_file.view()

process align {
//      maxForks 10

//      publishDir params.align_dir, mode: 'copy', overwrite: 'true'

        errorStrategy 'ignore'
        tag "${name}"

        input:
        set val(name), file(fastq) from in_file
        output:
        set val(name), file("${name}_sorted.bam"), file("${name}_sorted.bam.bai") into align_ch

        script:
        """
        bwa mem \
        -R "@RG\\tID:HHTN2BBXX.0\\tLB:LIBA\\tSM:${name}\\tPL:Illumina" \
        -t 31 -K 100000 -Y $ref_file_38 \
        ${fastq[0]} ${fastq[1]} | samtools sort -@ 31 -m 500M - > ${name}_sorted.bam
        samtools index ${name}_sorted.bam
        """
}

process align2 {


      input:
         set val(name), file(fastq) from in_file2
      output:
         set val(name), file("${name}_sorted2.bam"), file("${name}_sorted2.bam.bai") into align_ch2

         script:
         """
         bwa mem \
         -R "@RG\\tID:HHTN2BBXX.0\\tLB:LIBA\\tSM:${name}\\tPL:Illumina" \
         -t 31 -K 100000 -Y $ref_file_38 \
         ${fastq[0]} ${fastq[1]} | samtools sort -@ 31 -m 500M - > ${name}_sorted2.bam
         samtools index ${name}_sorted2.bam
         """

}

//align_ch.join(align_ch2).set {mergealign_ch}

//mergealign_ch.view()
process mergebams {

      publishDir params.align_dir, mode: 'copy', overwrite: 'true'


      input:
        set val(name), file(bam1), file(index1), file(bam2), file(index2) from mergealign_ch
      output:
        set val(name), file("${name}_merged.bam"), file("${name}_merged.bam.bai") into mergedbam_ch

      script
      """
      samtools merge ${name}_merged.bam ${bam1} ${bam2}
      samtools index ${name}_merged.bam
      """

}
process VarCall {
// // MaxForks 10

      publishDir params.varcall_dir, mode: 'copy', overwrite: 'true'

      errorStrategy 'ignore'
      tag "${name}"

      input:
        set val(name), file(bam), file(bai) from align_ch
      output:
        set val(name), file("*") into gvcf_ch

      script:
      """
      gatk HaplotypeCaller -R $ref_file_38 -L $intervals -I $bam -ERC GVCF -stand-call-conf 10 -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth     -A RMSMappingQuality -A ReadPosRankSumTest --allow-old-rms-mapping-quality-annotation-data -O ${name}.g.vcf
      """

}
process createDB {

//      publishDir "$gvcfs_dir", mode: 'copy', overwrite: 'true'

      input:
        path gvcfs_dir
      output:
        path ("gvcfsdb") into gvcfs_ch

      script:
      """
      gatk GenomicsDBImport --variant $gvcfs_list --genomicsdb-workspace-path gvcfsdb --intervals $intervals
      """
}

process genotypevcfs {

      input:
        path gvcfsdb_new
      output:
        set file ("pgxjointcall.vcf"), file ("pgxjointcall.vcf.idx") into jointcall_ch

      script:
      """
      gatk GenotypeGVCFs -R $ref_file_38 -V gendb://$gvcfsdb_new -G StandardAnnotation -O pgxjointcall.vcf
      """

}
process run_select_snps {

      input:
        set file (vcf), file(vcf_index) from jointcall_ch1
      output:
        file ("*") into snps_ch

      script:
      """
      gatk SelectVariants -R $ref_file_38 -select-type SNP -V $vcf -O pgxjointcall_SNPS.vcf
      """

}

process run_select_indels {

      input:
        set file (vcf), file(vcf_index) from jointcall_ch2
      output:
        file ("*") into indels_ch

      script:
      """
      gatk SelectVariants -R $ref_file_38 -select-type INDEL -V $vcf -O pgxjointcall_INDEL.vcf
      """
}
process filterSNPs {

      input:
        set file(vcf), file(vcf_index) from snps_ch
      output:
        file ("*") into snpsfiltered_ch

      script:
      """
      gatk VariantFiltration \
      -R ${ref_file_38} \
      --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
      --filter-expression "FS > 60.0" --filter-name "FS_gt_60" \
      --filter-expression "MQ < 40.0" --filter-name "MQ_lt_40" \
      --filter-expression "MQRankSum < -12.5" --filter-name "MQRS_lt_n12.5" \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS_lt_n8" \
      --filter-expression "SOR > 3.0" --filter-name "SOR_gt_3" \
      -V ${vcf} \
      -O pgx.snps.filtered.vcf
      """

}
process filterINDELS {

      input:
        set file(vcf), file(vcf_index) from indels_ch
      output:
        file ("*") into indelsfiltered_ch

      script:
      """
      gatk VariantFiltration \
      -R ${ref_file_38} \
      --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
      --filter-expression "FS > 200.0" --filter-name "FS_gt_200" \
      --filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRS_lt_n20" \
      --filter-expression "SOR > 10.0" --filter-name "SOR_gt_10" \
      -V ${vcf} \
      -O pgx.INDELs.filtered.vcf
      """

}
process merge_snps_indels {

      input:
        set file(vcf_s), file(vcf_s_index) from snpsfiltered_ch
        set file(vcf_i), file(vcf_i_index) from indelsfiltered_ch
      output:
        file("*") into filtered_ch

      script:
      """
      gatk SortVcf \
      -I ${vcf_s} \
      -I ${vcf_i} \
      -O pgx.filtered.vcf
      """
}

process select_pass {

//      publishDir "$PWD/filteredvcf", mode: 'copy', overwrite: false
      input:
        set file(vcf), file(vcf_index) from filtered_ch
      output:
        file("*") into pass_ch

      script:
      """
      gatk SelectVariants \
      -R ${ref_file_38} \
      --exclude-filtered \
      -V ${vcf} \
      -O pgx.filtered.pass.vcf
      """

}

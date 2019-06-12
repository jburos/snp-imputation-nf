params.bedFile = 'TODO'
params.famFile = 'TODO'
params.bimFile = 'TODO'
params.outdir = "imputation-results" 
params.chromosomeSizesFile = 'b37.chrom.sizes'
params.referenceBase = '1000GP_Phase3'
params.referenceHapsFilePattern = "1000GP_Phase3_chr%s.hap.gz"
params.referenceLegendFilePattern = "1000GP_Phase3_chr%s.legend.gz"
params.referenceGeneticMapPattern = "genetic_map_chr%s_combined_b37.txt"
params.referenceSample = "1000GP_Phase3.sample"

bedFileChan = Channel.fromPath(params.bedFile)
famFileChan = Channel.fromPath(params.famFile)
bimFileChan = Channel.fromPath(params.bimFile)

chromosomesList = 1..22

refChannel = Channel.fromPath(params.referenceBase)

println """\
         I M P U T E 2 - N F   P I P E L I N E    
         ===================================
         bedfile      : ${params.bedFile}
         geneticMap   : ${params.referenceBase}
         outdir       : ${params.outdir}
         """
         .stripIndent()

process rmDups {

  container 'jackinovik/docker-impute2'

  input:
  file bedFile from bedFileChan
  file famFile from famFileChan
  file bimFile from bimFileChan 

  output:
  file "cleaned_nodups.fam" into cleaned_fam
  file "cleaned_nodups.bed" into cleaned_bed
  file "cleaned_nodups.bim" into cleaned_bim
  
  """
  # remove duplicate SNPs based on position
  plink2 --bim ${bimFile} --bed ${bedFile} --fam ${famFile} --set-all-var-ids @:#[b37] --make-bed --out step0_set_var_ids
  plink2 --bfile step0_set_var_ids --rm-dup exclude-mismatch list --make-bed --out step1_excl_mismatch
  # prepare text file to use when restoring original snp names
  sort -k4 ${bimFile} | cut -f 4,2 > orig_names.txt 
  sort -k4 step1_excl_mismatch.bim | cut -f 4,2 > new_names.txt
  # keep only rsids & only one rsid per position
  join -1 2 -2 2 new_names.txt orig_names.txt -o 1.1,2.1 | grep rs | sort -k1,1 -u > restore_names.txt
  plink2 --bfile step1_excl_mismatch --update-name restore_names.txt --make-bed --out cleaned_nodups
  """
}

process splitChrs {

  container 'jackinovik/docker-impute2'

  input:
  file fam from cleaned_fam
  file bed from cleaned_bed
  file bim from cleaned_bim
  each chromosome from chromosomesList 

  output:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") into perChrom

  """
  plink2 --fam $fam --bed $bed --bim $bim --chr $chromosome --make-bed --out chr${chromosome}
  """
}

process harmonizeGenotypes {

  container 'jackinovik/docker-impute2'
  
  input:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") from perChrom 
  
  output:
  set val(chromosome), file("chr${chromosome}_aligned.bed"), file("chr${chromosome}_aligned.fam"), file("chr${chromosome}_aligned.bim") into perChromAligned

  """
  # download VCF file from 1000GP phase3
  wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
  # run genotypeHarmonizer
  java -Xmx40g -jar /install/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar --inputType PLINK_BED --input chr${chromosome} --update-id --outputType PLINK_BED --output chr${chromosome}_aligned --refType VCF --ref ALL.chr${chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes --ambiguousSnpFilter --check-ld --update-reference-allele
  # remove vcf to conserve disk space
  rm *.vcf.gz*
  """
}

process shapeitCheck {
  //validExitStatus 0,1,2
  //errorStrategy 'ignore'

  container 'jackinovik/docker-impute2'

  input:
  set val(chromosome), file("chr${chromosome}_aligned.bed"), file("chr${chromosome}_aligned.bim"), file("chr${chromosome}_aligned.fam") from perChromAligned
  file db_path from refChannel

  output:
  set val(chromosome), file("chr${chromosome}.alignments.log"), file("chr${chromosome}.alignments.snp.strand.exclude"), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") into shapitCheckChan

  script:
  hapFile = file( db_path.name + "/" + sprintf(params.referenceHapsFilePattern, chromosome) )
  legendFile = file( db_path.name + "/" + sprintf(params.referenceLegendFilePattern, chromosome) )
  sampleFile = file( db_path.name + "/" + params.referenceSample )

  """
  shapeit -check --input-bed chr${chromosome}.bed chr${chromosome}.bim chr${chromosome}.fam --input-ref $hapFile $legendFile $sampleFile --output-log chr${chromosome}.alignments
  """

}




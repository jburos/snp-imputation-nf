params.bedFile = 'TODO'
params.famFile = 'TODO'
params.bimFile = 'TODO'
params.outdir = "imputation-results" 
params.chromosomeSizesFile = 'b37.chrom.sizes'
params.referenceBase = 'ALL.integrated_phase1_SHAPEIT_16-06-14.nomono'
params.referenceHapsFilePattern = "ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz"
params.referenceLegendFilePattern = "ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz"
params.referenceGeneticMapPattern = "genetic_map_chr%s_combined_b37.txt"
params.referenceSample = "ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample"

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

process splitChrs {

  container 'jackinovik/docker-impute2'

  input:
  file bedFile from bedFileChan
  file famFile from famFileChan
  file bimFile from bimFileChan 
  each chromosome from chromosomesList 

  output:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") into perChromChan

  """
  plink2 --bim ${bimFile} --bed ${bedFile} --fam ${famFile} --chr $chromosome --make-bed --out chr${chromosome}_pre1
  plink2 -bfile chr${chromosome}_pre1 --rm-dup exclude-mismatch list --make-bed --out chr${chromosome}
  """

}

process shapeitCheck {
  //validExitStatus 0,1,2
  //errorStrategy 'ignore'

  container 'insilicodb/docker-impute2'

  input:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") from perChromChan
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




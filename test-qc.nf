params.bedFile = 'TODO'
params.famFile = 'TODO'
params.bimFile = 'TODO'
params.geneticMapTgz = 'TODO'
params.outDir = "imputation-results" 
params.chromosomeSizesFile = 'b37.chrom.sizes'
params.referenceHapsFilePattern = "ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz"
params.referenceLegendFilePattern = "ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz"
params.referenceGeneticMapPattern = "genetic_map_chr%s_combined_b37.txt"
params.referenceSample = "ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample"

bedFileChan = Channel.fromPath(params.bedFile)
famFileChan = Channel.fromPath(params.famFile)
bimFileChan = Channel.fromPath(params.bimFile)

chromosomesList = 1..22

genetic_map_tgz = file(params.geneticMapTgz)

println """\
         I M P U T E 2 - N F   P I P E L I N E    
         ===================================
         bedfile      : ${params.bedFile}
         geneticMap   : ${params.geneticMapTgz}
         outdir       : ${params.outDir}
         """
         .stripIndent()

process getGeneticMap {

   container 'jackinovik/docker-impute2'
   
   maxForks 1
   
   input:
   file genetic_map_tgz
   
   output:
   file '*.nomono' into genetic_map_dir
   
   """
   tar xfz ${genetic_map_tgz}
   """
}

process qcFilters {

  container 'jackinovik/docker-impute2'

  maxForks 6
  
  input:
  file bedFile from bedFileChan
  file famFile from famFileChan
  file bimFile from bimFileChan
  
  output:
  set file("qcresult.bed"), file("qcresult.fam"), file("qcresult.bim") into qcResultChan

  """
  plink2 --bim ${bimFile} --bed ${bedFile} --fam ${famFile} --geno 0.2 --make-bed --out qc01_geno
  plink2 --bfile qc01_geno --mind 0.2 --make-bed --out qc02_mind
  ## Check sex
  plink2 --bfile qc02_mind --split-x --make-bed --out qc03_splitx
  plink2 --bfile qc03_splitx --check-sex --out qc_check_sex
  ## TODO report on mismatching sex, possibly remove
  plink2 --bfile qc03_splitx --maf 0.05 --make-bed --out qc04_maf 
  ## TODO filter based on hwe among controls
  # plink2 --bfile qc04_maf --keep-if  == control
  """
}

process splitChrs {

  container 'jackinovik/docker-impute2'

  maxForks 6

  input:
  file bedFile from bedFileChan
  file famFile from famFileChan
  file bimFile from bimFileChan 
  each chromosome from chromosomesList 

  output:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") into plinkOutChan

  """
  plink2 --noweb --bim ${bimFile} --bed ${bedFile} --fam ${famFile} --chr $chromosome --make-bed --out chr${chromosome}
  plink2 -bfile chr${chromosome} --rm-dup exclude-mismatch list
  """

}

//plinkOutChan.subscribe{  println "${it.name}"}



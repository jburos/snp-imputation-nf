params.bedFile = 'TODO'
// Assume bim/fam files are named along with bedfile
params.geneticMapTgz = 'TODO'
params.outDir = "imputation-results" 
params.chromosomeSizesFile = 'b37.chrom.sizes'
params.referenceHapsFilePattern = "ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz"
params.referenceLegendFilePattern = "ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz"
params.referenceGeneticMapPattern = "genetic_map_chr%s_combined_b37.txt"
params.referenceSample = "ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample"

bedFileChan = Channel.fromPath(params.bedFile)
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

   container 'insilicodb/docker-impute2'
   
   maxForks 1
   
   input:
   file genetic_map_tgz
   
   output:
   file 'genetic_map' into genetic_map_dir
   
   """
   mkdir genetic_map
   tar xfz ${genetic_map_tgz} -C genetic_map
   """

}

process plink {

  container 'insilicodb/docker-impute2'

  maxForks 6

  input:
  file bedFile from bedFileChan
  each chromosome from chromosomesList 

  output:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") into plinkOutChan

  """
  plink --noweb --bfile ${bedFile.baseName} --chr $chromosome --make-bed --out chr${chromosome}
  plink -bfile chr${chromosome} --list-duplicate-vars ids-only suppress-first
  [[ -e "plink.dupvar" ]] && plink --bfile chr${chromosome} --exclude plink.dupvar --make-bed --out chr${chromosome}
  """

}

//plinkOutChan.subscribe{  println "${it.name}"}



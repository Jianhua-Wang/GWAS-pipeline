#!/usr/bin/env nextflow

/*
 * Author    :   Jianhua Wang
 * Date      :   04-09-2019
 *
 *
 * This is a pipeline for imputation
 * 
 * 
 * INPUT     :   plink format (.bed,.bim,.fam) after QC
 * OUTPUT    :   imputation result for per chromosome (.impute2)
 * REF PANEL :   1. (chr*.hap, chr*.leg) from 1000G: https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html
 * 
 *
 */

/*
 * 
 * 1. Split plink format into .gen by chromosomes
 * 2. Imputation with every 5,000,000 BP
 * 3. Concat imputation results
 * 4. Convert gen file to bgen file
 *
 */

/************************************************
*****        Input Preparation        ***********
*************************************************/
log.info '''
==================================================================
     \033[1;33m/\\\033[0m
    \033[1;33m/__\\\033[0m\033[1;31m\\\033[0m         This is a pipeline for imputation
   \033[1;33m/\033[0m  \033[1;31m---\\\033[0m        Author: Jianhua Wang
  \033[1;33m/\\\033[0m      \033[1;31m\\\033[0m       Date:   04-09-2019
 \033[1;33m/\033[0m\033[1;32m/\\\033[0m\033[1;33m\\\033[0m     \033[1;31m/\\\033[0m
 \033[1;32m/  \\   /\033[0m\033[1;31m/__\\\033[0m
\033[1;32m`----`-----\033[0m
==================================================================
    '''.stripIndent()
chromosomes_List = params.chromosomes_List

// Checks if the file exists
checker = { fn ->
    if (fn.exists())
        return fn;
    else
        error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}

// Input to imputation
plink_pat = "${params.input_dir}/${params.input_pat}"
plink_prefix = params.input_pat
plink_ch = Channel.fromFilePairs("${plink_pat}.{bed,bim,fam}", size:3, flat : true){ file -> file.baseName }\
                  .ifEmpty { error "No matching plink files" }\
                  .map { a -> [checker(a[1]), checker(a[2]), checker(a[3])] }

// Genetic maps
map_dir               = params.map_dir
map_pattern           = params.map_pattern

// 1000G Reference Panel
ref_panel_dir         = params.ref_panel_dir
ref_hap_pattern       = params.ref_hap_pattern
ref_leg_pattern       = params.ref_leg_pattern
ref_sample            = params.ref_sample

/*
 * start pipe
 * 1. Split plink format into .gen by chromosomes
 */

process plink {
  //maxForks 4

  input:
  each chromosome from chromosomes_List
  set file(bed), file(bim), file(fam) from plink_ch

  output:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") into plinkOutChan
  
  script:
  """
  plink --bfile $plink_prefix \
  --chr $chromosome \
  --make-bed \
  --out chr$chromosome
  """
}

process shapeitCheck {
  validExitStatus 0,1,2
  errorStrategy 'ignore'

  input:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") from plinkOutChan

  output:
  set val(chromosome), file("chr${chromosome}.alignments.log"), file("chr${chromosome}.alignments.snp.strand.exclude"), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") into shapitCheckChan

  script:
  hapFile = file( ref_panel_dir + sprintf(ref_hap_pattern, chromosome) )
  legendFile = file( ref_panel_dir + sprintf(ref_leg_pattern, chromosome) )
  sampleFile = file( params.ref_panel_dir + params.ref_sample )
  mapFile    = file( params.map_dir + sprintf(params.map_pattern, chromosome) )

  """
  shapeit \
  -check \
  -M ${mapFile} \
  --input-bed chr${chromosome}.bed chr${chromosome}.bim chr${chromosome}.fam \
  --input-ref $hapFile $legendFile $sampleFile \
  --output-log chr${chromosome}.alignments
  """

}

process shapeit {

  maxForks params.maxForks_shapeit

  input:
  set val(chromosome), file("chr${chromosome}.alignments.log"), file("chr${chromosome}.alignments.snp.strand.exclude"), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") from shapitCheckChan

  output:
  set val(chromosome), file("chr${chromosome}.phased.haps"), file("chr${chromosome}.phased.sample") into shapeitChan

  script:
  hapFile = file( ref_panel_dir + sprintf(ref_hap_pattern, chromosome) )
  legendFile = file( ref_panel_dir + sprintf(ref_leg_pattern, chromosome) )
  sampleFile = file( params.ref_panel_dir + params.ref_sample )
  mapFile    = file( params.map_dir + sprintf(params.map_pattern, chromosome) )
  excludeFile = "chr${chromosome}.alignments.snp.strand.exclude"

  """
  shapeit \
  --input-bed chr${chromosome}.bed chr${chromosome}.bim chr${chromosome}.fam \
  --input-ref $hapFile $legendFile $sampleFile \
  --exclude-snp $excludeFile \
  --input-map $mapFile \
  -O chr${chromosome}.phased \
  --thread 20 \
  --force
  """

}

def getChromosomeSize( chromosomeSizesFile, chromosome ) {

    def result = 0

    chromosomeSizesFile.splitEachLine("\t") {fields ->

        def genomeId
        def path
	if ( fields[0].trim() == "${chromosome}" ) {
          //println "in if"
          result = fields[1].trim().toInteger()
          return 
        }

    }

    result
}

def getChromosomeChunkPairs ( chromosomeSize, chunkSize=5000000 ) {

  def result = []
  def numberOfChunks = chromosomeSize / chunkSize
  def remainder = chromosomeSize % chunkSize
 
  1.upto(numberOfChunks) {
    endPosition = it * chunkSize
    startPosition = (endPosition - chunkSize) + 1
    result = result + [[startPosition, endPosition ]]
  }

  if ( remainder > 0 ) {
    result = result + [[endPosition + 1 , endPosition + remainder ]]
  }
  
  result
}

imputeChromChunckChannel = shapeitChan.flatMap { chromosome, gensFile, sampleFile ->
   def results = []
   
   def chunks = getChromosomeChunkPairs(getChromosomeSize(file(params.chromosomeSizesFile), chromosome))
   
   chunks.each { chunkStart, chunkEnd -> 
     results.push( [ chromosome, gensFile, sampleFile, chunkStart, chunkEnd] )
   }
   
   return results 
}

// 2. imputation
process impute2 {
  
  maxForks params.maxForks_impute2
  validExitStatus 0,1,2
  errorStrategy 'ignore'
  
  input:
  set val(chromosome), file("chr${chromosome}.phased.haps"), file("chr${chromosome}.phased.sample"), val(chunkStart), val(chunkEnd) from imputeChromChunckChannel
  
  output:
  set val(chromosome), file("chr${chromosome}-${chunkStart}-${chunkEnd}.imputed") into impute2Chan
  set val(chromosome), file("chr${chromosome}-${chunkStart}-${chunkEnd}.imputed_info") into infoChan

  script:
  mapFile    = file( map_dir + sprintf(map_pattern, chromosome) )
  hapFile    = file( ref_panel_dir + sprintf(ref_hap_pattern, chromosome) )
  legendFile = file( ref_panel_dir + sprintf(ref_leg_pattern, chromosome) )
  """
  impute2 \
  -use_prephased_g \
  -known_haps_g chr${chromosome}.phased.haps \
  -m ${mapFile} \
  -h ${hapFile} \
  -l ${legendFile} \
  -int $chunkStart $chunkEnd \
  -Ne 20000 \
  -o chr${chromosome}-${chunkStart}-${chunkEnd}.imputed
  if [ ! -f "chr${chromosome}-${chunkStart}-${chunkEnd}.imputed" ]; then
    touch "chr${chromosome}-${chunkStart}-${chunkEnd}.imputed";
    touch "chr${chromosome}-${chunkStart}-${chunkEnd}.imputed_info";
  fi

  """
}

impute2List = impute2Chan.toSortedList() 
impute2List = impute2List.val

impute2Map = [:]

impute2List.each { chrom, file ->

  if ( !impute2Map.containsKey(chrom) ) {
   impute2Map.put(chrom, [])
  }
  impute2Map.get(chrom).add(file)

}

impute2MapChannel = Channel.create()

impute2Map.each { chrom, fileList ->
  impute2MapChannel.bind([chrom, fileList])
}

impute2MapChannel.close()

// 3. concat impute2 results for per chromosome
process impute2Concat {
  
 
  input:
  set val(chromosome), file(imputedFiles) from impute2MapChannel

  output:
  set val(chromosome), file("chr${chromosome}_1KG.gen") into impute2ConcatChan

  """
  cat $imputedFiles > chr${chromosome}_1KG.gen
  """
}

// 4. Convert gen file to bgen file

impute2ConcatChan
 .toSortedList( { a, b -> a[0] <=> b[0] } ) 
 .flatten()
 .buffer( size:1, skip:1 )
 .flatMap{it.get(0)}
 .toList()
 .set{res}

process gen2bgen {
  
  publishDir params.output_dir, overwrite:true, mode:'link'
 
  input:
  file(genFiles) from res

  output:
  file("${params.input_pat}.bgen") into bgenChan

  """
  qctool -g chr#_1KG.gen -og ${params.input_pat}.bgen
  """
}
#!/usr/bin/env nextflow

/*
 * Author    :   Jianhua Wang
 * Date      :   10-27-2018
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
 * 1. Convert Unphased VCF to .gen.gz by chromosomes
 * 2. Split plink format into .gen by chromosomes
 * 3. Imputation with every 5,000,000 BP
 * 4. Concat imputation results
 *
 */

/************************************************
*****        Input Preparation        ***********
*************************************************/

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

/*
 * start pipe
 * 1. Convert Unphased VCF to .gen.gz by chromosomes
 * 2. Split plink format into .gen by chromosomes
 */

process plink_to_gen {
  //maxForks 4
  input:
  each chromosome from chromosomes_List
  set file(bed), file(bim), file(fam) from plink_ch

  output:
  set val(chromosome), file("${plink_prefix}_chr${chromosome}.gen"), file("${plink_prefix}_chr${chromosome}.sample") into plink_gen_ch
  
  script:
  """
  plink --bfile $plink_prefix \
  --chr $chromosome \
  --recode oxford \
  --out ${plink_prefix}_chr$chromosome
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

imputeChromChunckChannel = plink_gen_ch.flatMap { chromosome, gensFile, sampleFile ->
   def results = []
   
   def chunks = getChromosomeChunkPairs(getChromosomeSize(file(params.chromosomeSizesFile), chromosome))
   
   chunks.each { chunkStart, chunkEnd -> 
     results.push( [ chromosome, gensFile, sampleFile, chunkStart, chunkEnd] )
   }
   
   return results 
}

// 3. imputation
process impute2 {
  maxForks params.maxForks
  validExitStatus 0,1,2
  errorStrategy 'ignore'

  input:
  set val(chromosome), file("${plink_prefix}_chr${chromosome}.gen"), file("${plink_prefix}_chr${chromosome}.sample"), val(chunkStart), val(chunkEnd) from imputeChromChunckChannel
  
  output:
  set val(chromosome), file("chr${chromosome}-${chunkStart}-${chunkEnd}.imputed"), file("chr${chromosome}-${chunkStart}-${chunkEnd}.imputed_info") into impute2Chan

  script:
  mapFile    = file( map_dir + sprintf(map_pattern, chromosome) )
  hapFile    = file( ref_panel_dir + sprintf(ref_hap_pattern, chromosome) )
  legendFile = file( ref_panel_dir + sprintf(ref_leg_pattern, chromosome) )
  """
  impute2 \
  -m ${mapFile} \
  -h ${hapFile} \
  -l ${legendFile} \
  -g ${plink_prefix}_chr${chromosome}.gen \
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

// 4. concat impute2 results for per chromosome
process impute2Concat {
  
  publishDir params.output_dir
 
  input:
  set val(chromosome), file(imputedFiles), file(infoFiles) from impute2MapChannel

  output:
  set val(chromosome), file("chr${chromosome}_1KG.imputed") into impute2ConcatChan

  """
  cat $imputedFiles > chr${chromosome}_1KG.imputed
  cat $infoFiles > chr${chromosome}_1KG.imputed_info
  """
}
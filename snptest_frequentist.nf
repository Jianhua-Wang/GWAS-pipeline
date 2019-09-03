#!/usr/bin/env nextflow

/*
 * Author    :   Jianhua Wang
 * Date      :   10-27-2018
 *
 *
 * This is a pipeline for association tese using SNPTEST
 * 
 * 
 * INPUT     :   bgen format (.bgen,.sample) after imputation
 * OUTPUT    :   association results of per chromosome
 * 
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
bgen_pat = "${params.input_dir}/${params.input_pat}"
bgen_prefix = params.input_pat
bgen_ch = Channel.fromFilePairs("${bgen_pat}.{bgen,sample}", size:2, flat : true){ file -> file.baseName }\
                  .ifEmpty { error "No matching files" }\
                  .map { a -> [checker(a[1]), checker(a[2])] }

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

bgenChromChannel = bgen_ch.flatMap { bgenFile, sampleFile ->
    def results = []
    
    chromosomes_List.each { chromosome -> 
      results.push( [ chromosome, bgenFile, sampleFile] )
    }
    
    return results
}

bgenChromChunckChannel = bgenChromChannel.flatMap { chromosome, bgenFile, sampleFile ->
    def results = []
    
    def chunks = getChromosomeChunkPairs(getChromosomeSize(file(params.chromosomeSizesFile), chromosome))
    
    chunks.each { chunkStart, chunkEnd -> 
      results.push( [ chromosome, bgenFile, sampleFile, chunkStart, chunkEnd] )
    }
    
    return results
}

process snptest {
  
  maxForks params.maxForks_SNPTEST
  // publishDir params.output_dir, overwrite:true, mode:'link'
  
  input:
  set val(chromosome), file(bgenFile), file(sampleFile), val(chunkStart), val(chunkEnd) from bgenChromChunckChannel
  
  output:
  set val(chromosome), file("chr${chromosome}-${chunkStart}-${chunkEnd}.score") into snptestChan

  script:
  """
  snptest \
  -data ${bgenFile} ${sampleFile} \
  -range ${chromosome}:${chunkStart}-${chunkEnd} \
  -o chr${chromosome}-${chunkStart}-${chunkEnd}.score \
  -method score \
  -frequentist 1 \
  -pheno tiv \
  -cov_names age height sex
  sed -i '1,13d' chr${chromosome}-${chunkStart}-${chunkEnd}.score
  sed -i '\$d' chr${chromosome}-${chunkStart}-${chunkEnd}.score
  """
}

snptestList = snptestChan.toSortedList() 
snptestList = snptestList.val

snptestMap = [:]

snptestList.each { chrom, file ->

  if ( !snptestMap.containsKey(chrom) ) {
   snptestMap.put(chrom, [])
  }
  snptestMap.get(chrom).add(file)

}

snptestMapChannel = Channel.create()

snptestMap.each { chrom, fileList ->
  snptestMapChannel.bind([chrom, fileList])
}

snptestMapChannel.close()

process snptestConcatChunk {
  
  // publishDir params.output_dir, overwrite:true, mode:'link'
 
  input:
  set val(chromosome), file(snptestFiles) from snptestMapChannel

  output:
  set val(chromosome), file("chr${chromosome}.score") into snptestConcatChunkChan

  script:
  """
  cat $snptestFiles > chr${chromosome}.score
  """
}

process snptestConcatChro {
  
  publishDir params.output_dir, overwrite:true, mode:'link'
 
  input:
  file(chroFiles) from snptestConcatChunkChan.collect()

  output:
  file("${bgen_prefix}.score.txt.gz") into snptest_out_ch

  script:
  """
  echo "alternate_ids rsid chromosome position alleleA alleleB index average_maximum_posterior_call info cohort_1_AA cohort_1_AB cohort_1_BB cohort_1_NULL all_AA all_AB all_BB all_NULL all_total all_maf missing_data_proportion frequentist_add_pvalue frequentist_add_info frequentist_add_beta_1 frequentist_add_se_1 comment" > head
  cat head *.score > ${bgen_prefix}.score.txt
  gzip -f ${bgen_prefix}.score.txt
  """
}

process manhattan_qq {

  publishDir params.output_dir, overwrite:true, mode:'link'

  input:
  file(bgenie_out) from snptest_out_ch
  
  output:
  file("${base}_Manhattan_QQ.png")
  
  script:
  base = bgen_prefix
  template "manhattan_qq_bgenie.py"
}
#!/usr/bin/env nextflow

/*
 * Author    :   Jianhua Wang
 * Date      :   10-27-2018
 *
 *
 * This is a pipeline for association tese using BGENIE
 * 
 * 
 * INPUT     :   bgen format (.bgen,.bgen.bgi) after imputation
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
pheno = params.pheno
covar = params.covar
output_pat = params.output_pat
bgen_ch = Channel.fromFilePairs("${bgen_pat}.{bgen,bgen.bgi}", size:2, flat : true)\
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

process bgenie {
  
  maxForks params.maxForks_SNPTEST
  validExitStatus 0,1,2
  errorStrategy 'ignore'
  // publishDir params.output_dir, overwrite:true, mode:'link'
  
  input:
  set val(chromosome), file(bgenFile), file(sampleFile), val(chunkStart), val(chunkEnd) from bgenChromChunckChannel
  // file bgen_index from bgen_index_ch
  
  output:
  set val(chromosome), file("chr${chromosome}-${chunkStart}-${chunkEnd}.score") into snptestChan

  script:
  """
  ~/software/bgenie/v1.3/bgenie_v1.3_dynamic2 \
  --bgen ${bgenFile} \
  --pheno ${pheno} \
  --range ${chromosome} ${chunkStart} ${chunkEnd} \
  --covar ${covar} \
  --out chr${chromosome}-${chunkStart}-${chunkEnd}.score \
  --thread 1 \
  --pvals
  gunzip chr${chromosome}-${chunkStart}-${chunkEnd}.score.gz
  sed -i '1,1d' chr${chromosome}-${chunkStart}-${chunkEnd}.score
  if [ ! -f "chr${chromosome}-${chunkStart}-${chunkEnd}.imputed" ]; then
  touch "chr${chromosome}-${chunkStart}-${chunkEnd}.score";
  fi
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

process bgenieConcatChunk {
  
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

process bgenieConcatChro {
  
  publishDir params.output_dir, overwrite:true, mode:'link'
 
  input:
  file(chroFiles) from snptestConcatChunkChan.collect()

  output:
  file("${output_pat}.txt.gz") into bgenie_out_ch

  script:
  """
  echo "chr rsid pos a_0 a_1 af info beta se t -log10p" > head
  cat head *.score > ${output_pat}.txt
  gzip -f ${output_pat}.txt
  """
}

process manhattan_qq {

  publishDir params.output_dir, overwrite:true, mode:'link'

  input:
  file(bgenie_out) from bgenie_out_ch
  
  output:
  file("${base}_Manhattan_QQ.png")
  
  script:
  base = output_pat
  template "manhattan_qq_bgenie.py"
}
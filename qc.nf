#!/usr/bin/env nextflow

/*
 * Author:   Jianhua Wang
 * Date  :   09-29-2018
 *
 *
 * This is a pipeline for GWAS QC
 */



inpat = "${params.input_dir}/${params.input_pat}"

f_lo_male = params.f_lo_male
f_hi_female = params.f_hi_female
cut_mind = params.cut_mind
cut_geno = params.cut_geno
cut_maf = params.cut_maf
cut_hwe = params.cut_hwe
pi_hat = params.pi_hat
times_of_meanhet = params.times_of_meanhet

raw_ch       = Channel.create()
bim_ch       = Channel.create()

report = new LinkedHashMap()
repnames = ["dups","basic","snpmisspng","indmisspng","initmaf","inithwe","mafpng","hwepng","misshet","snpmiss","failedsex","misshetremf","pca","related","qc1","qc2"]
repnames.each { report[it] = Channel.create() }

// Checks if the file exists
checker = { fn ->
    if (fn.exists())
        return fn;
    else
        error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}

Channel
    .fromFilePairs("${inpat}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName }\
    .ifEmpty { error "No matching plink files" }\
    .map { a -> [checker(a[1]), checker(a[2]), checker(a[3])] }\
    .separate(raw_ch, bim_ch) { a -> [a,a] }


/*
 * Start pipe
 * 
 * Phase 1: Romove the duplicated variants
 *
 */

// Get the duplicated variants from the .bim files
process getDuplicateMarkers {
    publishDir params.output_dir, pattern: "*dups", overwrite:true, mode:'copy'
    input:
        set file(bed), file(bim), file(fam) from bim_ch

    output:
        file("${base}.dups") into duplicates_ch
        file("${base}.dups") into report["dups"]

    script:
        base = bed.baseName
        outfname = "${base}.dups"
        template "getdups.py"
}

// Romove the duplicated variants and generate basename-nd{.bed,.bim,.fam}
process removeDuplicateSNPs {
    input:
        set file(bed), file(bim), file(fam) from raw_ch
        file(dups) from  duplicates_ch

    output:
        set  file("${nodup}.bed"),file("${nodup}.bim"),file("${nodup}.fam")\
        into (qc1A_ch,qc1B_ch,qc1C_ch)
        file("${base}.orig") into report["basic"]
        file ("${nodup}.lmiss") into snp_miss_ch
        file ("${nodup}.imiss") into (ind_miss_ch1, ind_miss_ch2)

    script:
        base = bed.baseName
        nodup = "${base}-nd"
        """
        plink --keep-allele-order --bfile $base --must-have-sex --exclude $dups --missing --make-bed --out $nodup
        wc -l ${base}.bim > ${base}.orig
        wc -l ${base}.fam >> ${base}.orig
        """
}

/*
 * 
 * Phase 2: check out bad SNPs and sex
 *
 */

// Identify individual discordant sex information
process identifyIndivDiscSexinfo {
    input:
        file(plinks) from qc1B_ch

    publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.badsex"

    output:
        file(logfile) into  report["failedsex"]
        file(logfile) into  failed_sex_ch1
        file("${base}.hwe") into hwe_stats_ch

    script:
        base = plinks[0].baseName
        logfile= "${base}.badsex"
        sexcheck_report = "${base}.sexcheck"
        imiss  = "${base}.imiss"
        lmiss  = "${base}.lmiss"
        """
        plink --keep-allele-order --bfile $base --hardy --check-sex $f_hi_female $f_lo_male --out $base
        head -n 1 ${base}.sexcheck > $logfile
        grep  'PROBLEM' ${base}.sexcheck >> $logfile
        """
}

// sample.lmiss for snp missing rate
process generateSnpMissingnessPlot {
    input:
        file(lmissf) from snp_miss_ch

    output:
        file(snpmiss_plot) into report["snpmisspng"]

    script:
        input  = lmissf
        base   = lmissf.baseName
        label  = "SNPs"
        snpmiss_plot = "${base}-snpmiss_plot".replace(".","_")+".png"
        template "snpmiss_plot.py"
}

// sample.imiss for individual missing rate
process generateIndivMissingnessPlot {
    input:
        file(imissf) from ind_miss_ch1

    output:
        file(indmiss_plot) into report["indmisspng"]

    script:
        input  = imissf
        base   = imissf.baseName
        label  = "samples"
        indmiss_plot = "${base}-indmiss_plot".replace(".","_")+".png"
        template "indmiss_plot.py"
}

// Get MAF of each SNP
process getInitMAF {
    input:
        file(plink) from qc1C_ch

    output:
        file("${newbase}.frq") into init_freq_ch

    script:
        base = plink[0].baseName
        newbase = base.replace(".","_")
        """
        plink --bfile $base --freq --out $newbase
        """
}

// MAF plot
process showInitMAF {
    input:
        file(freq) from init_freq_ch

    output:
        set file("${base}.png"), file("${base}.txt") into report["initmaf"]

    script:
        base = freq.baseName+"-initmaf"
        base = base.replace(".","_")
        template "showmaf.py"
}

// HWE plots
process showHWEStats {
    input:
        file(hwe) from hwe_stats_ch

    output:
        set file("${base}.png"), file("${base}-qq.png") into report["inithwe"]

    script:
        base = hwe.baseName+"-inithwe"
        base = base.replace(".","_")
        template "showhwe.py"
}

// remove really realy bad SNPs and really bad individuals
process removeQCPhase1 {
    input:
        set file(bed), file(bim), file(fam) from qc1A_ch
    publishDir params.output_dir, overwrite:true, mode:'copy', pattern:"*.irem"

    output:
        file("${output}*.{bed,bim,fam}") into (qc2A_ch,qc2B_ch,qc2C_ch,qc2D_ch)
        set file("${base}-QCphase2.out"), file("${output}.irem") into report["qc1"]

    script:
        base=bed.baseName
        output = "${base}-c".replace(".","_")
        """
        plink --keep-allele-order --bfile $base  --mind $cut_mind --make-bed --out temp1
        plink --keep-allele-order --bfile temp1  --geno $cut_geno --make-bed --out temp2
        plink --keep-allele-order --bfile temp2  --maf $cut_maf --make-bed --out temp3
        plink --keep-allele-order --bfile temp3  --hwe $cut_hwe --make-bed  --out $output 
        cat *log > logfile
        touch tmp.irem
        cat *.irem > ${output}.irem
        qc1logextract.py logfile ${output}.irem $cut_mind $cut_geno $cut_maf $cut_hwe ${base}
        """
}

/*
 *
 * Phase 3: check out bad individuals
 *
 */

// PCA
process compPCA {
    input:
        file plinks from qc2A_ch

    output:
        set file ("${prune}.eigenval"), file("${prune}.eigenvec") into pcares

    script:
        base = plinks[0].baseName
        prune= "${base}-prune".replace(".","_")
        """
        plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check
        plink --bfile ${base} --extract check.prune.in --make-bed --out $prune
        plink --bfile ${prune} --pca --out $prune  
        """
}

// PCA plot
process drawPCA {
    input:
        set file(eigvals), file(eigvecs) from pcares

    output:
        file (output) into report["pca"]
        
    script:
        base=eigvals.baseName
        // also relies on "col" defined above
        output="${base}-pca".replace(".","_")+".png"
        template "drawPCA.py"
}

// Get which SNPs should be pruned for IBD
process pruneForIBD {
    input:
        file plinks from qc2B_ch

    output:
        file "${outf}.genome" into find_rel_ch

    script:
        base   =  plinks[0].baseName
        outf   =  base.replace(".","_")
        """
        plink --bfile $base --threads 4 --autosome --indep-pairwise 60 5 0.2 --out ibd
        plink --bfile $base --threads 4 --autosome --extract ibd.prune.in --genome --out $outf
        """
}

// run script to find a set of individuals we can remove to ensure no relatedness
//  Future - perhaps replaced with Primus
process findRelatedIndiv {
    input:
        file (missing) from ind_miss_ch2
        file (ibd_genome) from find_rel_ch

    output:
        file(outfname) into related_indivs_ch1
        file(outfname) into report["related"]
    publishDir params.output_dir, overwrite:true, mode:'copy',pattern: "*.txt"

    script:
        base = missing.baseName
        outfname = "${base}-fail_IBD".replace(".","_")+".txt"
        template "removeRelInds.py"
}

process calculateSampleHeterozygosity {
    input:
        file(nodups) from qc2C_ch

    output:
        set file("${hetf}.het"), file("${hetf}.imiss") into (hetero_check_ch, plot1_het_ch)

    script:
        base = nodups[0].baseName
        hetf = "${base}".replace(".","_")
        """
        plink --bfile $base --het --missing --out $hetf
        """
}

process generateMissHetPlot {
    input:
        set file(het), file(imiss) from plot1_het_ch

    output:
        file(output) into report["misshet"]

    script:
        base = imiss.baseName
        output  = "${base}-imiss-vs-het".replace(".","_")+".png"
        template "missHetPlot.py"
}

// Find those who have bad heterozygosity
process getBadIndivsMissingHet {
    input:
        set file(het), file(imiss) from hetero_check_ch

    output:
        file(outfname) into failed_miss_het
        file(outfname) into report["misshetremf"]
    publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.txt"

    script:
        base = het.baseName
        outfname = "${base}-fail_het".replace(".","_")+".txt"
        template "select_miss_het_qcplink.py"
}

process removeQCIndivs {
    input:
        file(f_miss_het)     from failed_miss_het
        file(rel_indivs)     from related_indivs_ch1
        file (f_sex_check_f) from failed_sex_ch1
        set file(bed), file(bim), file(fam) from qc2D_ch

    publishDir params.output_dir, overwrite:true, mode:'copy',pattern: "*.{bed,bim,fam,irem}"

    output:
        file("${out}.{bed,bim,fam}") into (qc3A_ch,qc3B_ch)
        set file("${out}-QCphase3.out"), file("${out}.irem") into report["qc2"]

    script:
        base = bed.baseName
        out  = "${base}-c".replace(".","_")
        """
        qc2logextract.py $f_sex_check_f $rel_indivs $f_miss_het $out $f_lo_male $f_hi_female $pi_hat $times_of_meanhet
        plink --keep-allele-order --bfile $base --remove ${out}.irem --make-bed --out $out
        """
}

/*
 *
 * Phase 4: Review MAF & HWE and generate report
 *
 */

process calculateMaf {
    input:
        file(plinks) from qc3A_ch

    output:
        file "${base}.frq" into maf_plot_ch
        file "${base}.hwe" into hwe_scores_ch

    script:
        base = plinks[0].baseName
        out  = base.replace(".","_")
        """
        plink --bfile $base --hardy --freq --out $out
        """
}



process generateMafPlot {
    input:
        file input from maf_plot_ch

    output:
        file("${base}-maf_plot.png") into report["mafpng"]

    script:
        base    = input.baseName
        output  = "${base}-maf_plot.png"
        template "mafplot.py"
}

process generateHwePlot {
    input:
        file hwe from hwe_scores_ch

    output:
        file("${base}-hwe_plot.png") into report["hwepng"]

    script:
        input  = hwe
        base   = hwe.baseName.replace(".","_")
        output = "${base}-hwe_plot.png"
        template "hweplot.py"
}

process produceReports {
    input:
        set file(newbed), file(newbim), file(newfam) from qc3B_ch

        file(orig) from report["basic"]
        file(dup_marks) from report["dups"]

        file(snpmisspng) from report["snpmisspng"]
        file(indmisspng) from report["indmisspng"]
        set file(maf_png), file(maf_txt) from report["initmaf"]
        set file(hwe_png), file(hwe_qq_png) from report["inithwe"]
        set file(qc1_out), file(qc1_irem) from report["qc1"]

        file(badsex) from  report["failedsex"]
        file(pca_png) from report["pca"]
        file(fail_IBD) from report["related"]
        file(miss_het_png) from report["misshet"]
        file(fail_het) from report["misshetremf"]
        set file(qc2_out), file(qc2_irem) from report["qc2"]

        file(re_maf) from report["mafpng"]
        file(re_hwe) from report["hwepng"]

    publishDir params.output_dir, overwrite:true, mode:'copy',pattern: "*.html"
    
    output:
        file("${base}-GWAS-QC_report.html") into final_ch

    script:
        base = params.input_pat
        template "generate_report.py"
}

final_ch.subscribe { b=it.baseName; println "The output report is called ${params.output_dir}/${b}.html"}
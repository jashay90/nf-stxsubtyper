#!/usr/bin/env nextflow
nextflow.preview.topic=true

process BLASTN {
	tag "blastn ${query.simpleName} against stx${gene}${subunit}"
	label "multithread"
	input:
	tuple path(query), val(gene), val(subunit)
	
	output:
	tuple val(query.simpleName), val(gene), val(subunit), path("${query.simpleName}_blastout_${gene}${subunit}.tsv"), emit: results
	tuple val("${task.process}"), val('blastn'), eval("blastn -version | head -1 | sed 's/^blastn:\s*//'"), topic: versions

	script:
	"""
	header="qacc sacc qstart qend sstart send qframe sframe slen length mismatch qseq"
	echo \$header | sed 's/\s/\t/g' > "${query.simpleName}_blastout_${gene}${subunit}.tsv"
	blastn -db ${params.basedir}/${params.blastdbprefix}$gene$subunit${params.blastdbsuffix} -query $query -perc_identity ${params.perc_identity} -max_target_seqs ${params.max_target_seqs} -num_threads ${task.cpus} -outfmt "6 \${header}" >> "${query.simpleName}_blastout_${gene}${subunit}.tsv"
	"""
}

process PARSEBLASTN {
	tag "parsing blastn output for ${id} stx${stx}"
	publishDir "${params.outdir}/blast_processed", pattern: "*properpairs.tsv", mode: 'copy'
	input:
	tuple val(id), val(stx), path(ablast), path(bblast)

	output:
	path "${id}_stx${stx}_filtered_A.tsv", emit: filtered // just confirming this exists
	path "${id}_stx${stx}_properpairs.tsv", optional: true, emit: properpairs
	tuple val(id), val(stx), path("${id}_stx${stx}_prot.fasta"), optional: true, emit: protfasta

	script:
	"""
	blastgrab.py -a ${ablast} -b ${bblast} -c ${params.gapcutoff} -l ${params.minlen} -o ${id}_stx${stx}
	"""
}

process ALIGN {
	tag "aligning ${id} Stx${stx} to reference sequences"
	publishDir "${params.outdir}/alignments", pattern: "*_aligned.fasta", mode: 'copy'
	input:
	tuple val(id), val(stx), val(protseq)

	output:
	tuple val(id), val(stx), val(protseq.id), val(protseq.desc), path("${id}_Stx${stx}_${protseq.id}_${protseq.desc}_aligned.fasta"), emit: results
	tuple val("${task.process}"), val('muscle'), eval("muscle -version | head -1 | sed 's/^MUSCLE[[:space:]]//I; s/[[:space:]].*//'"), topic: versions

	script:
	"""
	echo ">${params.qname}" > merged.fa
	echo ${protseq.seqString} >> merged.fa
	cat ${params.basedir}/${params.refprefix}$stx${params.refsuffix} >> merged.fa
	muscle -align merged.fa -output ${id}_Stx${stx}_${protseq.id}_${protseq.desc}_aligned.fasta
	"""
}

process MAKETREE {
	tag "building tree for ${id} Stx${stx} encoded at ${contig} ${loc}"
	publishDir "${params.outdir}/trees", pattern: "*.tree", mode: 'copy'
	input:
	tuple val(id), val(stx), val(contig), val(loc), path(alignment)

	output:
	tuple val(id), val(stx), val(contig), val(loc), path("${alignment.baseName}.tree"), emit: results
	tuple val("${task.process}"), val('fasttree'), eval('fasttree -help 2>&1 | head -1 | sed \'s/^FastTree \\([0-9.]*\\) .*$/\\1/\''), topic: versions

	script:
	"""
	FastTree ${alignment} > ${alignment.baseName}.tree
	"""
}

process CLOSESTLEAF {
	tag "identifying closest leaf for ${id} Stx${stx} encoded at ${contig} ${loc}"
	input:
	tuple val(id), val(stx), val(contig), val(loc), path(tree)

	output:
	tuple val(id), val(stx), val(contig), val(loc), path("${tree.baseName}_closestleaf.txt")

	script:
	"""
	getclosestleaf.py -t ${tree} -l ${params.qname} > ${tree.baseName}_closestleaf.txt
	"""
}

process CHECKMOTIF {
	tag "checking motif of ${id} Stx${stx} encoded at ${contig} ${loc}"
	publishDir "${params.outdir}/motifs", pattern: "*_motif.txt", mode: 'copy'
	input:
	tuple val(id), val(stx), val(contig), val(loc), path(alignment)

	output:
	tuple val(id), val(stx), val(contig), val(loc), path("${alignment.baseName}_motif.txt")

	script:
	"""
	identifymotif.py -i ${alignment} -s ${params.qname} -k ${params.knownseq} > ${alignment.baseName}_motif.txt
	"""	
}

process STXTYPER {
	tag "running stxtyper on ${query.baseName}"
	label "multithread"
	publishDir "${params.outdir}/stxtyper", pattern: "*_stxtyper.tsv", mode: 'copy'
	input:
	path query

	output:
	tuple val(query.simpleName), path("${query.simpleName}_stxtyper.tsv"), emit: results
	tuple val("${task.process}"), val('stxtyper'), eval("stxtyper --version 2>&1"), topic: versions

	script:
	"""
	stxtyper -n $query --name ${query.simpleName} -o ${query.simpleName}_stxtyper.tsv --threads ${task.cpus}
	"""
}

process KMA {
	tag "running KMA on ${pair_id}"
	label "multithread"
	publishDir "${params.outdir}/kma", pattern: "*.res", mode: 'copy'
	input:
	tuple val(pair_id), path(reads)

	output:
	tuple val(pair_id), path("${pair_id}.res"), emit: results
	tuple val("${task.process}"), val('kma'), eval("kma -v 2>&1 | sed 's/^KMA-//'"), topic: versions

	script:
	"""
	kma ${params.kmaopts} -ipe ${reads[0]} ${reads[1]} -t_db ${params.basedir}/${params.kmadb} -o ${pair_id}
	"""
}

process REPORT {
        tag "generating summary report for ${id}"
	publishDir "${params.outdir}/report", pattern: "*.pdf"
        input:
        tuple val(id), path(kma), path(stxtyper), val(stx), val(contig), val(loc), path(closestleaf), path(motif) 

        output:
        tuple val(id), path("${id}_report.pdf")

        script:
        """
        args=(--output "${id}_report" --id "${id}")
        # kma
        if [ -n "$kma" ]; then
                args+=(--kma "$kma")
        fi

        # stxtyper
        if [ -n "$stxtyper" ]; then
                args+=(--stxtyper "$stxtyper")
        fi

	# custom pipeline
        for v in closestleaf motif stx contig loc; do
                case \$v in
                        closestleaf) value="$closestleaf" ;;
                        motif)       value="$motif" ;;
                        stx)         value="$stx" ;;
                        contig)     value="$contig" ;;
                        loc)         value="$loc" ;;
                esac
                if [ -n "\$value" ]; then
                        args+=(--\$v)
                        for f in \$value; do
                                clean=\$(echo "\$f" | sed 's/^\\[//; s/\\]\$//; s/,\$//')
                                args+=( "\$clean" )
                        done
                fi
        done

        latexreport.py "\${args[@]}"
	echo "\${args[@]}" > tmptest
        """
}


// Replaces REPORT: produces a single-row TSV for this sample
process MAKETSV {
    label "pandas"
    tag "making summary row for ${id}"
    publishDir "${params.outdir}/report", pattern: "*_row.tsv", mode: 'copy'
    input:
    tuple val(id), path(kma), path(stxtyper), val(stx), val(contig), val(loc), path(closestleaf), path(motif), val(samplename), val(coverage) 

    output:
    path "${id}_row.tsv"

    script:
    """
    args=(--id "${id}" --output "${id}_row.tsv")

    for v in kma stxtyper samplename coverage; do
        case \$v in
            kma)        value="$kma"        ;;
            stxtyper)   value="$stxtyper"   ;;
            samplename) value="$samplename" ;;
            coverage)   value="$coverage"   ;;
        esac
        if [ -n "\$value" ]; then
            args+=(--\$v \$value)
        fi
    done

    for v in closestleaf motif stx contig loc; do
        case \$v in
            closestleaf) value="$closestleaf" ;;
            motif)       value="$motif"       ;;
            stx)         value="$stx"         ;;
            contig)      value="$contig"      ;;
            loc)         value="$loc"         ;;
        esac
        if [ -n "\$value" ]; then
            args+=(--\$v)
            for f in \$value; do
                clean=\$(echo "\$f" | sed 's/^\\[//; s/\\]\$//; s/,\$//')
                args+=( "\$clean" )
            done
        fi
    done

    xlsxreport.py "\${args[@]}"
    """
}

process WRITE_JSONS {
	input:
	val versions

	output:
	path vjson, emit: vjson
	path pjson, emit: pjson

	exec:
	vjson = task.workDir.resolve('versions.json')
	pjson = task.workDir.resolve('params.json')
	jsonVersions = new groovy.json.JsonBuilder(versions).toPrettyString()
	file(vjson).text = jsonVersions
	jsonParams = new groovy.json.JsonBuilder(params).toPrettyString()
	file(pjson).text = jsonParams
}	
	

// Runs once after all samples are done; assembles the Excel file
process COMBINE_XLSX {
    label "pandas"
    tag "combining all samples into Excel report"
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path(tsvs)
    path(vjson)
    path(pjson)

    output:
    path "summary_report.xlsx"

    script:
    """
    combine_xlsx.py ${tsvs} --params ${pjson} --versions ${vjson} --output summary_report.xlsx
    """
}

params.genomes = null
params.reads = null
params.outdir = "results"
params.makepdfs = false
params.samplenames = null
params.coverage = null

workflow {
	if (params.genomes != null) {
		genome_ch = Channel.fromPath(params.genomes)

		// custom pipeline for checking motif, comparing against our database
		gene_ch = Channel.of('1', '2')
		subunit_ch = Channel.of('A', 'B')
		// run blastn against databases of 1a, 1b, 2a, and 2b subunits
		BLASTN(genome_ch.combine(gene_ch).combine(subunit_ch))
		blastn_ch = BLASTN.out.results.groupTuple(by: [0,1], size: 2).map { id, gene, subunits, files ->
    			def pairs = [subunits, files].transpose()  // [["B", fileB], ["A", fileA]]
    			def submap = pairs.collectEntries { s, f -> [(s): f] }  // makes ["B": fileB, "A": fileA]
    			tuple(id, gene, submap['A'], submap['B'])
		}
		// run summary script for stx1 and stx2
		PARSEBLASTN(blastn_ch)
		prot_ch = PARSEBLASTN.out.protfasta.splitFasta(record: ['id': true, 'desc': true, 'seqString': true])
		// for each protein sequence, align it against reference sequences
		ALIGN(prot_ch)
		// make a tree, identify closest leaf
		MAKETREE(ALIGN.out.results)
		CLOSESTLEAF(MAKETREE.out.results)
		// for stx2, check the motif
		CHECKMOTIF(ALIGN.out.results.filter {id, stx, contig, loc, alignment -> stx == '2'})
		// group results from custom pipeline
		custom_grouped = CLOSESTLEAF.out.join(CHECKMOTIF.out, by: [0, 1, 2, 3], remainder: true).groupTuple(by: 0).map { tup ->
			tup[5] = tup[5].findAll { it != null }
			tup
		}

		// NCBI's stxtyper
		STXTYPER(genome_ch)
		stxtyper = STXTYPER.out.results
	} else {
		stxtyper = Channel.empty()
		custom_grouped = Channel.empty()
	}

	if (params.reads != null) {
		// Read-based detection with kma
		read_ch = Channel.fromFilePairs(params.reads)
		KMA(read_ch)
		kma = KMA.out.results
	} else {
		kma = Channel.empty()
	}

	// optional metadata added from command line
	if (params.samplenames != null) {
		snames = Channel.of(params.samplenames.split(';')).map { row ->
		def cols = row.split(',')
		tuple(cols[0], cols[1])
		}
	} else {
		snames = Channel.empty()
	}
        if (params.coverage != null) {
                cov = Channel.of(params.coverage.split(';')).map { row ->
                        def cols = row.split(',')
                        tuple(cols[0], cols[1])
                }
	} else {
		cov = Channel.empty()
	}

	// combine everything
	allresults = kma.join(stxtyper, by: 0, remainder: true).join(custom_grouped, by: 0, remainder: true).map { tup ->
                // if missing elements because there were no pipeline hits, add empty elements
                def expected = 8
                def padded = tup.size < expected ? tup + Collections.nCopies(expected - tup.size(), null) : tup
                // replace nulls with empty lists (makes nextflow happy)
                padded.collect { it == null ? [] : it }
        }
	everything = allresults.join(snames, by: 0, remainder: true).join(cov, by: 0, remainder: true).map { items ->
        	def sname    = items[-2] != null ? items[-2] : "[]"
        	def coverage = items[-1] != null ? items[-1] : "[]"
        	tuple(*items[0..-3], sname, coverage)
        }

	// make json metadata files
	ch_versions = Channel.topic('versions').unique()
	WRITE_JSONS(ch_versions.collect(flat: false))

	// make a summary Excel file
	MAKETSV(everything)
	COMBINE_XLSX(MAKETSV.out.collect(), WRITE_JSONS.out.vjson, WRITE_JSONS.out.pjson)
	if (params.makepdfs) {
		// generate a summary PDF report
		REPORT(allresults)
	}
}

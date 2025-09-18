#!/usr/bin/env nextflow

process BLASTN {
	tag "blastn ${query.simpleName} against stx${gene}${subunit}"
	label "multithread"
	input:
	tuple path(query), val(gene), val(subunit)
	
	output:
	tuple val(query.simpleName), val(gene), val(subunit), path("${query.simpleName}_blastout_${gene}${subunit}.tsv")
	
	script:
	"""
	header="qacc sacc qstart qend sstart send qframe sframe slen length mismatch qseq"
	echo \$header | sed 's/\s/\t/g' > "${query.simpleName}_blastout_${gene}${subunit}.tsv"
	blastn -db ${projectDir}/${params.blastdbprefix}$gene$subunit${params.blastdbsuffix} -query $query -perc_identity ${params.perc_identity} -max_target_seqs ${params.max_target_seqs} -num_threads ${task.cpus} -outfmt "6 \${header}" >> "${query.simpleName}_blastout_${gene}${subunit}.tsv"
	"""
}

process PARSEBLASTN {
	tag "parsing blastn output for ${id} stx${stx}"
	publishDir "${params.outdir}/blast_processed", pattern: "*properpairs.tsv"
	input:
	tuple val(id), val(stx), path(ablast), path(bblast)

	output:
	path "${id}_stx${stx}_filtered_A.tsv", emit: filtered // just confirming this exists
	tuple val(id), val(stx), path("${id}_stx${stx}_prot.fasta"), optional: true, emit: protfasta

	script:
	"""
	blastgrab.py -a ${ablast} -b ${bblast} -c ${params.gapcutoff} -l ${params.minlen} -o ${id}_stx${stx}
	"""
}

process ALIGN {
	tag "aligning ${id} Stx${stx} to reference sequences"
	input:
	tuple val(id), val(stx), val(protseq)

	output:
	tuple val(id), val(stx), val(protseq.id), val(protseq.desc), path("${id}_Stx${stx}_${protseq.id}_${protseq.desc}_aligned.fasta")

	script:
	"""
	echo ">${params.qname}" > merged.fa
	echo ${protseq.seqString} >> merged.fa
	cat ${projectDir}/${params.refprefix}$stx${params.refsuffix} >> merged.fa
	muscle -align merged.fa -output ${id}_Stx${stx}_${protseq.id}_${protseq.desc}_aligned.fasta
	"""
}

process MAKETREE {
	tag "building tree for ${id} Stx${stx} encoded at ${contig} ${loc}"
	input:
	tuple val(id), val(stx), val(contig), val(loc), path(alignment)

	output:
	tuple val(id), val(stx), val(contig), val(loc), path("${alignment.simpleName}.tree")

	script:
	"""
	FastTree ${alignment} > ${alignment.simpleName}.tree
	"""
}

process CLOSESTLEAF {
	tag "identifying closest leaf for ${id} Stx${stx} encoded at ${contig} ${loc}"
	publishDir "${params.outdir}/trees", pattern: "*_closestleaf.txt"
	input:
	tuple val(id), val(stx), val(contig), val(loc), path(tree)

	output:
	tuple val(id), val(stx), val(contig), val(loc), path("${tree.simpleName}_closestleaf.txt")

	script:
	"""
	getclosestleaf.py -t ${tree} -l ${params.qname} > ${tree.simpleName}_closestleaf.txt
	"""
}

process CHECKMOTIF {
	tag "checking motif of ${id} Stx${stx} encoded at ${contig} ${loc}"
	publishDir "${params.outdir}/motifs", pattern: "*_motif.txt"
	input:
	tuple val(id), val(stx), val(contig), val(loc), path(alignment)

	output:
	tuple val(id), val(stx), val(contig), val(loc), path("${alignment.simpleName}_motif.txt")

	script:
	"""
	identifymotif.py -i ${alignment} -s ${params.qname} -k ${params.knownseq} > ${alignment.simpleName}_motif.txt
	"""	
}

process STXTYPER {
	tag "running stxtyper on ${query.simpleName}"
	label "multithread"
	publishDir "${params.outdir}/stxtyper", pattern: "*_stxtyper.tsv"
	input:
	path query

	output:
	tuple val(query.simpleName), path("${query.simpleName}_stxtyper.tsv")

	script:
	"""
	stxtyper -n $query --name ${query.simpleName} -o ${query.simpleName}_stxtyper.tsv --threads ${task.cpus}
	"""
}

process KMA {
	tag "running KMA on ${pair_id}"
	label "multithread"
	publishDir "${params.outdir}/kma", pattern: "*.res"	
	input:
	tuple val(pair_id), path(reads)

	output:
	tuple val(pair_id), path("${pair_id}.res")

	script:
	"""
	kma -ID 90 -mem_mode -ef -ipe ${reads[0]} ${reads[1]} -t_db ${projectDir}/${params.kmadb} -o ${pair_id}
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

	# custome pipeline
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


params.genomes = null
params.reads = "data/input/*_R{1,2}.fastq.gz"
params.outdir = "results"

workflow {
	if (params.genomes != null) {
		genome_ch = Channel.fromPath(params.genomes)

		// custom pipeline for checking motif, comparing against our database
		gene_ch = Channel.of('1', '2')
		subunit_ch = Channel.of('A', 'B')
		// run blastn against databases of 1a, 1b, 2a, and 2b subunits
		BLASTN(genome_ch.combine(gene_ch).combine(subunit_ch))
		blastn_ch = BLASTN.out.groupTuple(by: [0,1], size: 2).map { id, gene, subunits, files ->
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
		CLOSESTLEAF(MAKETREE(ALIGN.out))
		// for stx2, check the motif
		CHECKMOTIF(ALIGN.out.filter {id, stx, contig, loc, alignment -> stx == '2'})
		// group results from custom pipeline
		custom_grouped = CLOSESTLEAF.out.join(CHECKMOTIF.out, by: [0, 1, 2, 3], remainder: true).groupTuple(by: 0).map { tup ->
			tup[5] = tup[5].findAll { it != null }
			tup
		}

		// NCBI's stxtyper
		STXTYPER(genome_ch)
		stxtyper = STXTYPER.out
	} else {
		stxtyper = Channel.empty()
		custom_grouped = Channel.empty()
	}

	if (params.reads != null) {
		// Read-based detection with kma
		read_ch = Channel.fromFilePairs(params.reads)
		KMA(read_ch)
		kma = KMA.out
	} else {
		kma = Channel.empty()
	}

	// combine everything
	everything = kma.join(stxtyper, by: 0, remainder: true).join(custom_grouped, by: 0, remainder: true).map { tup ->
                // if missing elements because there were no pipeline hits, add empty elements
                def expected = 8
                def padded = tup.size < expected ? tup + Collections.nCopies(expected - tup.size(), null) : tup
                // replace nulls with empty lists (makes nextflow happy)
                padded.collect { it == null ? [] : it }
        }
	// generate a summary report
	REPORT(everything)
}

#!/usr/bin/env nextflow

process BLASTN {
	tag "blastn ${query.simpleName} against stx${gene}${subunit}"
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
	input:
	path query

	output:
	path "${query.simpleName}_stxtyper.tsv"

	script:
	"""
	stxtyper -n $query --name ${query.simpleName} -o ${query.simpleName}_stxtyper.tsv --threads ${task.cpus}
	"""
}

workflow {
	genome_ch = Channel.fromPath('data/input/*.fasta')

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
	prot_ch.view()
	// for each protein sequence, align it against reference sequences
	ALIGN(prot_ch)
	// make a tree, identify closest leaf
	CLOSESTLEAF(MAKETREE(ALIGN.out))
	// for stx2, check the motif
	CHECKMOTIF(ALIGN.out.filter {id, stx, contig, loc, alignment -> stx == '2'})

	// NCBI's stxtyper
	STXTYPER(genome_ch)
}

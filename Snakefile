rule all:
    input:
        auspice_tree = "auspice/enterovirus_d68_vp1_tree.json",
        auspice_meta = "auspice/enterovirus_d68_vp1_meta.json"


rule files:
    params:
        raw_vipr = "data/allEntero-31Jan19.tsv", #raw VIPR download!

        #samples sequenced in Sweden
        swedish_seqs = "data/ev_d68_genomes_2018_vp1.fasta",
        swedish_meta = "data/20190128_Karolinska-region.csv",
        #samples added manually - not from ViPR
        manual_seqs = "data/manual-seqs.fasta",
        manual_meta = "data/manual-meta.csv",

        dropped_strains = "config/dropped_strains.txt",
        reference = "config/ev_d68_reference_vp1.gb",
        blast_ref = "config/ev_d68_reference_vp1.fasta",
        colors = "config/colors.tsv",
        clades = "config/clades.tsv",
        auspice_config = "config/auspice_config.json",
        regions = "config/geo_regions.tsv"

files = rules.files.params

import os.path
RERUN = True if os.path.isfile("genbank/current_vipr_download.tsv") else False

#always parse new genbank files - add regions if wanted
#this makes new strain names from strain+accession to be unique
rule parse_vipr_meta:
    input:
        meta = files.raw_vipr,
        regions = files.regions
    output:
        out = "temp/current_vipr_download.tsv"
    message:
        "This rerun will use existing GenBank files! Only new accession numbers will be downloaded" if RERUN else "Starting new run from scratch. All VIPR samples will be downloaded."
    shell:
        "python vipr_parse.py --input {input.meta} --output {output.out} --regions {input.regions}"


#if new run, ensure to exclude any sequences we already have in 'swedish' or 'manual'
rule remove_dupes:
    input:
        genbank_meta = rules.parse_vipr_meta.output,
        swed_meta = files.swedish_meta,
        man_meta = files.manual_meta,
    output:
        "temp/meta_to_download.tsv"
    shell:
        "python find_new.py --input-new {input.genbank_meta} --exclude {input.swed_meta} {input.man_meta} --output {output}"


#if rerun, find only the new ones to download - exclude old, swedish, & manual samples we already have
rule find_new:
    input:
        old_meta = ancient("genbank/current_vipr_download.tsv"),
        swed_meta = files.swedish_meta,
        man_meta = files.manual_meta,
        new_meta = rules.parse_vipr_meta.output.out
    output:
        out = "temp/new_meta_to_download.tsv"
    shell:
        """
        python find_new.py --input-new {input.new_meta} \
            --exclude {input.swed_meta} {input.old_meta} {input.man_meta} \
            --output {output.out}
        """

#download only new and non-duplicate sequences
rule download_seqs:
    input:
        meta = rules.find_new.output[0] if RERUN else rules.remove_dupes.output[0]
    output:
        sequences = "temp/downloaded_seqs.fasta", 
        meta = "temp/downloaded_meta.tsv"
    run:
        import pandas as pd
        from Bio import Entrez, SeqIO
        from augur.parse import forbidden_characters
        Entrez.email = "richard.neher@unibas.ch"

        meta = pd.read_csv(input.meta, sep='\t')
        additional_meta = {}
        with open(output.sequences, 'w') as fh:
            for ri, row in meta.iterrows():
                try:
                    handle = Entrez.efetch(db="nucleotide", id=row.accession, rettype="gb", retmode="text")
                except:
                    print(row.accession, "did not work")
                    continue
                print(row.strain, row.accession)
                rec = SeqIO.read(handle, 'genbank')
                try:
                    authors = rec.annotations['references'][0].authors
                    title = rec.annotations['references'][0].title
                except:
                    authors = ''
                    title = ''

                url = 'https://www.ncbi.nlm.nih.gov/nuccore/'+row.accession
                additional_meta[ri] = {'url':url, 'authors':authors, 'title':title}
                tmp = row.strain
                for c,r in forbidden_characters:
                    tmp=tmp.replace(c,r)
                rec.id = tmp
                rec.name = tmp
                rec.description = ''
                SeqIO.write(rec, fh, 'fasta')

        add_meta = pd.DataFrame(additional_meta).transpose()
        all_meta = pd.concat((meta, add_meta), axis=1)
        all_meta.to_csv(output.meta, sep='\t', index=False)

#blast all the downloaded sequences
rule blast:        
    input:
        blast_db_file = files.blast_ref,
        seqs_to_blast = rules.download_seqs.output.sequences
    output:
        blast_out = "temp/blast_out.csv"
    run:
        from subprocess import call
        import os
        comm = ["makeblastdb -in", input.blast_db_file, "-out temp/entero_db_vp1 -dbtype nucl"]
        cmd = " ".join(comm)
        os.system(cmd)
        comm2 = ["blastn -task blastn -query", input.seqs_to_blast, "-db temp/entero_db_vp1 -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out", output.blast_out, "-evalue 0.0005"]
        cmd2 = " ".join(comm2)
        os.system(cmd2)

#take only those downloaded that match VP1
rule blast_sort:
    input:
        blast_result = rules.blast.output.blast_out,
        input_meta = rules.download_seqs.output.meta,
        input_seqs = rules.download_seqs.output.sequences
    output:
        out_seqs = "temp/add_sequences_vp1.fasta" if RERUN else "temp/genbank_vp1_sequences.fasta",
        out_meta = "temp/add_meta_vp1.tsv" if RERUN else "temp/genbank_vp1_meta.tsv"
    shell:
        """
        python blast_sort.py --blast {input.blast_result} --meta {input.input_meta} \
            --seqs {input.input_seqs} --out_seqs {output.out_seqs} --out_meta {output.out_meta}
        """

#concatenate meta and sequences to existing genbank, if a rerun:
#concat meta
rule add_meta:
    input:
        metadata = [ancient("genbank/genbank_vp1_meta.tsv"), "temp/add_meta_vp1.tsv"]
    output:
        metadata = "temp/genbank_vp1_meta.tsv"
    run:
        import pandas as pd
        from augur.parse import fix_dates, forbidden_characters
        md = []
        for fname in input.metadata:
            tmp = pd.read_csv(fname, sep='\t' if fname.endswith('tsv') else ',')
            tmp_name = []
            for x in tmp.strain:
                f = x
                for c,r in forbidden_characters:
                    f=f.replace(c,r)
                tmp_name.append(f)
            tmp.strain = tmp_name
            md.append(tmp)
        all_meta = pd.concat(md)
        all_meta.to_csv(output.metadata, sep='\t', index=False)
#concat sequences
rule add_sequences:
    input:
        ancient("genbank/genbank_vp1_sequences.fasta"), "temp/add_sequences_vp1.fasta"
    output:
        "temp/genbank_vp1_sequences.fasta"
    shell:
        '''
        cat {input} > {output}
        '''

#concatenate genbank meta and seqs with Swedish & manual samples
rule concat_meta:
    input:
        metadata = [files.swedish_meta, files.manual_meta, "temp/genbank_vp1_meta.tsv"]
    output:
        metadata = "temp/metadata.tsv"
    run:
        import pandas as pd
        from augur.parse import fix_dates, forbidden_characters
        md = []
        for fname in input.metadata:
            tmp = pd.read_csv(fname, sep='\t' if fname.endswith('tsv') else ',')
            tmp_name = []
            for x in tmp.strain:
                f = x
                for c,r in forbidden_characters:
                    f=f.replace(c,r)
                tmp_name.append(f)
            tmp.strain = tmp_name
            md.append(tmp)
        all_meta = pd.concat(md)
        all_meta.to_csv(output.metadata, sep='\t', index=False)

#concatenate genbank seqs with Swedish & manual
rule concat_sequences:
    input:
        files.swedish_seqs, files.manual_seqs, "temp/genbank_vp1_sequences.fasta"
    output:
        "temp/sequences.fasta"
    shell:
        '''
        cat {input} > {output}
        '''

#if all has gone well up to here, copy relevant files to allow rerun next time!
rule make_database:
    input:
        gen_seqs = "temp/genbank_vp1_sequences.fasta",
        gen_meta = "temp/genbank_vp1_meta.tsv",
        download = "temp/current_vipr_download.tsv",
        seqs = rules.concat_sequences.output,
        meta = rules.concat_meta.output.metadata,
    output:
        gen_seqs = "genbank/genbank_vp1_sequences.fasta",
        gen_meta = "genbank/genbank_vp1_meta.tsv",
        download = "genbank/current_vipr_download.tsv",
        seqs = "results/sequences.fasta",
        meta = "results/metadata.tsv"
    message:
        "Genbank files updated with new sequences!" if RERUN else "Genbank files stored. Reruns will only download new accession numbers."
    shell:
        '''
        cp {input.gen_seqs} genbank
        cp {input.gen_meta} genbank
        cp {input.download} genbank
        cp {input.seqs} results
        cp {input.meta} results
        rm {input.gen_seqs}
        rm {input.gen_meta}
        rm {input.download}
        '''



##############################
# now run usual augur analysis
###############################

rule filter:
    input:
        sequences = rules.make_database.output.seqs, #rules.concat_sequences.output,
        metadata = rules.make_database.output.meta, #rules.concat_meta.output.metadata,
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered.fasta"
    params:
        sequences_per_category = 20,
        categories = "country year month",
        min_date = 1980
    shell:
        """
        augur filter --sequences {input.sequences} --metadata {input.metadata} \
            --output {output.sequences} \
            --group-by {params.categories} \
            --sequences-per-group {params.sequences_per_category} \
            --exclude {input.exclude}  --min-date {params.min_date}
        """

rule align:
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_vp1.fasta"
    shell:
        """
        augur align --sequences {input.sequences} --output {output.alignment} \
            --reference-sequence {input.reference} --fill-gaps
        """

rule tree:
    input:
        alignment = [rules.align.output.alignment] 
    output:
        tree = "results/raw_tree_vp1.nwk"
    shell:
        """
        augur tree --alignment {input.alignment} --output {output.tree}
        """

rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = rules.make_database.output.meta,
    output:
        tree = "results/tree_vp1.nwk",
        node_data = "results/branch_lengths_vp1.json"
    params:
        clock_filter_iqd = 5
    shell:
        """
        augur refine --tree {input.tree} --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} --output-node-data {output.node_data} \
            --timetree --date-confidence --date-inference marginal --coalescent opt \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment,
    output:
        nt_data = "results/nt_muts_vp1.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral --tree {input.tree} --alignment {input.alignment} \
            --output {output.nt_data} --inference {params.inference}
        """

rule translate:
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.nt_data,
        reference = files.reference
    output:
        aa_data = "results/aa_muts_vp1.json"
    shell:
        """
        augur translate --tree {input.tree} --ancestral-sequences {input.node_data} \
            --output {output.aa_data} --reference-sequence {input.reference}
        """

rule clades:
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.aa_data,
        nuc_muts = rules.ancestral.output.nt_data,
        clades = files.clades
    output:
        clade_data = "results/clades_vp1.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output {output.clade_data}
        """

rule traits:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.make_database.output.meta,
    output:
        node_data = "results/traits_vp1.json",
    params:
        columns = "country region"
    shell:
        """
        augur traits --tree {input.tree} --metadata {input.metadata} \
            --output {output.node_data} --confidence --columns {params.columns}
        """

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.make_database.output.meta,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.nt_data,
        aa_muts = rules.translate.output.aa_data,
        colors = files.colors,
        clades = rules.clades.output.clade_data,
        auspice_config = files.auspice_config
    output:
        auspice_tree = rules.all.input.auspice_tree,
        auspice_meta = rules.all.input.auspice_meta
    shell:
        """
        augur export --tree {input.tree} --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.clades}\
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} --output-meta {output.auspice_meta}
        """


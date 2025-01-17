lyve-KSNP
=========

A set of wrapper scripts for KSNP

Install
-------
1. Clone with git
2. Update your PATH to the git directory

Updating
--------
Use git pull -u

Usage
-----

1. Create a ksnp database with `lyve-manage-ksnp.pl`

        Usage: lyve-manage-ksnp.pl *.fastq[.gz] -d database.fasta
        NOTE: you can also add fasta files that have already been converted; however you cannot add already-merged fasta files.
        -d database Of merged fastas, produced by ksnp executable merge_fasta_reads
                    The database can also be a blank or nonexistant file that this
                    script can create.
        --action indicates one of the following: add, add-assemblies, remove, repair, query. Default: query
          add:            add reads
          add-assemblies: add assembly
          remove:         remove entry
          repair:         recreate the ksnp fasta file index files
          query:          find out if an entry is present
        -t tmp/

2. Run KSNP with `lyve-ksnp.pl`

          Usage: lyve-ksnp.pl -d database.fasta -o outdir/ -t tmpdir/
          -kmerlength kmer length (default: 25)
          -keep to keep temporary files
          -q qsubxopts Any extra option you want to pass to qsub, e.g. -q '-q long.q'

Output files
------------
1. The fasta file that you specify in `lyve-manage-ksnp.pl`, e.g., db.fasta
2. The database files for Lyve-KSNP
  3. `db.fasta.reads`: the raw reads filenames that KSNP needs
  4. `db.fasta.assemblies`: the finished assemblies filenames that KSNP needs
  5. `db.fasta.index`: the index file that Lyve-KSNP needs
  6. `db.fasta.locked`: This zero-byte file will appear when running Lyve-KSNP to avoid concurrency
7. out/  The KSNP output directory

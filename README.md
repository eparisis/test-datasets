# test-datasets: `cyanosuite` branch

Test data for [eparisis/CyanoSuite](https://github.com/eparisis/CyanoSuite). Each pipeline gets its
own branch. The pipeline reads these files with
`pipelines_testdata_base_path = 'https://raw.githubusercontent.com/eparisis/test-datasets/cyanosuite/'`.

Sources, regions and results are described in CyanoSuite's `tests/data/README.md`. The datasets
are rebuilt with `make_test_datasets.sh` (copied from there). MD5s are in `md5.txt`.

| File | Used by |
|---|---|
| `mini_samplesheet.csv`, `mini_coculture.fastq.gz`, `mini_cyano.fna.gz`, `mini_noncyano.fna.gz` | `-profile test_mini` |
| `mid_samplesheet.csv`, `mid_coculture.fastq.gz`, `mid_nostoc_ms1.fna.gz` | `-profile test_mid` |
| `host_human_chr1_fragment.fasta` | `decontamination_reference` in both profiles |
| `gunc-mock.dmnd` | `gunc_db` in `test_mini` |

Every file is under 100 MB, so everything lives in git on this branch.

# test-datasets: `cyanosuite` branch

Test data for [eparisis/CyanoSuite](https://github.com/eparisis/CyanoSuite). Each pipeline gets its
own branch. The pipeline reads these files with
`pipelines_testdata_base_path = 'https://raw.githubusercontent.com/eparisis/test-datasets/cyanosuite/'`.

Sources, regions and results are described in CyanoSuite's `tests/data/README.md`. The datasets
are rebuilt with `make_test_datasets.sh` (copied from there). MD5s are in `md5.txt`.

| File | Used by |
|---|---|
| `mini_samplesheet.csv`, `mini_coculture.fastq.gz`, `mini_cyano.fna.gz`, `mini_noncyano.fna.gz` | `-profile test_mini` |
| `mid_samplesheet.csv`, `mid_nostoc_ms1.fna.gz` | `-profile test_mid` |
| `host_human_chr1_fragment.fasta` | `decontamination_reference` in both profiles |
| `gunc-mock.dmnd` | `gunc_db` in `test_mini` |

`mid_coculture.fastq.gz` (313 MB) is over GitHub's 100 MB file limit, so it is not in git. It is
attached to the release `cyanosuite-mid-v0.1.0` and `mid_samplesheet.csv` points at its download
URL. Locally it sits in the git-ignored `release-assets/` folder until it is uploaded.

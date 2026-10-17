# Saved inputs for the worked examples (GENCODE 44 and 50)

For LEPR, NTRK2, NTRK3 and FLT1: the config `annotate --gtf` proposes from GENCODE 44 (Ensembl 110, the release of the
choroid-plexus study's index) and from GENCODE 50 (Ensembl 116), and the sequence `identifiability --gtf
--transcripts-fasta --save-inputs` saved from the same files, the configured transcripts and the gene background.
`LEPR.gencode44.study.config.json` holds the study's classes (LepRa as `short`: ENST00000371060, ENST00000616738; LepRb
as `long`: ENST00000349533) and runs on `LEPR.gencode44.inputs.json`, which holds all eight LEPR transcripts of GENCODE
44. Any of them repeats with no network:

```bash
isoform-dominance identifiability --config LEPR.gencode50.config.json --inputs LEPR.gencode50.inputs.json
```

Made on 2026-10-16 with isoform-dominance 2.7.0, per gene G and release N:

```bash
isoform-dominance annotate --gene G --gtf gencode.vN.annotation.gtf.gz --out G.gencodeN.config.json
isoform-dominance identifiability --config G.gencodeN.config.json --gtf gencode.vN.annotation.gtf.gz \
    --transcripts-fasta gencode.vN.transcripts.fa.gz --save-inputs G.gencodeN.inputs.json
```

from these GENCODE files (each inputs file records them under `annotation_source`):

| file | bytes | SHA-256 |
|---|---|---|
| `gencode.v44.annotation.gtf.gz` (2023-03-01) | 49,721,965 | `01f817afed65feee863361b4baf30a95e722ee5c5d508ff77b04106ef7ba20d3` |
| `gencode.v44.transcripts.fa.gz` | 82,601,810 | `50a4143885523829c2fc8b1399b97c136ead8f313457165b855d63751a4e63f4` |
| `gencode.v50.annotation.gtf.gz` (2026-04-08) | 124,527,720 | `83fba3e9b03f0b8c958f3595c6c350adc55f468abf8b0e47b6d5284cfe13a453` |
| `gencode.v50.transcripts.fa.gz` | 183,554,921 | `5a320f524d73b5793518eb19b118829033713443d0f42af20a67bb31cc06cf56` |

| file | holds | SHA-256 |
|---|---|---|
| `FLT1.gencode44.config.json` | config: `iso_733aa` (1) vs `iso_1338aa` (1) | `0bf39bf98d6bacb269d24e44c423488efcb9d5917a8151985cdbc84000970798` |
| `FLT1.gencode44.inputs.json` | saved inputs: 2 configured + 8 background transcripts, release 110 | `5602414a239b0af2fd66c3d06152540c5efce1e6dcd01552348b005bca4340fd` |
| `FLT1.gencode50.config.json` | config: `iso_733aa` (1) vs `iso_1338aa` (6) | `50d9082962e13056b60d9e15d0c372d41191e9d8b30ee1493baec286a4110e12` |
| `FLT1.gencode50.inputs.json` | saved inputs: 7 configured + 10 background transcripts, release 116 | `6301e569b7f6f2a05f7a1e0dffaf9a2340425afb8b276bb753ca8adfa146d16e` |
| `LEPR.gencode44.config.json` | config: `iso_906aa` (2) vs `iso_1165aa` (1) | `0c88c1c01b228931120b5f5b92f2e4a875202dbece6fcdbb4818e0c4fc8e9c03` |
| `LEPR.gencode44.inputs.json` | saved inputs: 3 configured + 5 background transcripts, release 110 | `e6c7c44d4b4dc2ff0ddead4a3d1535661b7a99eb849591326989844f15755133` |
| `LEPR.gencode44.study.config.json` | config: `short` (2) vs `long` (1) | `57fac9687b46485f81a3d16225c669fcaa02efd7191738fd9ad531af10552517` |
| `LEPR.gencode50.config.json` | config: `iso_896aa` (7) vs `iso_1165aa` (2) | `4e11c9cf38d59ef5ee1a4b6ef80a49837d5114bf18c00c5e20b59b8a2b9f6ddb` |
| `LEPR.gencode50.inputs.json` | saved inputs: 9 configured + 9 background transcripts, release 116 | `6239cdbb56288b324a58982d36398a699fde5686500869338d05756ed8946a1b` |
| `NTRK2.gencode44.config.json` | config: `iso_477aa` (18) vs `iso_822aa` (10) | `3bead80a76f0354f930a272d0457c840a0e0f7e9500e4e467331c4f3b5e731e1` |
| `NTRK2.gencode44.inputs.json` | saved inputs: 28 configured + 30 background transcripts, release 110 | `d023a6562a7bfe9b4b813fe6423f83c755bf71323c8c83a4c0f2ad4a233a5250` |
| `NTRK2.gencode50.config.json` | config: `iso_477aa` (22) vs `iso_822aa` (12) | `1b25c460365cd0e2c2054002497e79626aad54c1ede90ed3d6a2f0f61159f4a7` |
| `NTRK2.gencode50.inputs.json` | saved inputs: 34 configured + 30 background transcripts, release 116 | `0dcaf64c4c9b277f4b71da54026ed8294a41799c58b50f258f1c103004ed0619` |
| `NTRK3.gencode44.config.json` | config: `iso_617aa` (2) vs `iso_839aa` (6) | `29ab172b1e75ea03418362e6e097909c7d9cc9ef39891cb310193ea6888595d4` |
| `NTRK3.gencode44.inputs.json` | saved inputs: 8 configured + 11 background transcripts, release 110 | `e07ce5343ab02f22b3067f80f4ce91013972bd4c611a5dc3f8825f184e2812b4` |
| `NTRK3.gencode50.config.json` | config: `iso_604aa` (7) vs `iso_825aa` (20) | `e777bc52b82838a531e98c93fd9e3d468f9d31ef9d593e4fe03202cf822a5c22` |
| `NTRK3.gencode50.inputs.json` | saved inputs: 27 configured + 11 background transcripts, release 116 | `947886045b25272ca02d987c56b6f685d778ea2f5c9373e3c85cd0bcf7c12034` |

`tests/test_saved_identifiability.py` holds the `identifiability --json` report on each of the nine runs to the one
2.7.0 gave (`tests/data/v270_identifiability`), and `tests/test_readme_example.py` runs the README's identifiability
example on `LEPR.gencode50.*`. `isoform-dominance selftest` runs the study config on `LEPR.gencode44.inputs.json` from
byte-for-byte copies in `src/isoform_dominance/data/`, which the wheel carries; `tests/test_selftest.py` checks that
they are these files.

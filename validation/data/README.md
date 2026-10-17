# Stored simulation outputs and inputs (issue #2)

The outputs the README's simulation figures come from, and the inputs they were made from. They were made on 2026-09-11 to
2026-09-14 with `../run_simulation.py` and `../simulate_reads.py`, outside this repository, and are kept here unchanged.

## Design

- **Genes and sequence.** 49 genes. For each, the config `annotate` proposed against Ensembl 116 on 2026-09-11
  (`rest.ensembl.org`) and the cDNA of every transcript of the gene, configured or not, from the same release.
- **Reads.** One library holds every gene. Each transcript is given TPM 10 at 30M fragments (`--depth-scale 1`, mean
  effective length 1500), so a gene's fragment count grows with its number of transcripts. Fragment lengths are drawn
  from N(200, 60) and clipped to [100 nt, the transcript's length], fragments are placed along the transcript, 100 nt
  are read from each end, and each base is substituted with probability 0.002. Positional coverage skew `bias_5p`: above
  0 toward the 5' end, below 0 toward the 3' end. 30 replicates per design; the seed of replicate r is `r * 7919 + 13`,
  and the skew moves where fragments start, not how many there are, so every design draws the same fragment counts.
- **Quantification.** Salmon 2.7.0 (bioconda's Rust rewrite of Salmon; the version number is Salmon's, not this
  package's), `salmon quant -l A` (the library type was detected as ISF), against an index (k = 31, no decoys,
  duplicates removed: 453,553 of 465,769 records) of Ensembl 116's `Homo_sapiens.GRCh38.cdna.all.fa.gz` (183,898,799
  bytes, SHA-256 `683eb19310c40bf1396e4718f45afa2ce86755717c0990f47a171f535d248ea1`). At least 99.995% of fragments
  mapped in every replicate.
- **What is compared.** Per gene and replicate, log2 of the ratio of the two classes' summed `NumReads` against log2 of
  the ratio of the fragments actually drawn for them; a replicate counts when all four are positive. Per gene and design:
  `n` (replicates that count), `realised_sd` (SD of the estimated log2 ratio), `bias` (median error), `rmse`, and what
  the package predicted at the same design before any read existed: `verdict`, `estimable`, `conditioning`,
  `min_resolvable_log2fc`, `beyond_linear`, `predicted_se` (= `min_resolvable_log2fc` / 1.96). A gene with fewer than
  three replicates that count has no row.

## Files

| file | `bias_5p` | coverage | genes | SHA-256 |
|---|---|---|---|---|
| `sim/c1_bm20.json` | -2.0 | 3' | 47 | `6512902fad47b118e9932ddf8ea9a4858197a59f559512fe4a65dba0ebfd0cfb` |
| `sim/c1_bm10.json` | -1.0 | 3' | 48 | `d3338d3f7c86ee78e486d0f7115321bb9d80a32475e2de7ce158da1753f5e428` |
| `sim/c1_bm03.json` | -0.3 | 3' | 47 | `9cf7edc90872a95f29c6e614398256c4fe99af9160784524a613567549b2e826` |
| `sim/c1_d1.json` | 0 | uniform | 49 | `a5906a171869cbac8a5c2791591569be78a473281cb6d92c8e1bdd6e11756c3c` |
| `sim/c1_b03.json` | +0.3 | 5' | 49 | `381bef5c77b84b0710bd0a962026c67ab8bb67b7069d2aaad817e42918f23fe1` |
| `sim/c1_b10.json` | +1.0 | 5' | 48 | `80550bd8c3a33a66aa68868c43e1ccc66d1987c2dde76b9c3b64ed3e06e20bf4` |
| `sim/c1_bias5p.json` | +2.0 | 5' | 49 | `ec07605886db08b841931256da3e1e16fe6a0a1cfd9f97b22fb2255ee5e03a23` |
| `sim/chk_d1.json` | 0 | uniform, run again after the harness gained a second quantifier: every gene's row is identical to `c1_d1.json` | 49 | `6f0038cb4208e43e16622752289705cd5ef9fe5c575fb7d22efb58f9b8f766aa` |

The same reads quantified by kallisto 0.52.0 (`kallisto quant --fr-stranded --plaintext`, the library type asserted; an
index of the same FASTA at k = 31, the harness's default; at least 99.6% of fragments pseudoaligned). kallisto has no
positional or sequence bias model.

| file | `bias_5p` | coverage | genes | SHA-256 |
|---|---|---|---|---|
| `sim/k1_bm20.json` | -2.0 | 3' | 47 | `8db3a66bfe05d94612e4914e7718afae2c1702130dd209b7608b1ee11cc4f800` |
| `sim/k1_bm10.json` | -1.0 | 3' | 47 | `456c95456377f0a114c4368ebd3b582b0f05630afa7a8294dcba973772473963` |
| `sim/k1_bm03.json` | -0.3 | 3' | 48 | `94fdbc386fe5aacbb28598e94736b1558a1ae164d37b489df153fd3fad34a45e` |
| `sim/k1_d1.json` | 0 | uniform | 48 | `a21c917d2aef4218d4d200ce8cd022193ef9d1166971cebe3fa149bee29a2043` |
| `sim/k1_b03.json` | +0.3 | 5' | 48 | `fd87266a31d9581d92312e77e515df6858a8fd844ce1cbef8e9fe3fdbda6e16f` |
| `sim/k1_b10.json` | +1.0 | 5' | 48 | `bdc23e79af1c2c77ab7daf3c76d7f412464efd985049005292fe8cb3c90cac07` |
| `sim/k1_b20.json` | +2.0 | 5' | 48 | `f535c0027f21d2b2316789158286929d838bf6ae9de83c8539bb1f31c66125b7` |

| file | holds | SHA-256 |
|---|---|---|
| `sim/class_sizes.json` | per gene, over the 30 replicates of the `bias_5p` = -2.0 design: the median fragments drawn for each class of the primary comparison (`nA`, `nB`), the median true log2 ratio (`true`), the median fragments of the gene (`total`), and how many replicates had a defined true ratio (`n_true`) | `7c83c78a3a28f41a254f69905aa83958b17ca764935fc01deb96f973fff312dc` |
| `sweep116.json.gz` | gzipped JSON `{"format", "source", "genes": {gene: {"config", "sequences"}}}`: each gene's config and `{transcript_id: cdna}` exactly as the simulation read them, from one file per gene each (SHA-256 of those files below) | `ebbb3654cb1ec816eb71e561aa2481f681a15654db339eb6f4614fc7d565455b` |

`../summarize.py` reads them and recomputes every figure the documentation quotes from them, and `tests/test_validation.py`
holds the files to the SHA-256 above, the package's predictions to the ones the rows recorded, and the summary to them.
Salmon is not run in CI.

<details><summary>SHA-256 of the 98 per-gene files bundled in <code>sweep116.json.gz</code></summary>

```
8755a7f3312c73b25cf97f1944c90abba4610792d3510464cd4dbc9cc5748d5b  sweep/ALDOA.cfg.json
bbd1cbddad7c1238b008b07746bbaadfa6529093884226b1e0aea370737dad90  sweep/ALDOA.seqs.json
457a5da1daac83ca42955b4d61e4e60cc1b96d0004806e13c0a54aee91c1fcf4  sweep/AR.cfg.json
a7a8b3a8b564d2b47261703963cb44554b073e08c7b7fdcca096b7e292a2ce1d  sweep/AR.seqs.json
84745f95a9eb36f7c59bd91455a64d836ff94d6bfb250aeec820a1c11baea498  sweep/AXIN1.cfg.json
83defcd7658f478864268790272ad901af35cca6978226b54e0788340ab30a4f  sweep/AXIN1.seqs.json
a0e61b692f70f346b1d8e2fdcce478ff91da7d39770877aeb1c6d1a7b678791a  sweep/CANX.cfg.json
4ed62a31d9efbb4c47e1d3cb51c1fd4d79e962b47ef70f25208775af762df5a3  sweep/CANX.seqs.json
127ccf70924ede5e52cf1067d4ae5a4b217360b8689fa834771d8e5cf8a97cc0  sweep/CASP9.cfg.json
61b21afb1684c037ca57e2e4718bc2ec89ddcbbe0fd7f6d52d850c8b9f291c9d  sweep/CASP9.seqs.json
ca33b3198ced55b860962ce394f788e843ed3f5f1b4bf138e6b2c99fba120b29  sweep/CCND1.cfg.json
51b2d46989b71f789ed8b1eda2443f911ea3d036b72c218c80e791a35012879a  sweep/CCND1.seqs.json
c982fd6f86592916832f7ee8698bacf92f952f0792ad9b4a8a38a38dbda88911  sweep/CD40.cfg.json
3ce988ae81df63806c75d6bcda920d042d0b8979a062b03d51fd0654b3b98889  sweep/CD40.seqs.json
db99d6f06426b25323b75d45bb00191d43010e54f1dfc57c6d7c123b06c2aa02  sweep/CD44.cfg.json
7d4b00e650070224422c89fa6232dbff1890f826ea84e6e2ad55595fc4ec82e6  sweep/CD44.seqs.json
2090d7556413befba63ca47fd2243654cb3eb7fcc9d60ff02a5afc479ee79c47  sweep/COX5A.cfg.json
2dce7185807ccec6df81f2d84687dd087dcbc4a59a1b719e52f206922c9422f6  sweep/COX5A.seqs.json
63c1fa0f1ffc6ae672072885389951fb40197865e68e67cd497791d85c9eafa9  sweep/CS.cfg.json
ab3c7aae8c25716da02691e8f227f1667680b5f0c0483e93dd4173e4e166279f  sweep/CS.seqs.json
09d178da3aa302705c3fda39e027a38575ae9e4db79e768195d97be1437aa33f  sweep/CSF1R.cfg.json
5bfa557a5a62614dc367332db91b5b67ce729f00ace486d9c5d274cdd357987b  sweep/CSF1R.seqs.json
cded84d1fa8d9ec5223827ed4084d29582711a16feb1ad5e9ac8e6e2af28a8a1  sweep/DMD.cfg.json
70029bb535a90068acacc766f86e51df8f77ec2d85c11afb4ac9984200446f8b  sweep/DMD.seqs.json
92f7c5e1386561287c244bd2f84fa8336b702016313d5742c3589e10ab4b644a  sweep/ECHS1.cfg.json
a5ada04f779db8dec7cd26d277eb8f8f4469bba71c839affbae85e729be2b65e  sweep/ECHS1.seqs.json
ff1ef3f41e403a555316e50b357d0dfdb7c5c5ee7f00e3a9b4d78f3b4c1e4b63  sweep/EPOR.cfg.json
00317f627360484e4eb4dcb26cf6be3eb7ef49cd9dc02b7ac832046d91d090af  sweep/EPOR.seqs.json
21dad3d807874246d1e5af2bc4a2830d465cd467815823eeb7f32d83f23cde82  sweep/ERBB4.cfg.json
123293eb03d2c09344582f815cfae6e4d89391b8a903c8bdb10d615019c3d13e  sweep/ERBB4.seqs.json
678e581eb219fb050bd427b786d7b81ecc25d7dc688ae6f7df357ad0f450da5e  sweep/ESR1.cfg.json
170fb78bc854b0e1c99136f2f31cc708e924b619010f3d4b590a28ad821fcfa8  sweep/ESR1.seqs.json
723a03a16ed2f756e9fb8d8f994807b09be4e0f022cda4d79bd322a6e9e99730  sweep/FASN.cfg.json
05e54a68685dac3277bdf283ee8c405564bf0c8b9171deda650e21ff59c9a0a8  sweep/FASN.seqs.json
bd251143998b928123704577f21766f8b93428190c85ec2333fb7cdea2a61e69  sweep/FLT1.cfg.json
e25949bf50337e19f880f806c07a9f11d4c6b238ffb536f16f3ddc6a493e7130  sweep/FLT1.seqs.json
e7674929888ac0219c66fafa3cbff0358b53121720f5fb7bd354be2c7643fb41  sweep/FOXO1.cfg.json
db565918f6d166e0e3d4cb1bc797dc70bc923dbee0c35178c83e212062144f14  sweep/FOXO1.seqs.json
c813039fdefad55a9bc707027561421ec25e27af30e5ce1aab99b1ec2a96647b  sweep/GHR.cfg.json
40699110340a117bd7d177e60ec3a494737bf7769c15bac3c8adfd88d999689b  sweep/GHR.seqs.json
4e74a45cdb9d3c2643a410939720d5c43b1d4aa75e9633d801629c221e0fbe84  sweep/GNAS.cfg.json
fc8373cab7f181dd4ab01473d977ec30a6fee646ab90b4807ad5e643e016a0c6  sweep/GNAS.seqs.json
7d6263ba5e0e9f53740ee2272efdaaa0fd609241316111ef3fda479bd1b00601  sweep/HK1.cfg.json
7fa2a001f75bfc45b0cb6553b2c2a9abd4bc9fb65d94f1e03529310834f71275  sweep/HK1.seqs.json
7521444f33994ac156fe68f3f9e1b5e83b4ca40c2b8b9c808c6763769c01c6d8  sweep/IGF1R.cfg.json
8fe389bdf4ebde25ef32b220a40ad4b2f6f728652161b9b615193b4ec434cf8b  sweep/IGF1R.seqs.json
505ea2873b58c3b36b555f50cb574acc45c7e22e56aed1b9dd84f81046806a9e  sweep/IL1R2.cfg.json
f2248f20908088ded5fc436dfac3c2e7cb4e6af14c07a948ed9c2440930caaeb  sweep/IL1R2.seqs.json
23bfccdcacaa960b809de8775151f1434dd827ee5855a2b9a386fc5bb734275b  sweep/IL6R.cfg.json
d4f92ccb7df0869046f9226d87edcc089b81ec6cdece8ba7f9655316caef9ba5  sweep/IL6R.seqs.json
491b007e180861ebfbe0b4f1d4a0493b0822fea17588c24830aa075db4eaa9bc  sweep/INSR.cfg.json
114383fabb6ad8affa6c701e1d09adcfd1d1a4ca082a25e70ff26a821638bdf2  sweep/INSR.seqs.json
d95d615f149cecac3b0cafb71a98771dde646f6457d355208f0a6a6bf0c10a02  sweep/JAG1.cfg.json
4e11cb0dcea9cb99f13fa809d219eb60f764d2616d3deaa771daf6b3b30478f2  sweep/JAG1.seqs.json
953be824bb48eded9375cdb5158f5878b1b2989f83a9549f2da4b4d9c86cb286  sweep/LDHA.cfg.json
2e3fa777ea2cee9019219b7b0b2388cafa2513e0a519ba6fd3e4f7eeb193b18b  sweep/LDHA.seqs.json
5b586c5853b55e039573af9e29abcd04e5d046322e04751540318d6fd4c6375a  sweep/LEPR.cfg.json
0f0c5872616bf58f440a088ba7788942a14696aedf208f21f0e5e0efe3f340a4  sweep/LEPR.seqs.json
747b65eb12f14394547e91089ed01f017b2307dd9d10164aeb36b425947936e5  sweep/MCL1.cfg.json
fb24fc080be13050f6364b23a998fc022bf282e79cc95e8c22f2bb2f0ec4afdf  sweep/MCL1.seqs.json
d098e8ae886357084226d70be5afd925136972fdb0aa64802efcaf9e459b43aa  sweep/MET.cfg.json
d5f36b090cee54437626f40032f7e478119740b7f9b1d3ca6b4cb30fce8f3f42  sweep/MET.seqs.json
4baea3ea42bac1367cb2dccbd56536aa0ed603ddedcfafb82deef8579b3d4543  sweep/NDUFS1.cfg.json
de36c0124807c457141cf0375dd5c74652917b52647f6d0840db4302209a445f  sweep/NDUFS1.seqs.json
2e06ee66d45129e335f108dbca1e58695f6db01b72918732875336810d0cf00a  sweep/NR1H3.cfg.json
098ff7e5d77169117f0ff003cc5ba075216793995a7edf2fafc8bc6176ed928f  sweep/NR1H3.seqs.json
9679013ca14409cc13e5c554efb77e549c644653a1a530710536ee6c981c02f3  sweep/NTRK2.cfg.json
2e04bddb959dfb4b6fc50e8a7e430c8623c3b1eebd1c16fa9e6433347028a662  sweep/NTRK2.seqs.json
03e3d39f90ca3cf6185a8c020dee6735271d79c21f82a20f5df2ad60e84cf707  sweep/NTRK3.cfg.json
d0cf7dcc82bf5dd3760910e0a2dffc907c6f40e0e6b8a34f67d2e0765a6b1320  sweep/NTRK3.seqs.json
eb8ef4bad6fc17021086b7519cdf1670aaf304b38f2170a9f78fe99da4d9a638  sweep/PFKL.cfg.json
f0b179826bad39a38cab44262bdf33748f5b8a880756ff7960cbea807e0a52dd  sweep/PFKL.seqs.json
11c7998859ab2e50489bf4a80d4ba02d707651951150d387d0ced3bfd020a5a1  sweep/PIK3CA.cfg.json
31ec558702af39bf3588cc2a103af398499c6fced4cf21f5b574df4a8b0afa80  sweep/PIK3CA.seqs.json
8bf5e040cb75a51df2fb540fd55e444ed2b12ddaf8435c344adae00ef426d6b4  sweep/PPARG.cfg.json
e6a150c77edb0c6f85c52ce11455e7f54aeac66894e9cbc9a9c79f023e3f7142  sweep/PPARG.seqs.json
389acf16502c5573ec292d518af608b7a2b89ccd6ce436aa6cf20aac886ef8c2  sweep/PTPRC.cfg.json
965f810e8e98f8bb3b531592c9dde8004a1bf09813cafc9279fd121b32585e52  sweep/PTPRC.seqs.json
a69dd67259d3afc84f4fedba89ef21fb7f08ff280c5dcd2f38c8912f5b9249c7  sweep/SLC2A1.cfg.json
ea5cd25ff1f1777a17e62d859cb26847b731c512a4a0e77124ac64a3078ba0c1  sweep/SLC2A1.seqs.json
727d1fb59675dc4f63a1de17494d998a8256aec0db1f6c6d343f9cfb8e121559  sweep/SMN2.cfg.json
d26e852815ed97569d8282105a938a09a546b22d49176fb89b36c2c0390bc952  sweep/SMN2.seqs.json
48ad08476f4ded7907e81bb81011b83fd9614edb9d47fdbf76fca46e54fa2d8e  sweep/SREBF1.cfg.json
92ba141bbbff243eb9db0560ed964d3c4a3922dec64b4e014dfe380cb60a02d2  sweep/SREBF1.seqs.json
1a60a9c6cecb98a24f867f4ce362b7a36a9528627dec9989b55086cb0486a440  sweep/TGFB1.cfg.json
7442c54ce7eea67292610fb6d84c76cdb92277ddaa87bbd346fb765d7a8f654d  sweep/TGFB1.seqs.json
694f39a6af60faa6cf0a8542765e0117b35e21cba687ba8b4e9497cdc9388228  sweep/TNFRSF1A.cfg.json
7f70d306aeb7694568f0daac27a2618985bb442ec59e250aa012f73320b004c6  sweep/TNFRSF1A.seqs.json
9663a1591673f54d5d787ca66f04769c21634efb0abcabb89a5c39ae71e082bc  sweep/TNFRSF1B.cfg.json
7fb5e786f60fd496a671c85f2ca9268fb9c8a8617de36334c15f20bac4cf1236  sweep/TNFRSF1B.seqs.json
3efa27be2e4fb3f92123608d29a6d836020f5f7a48c177303fff63b934eff5d0  sweep/TP53.cfg.json
8c6e89853b9a86f6a7a0b16abde1144f57318ead286a233ec61581bb1649d616  sweep/TP53.seqs.json
4d87fd86d707a979f03879d1084996b33127ec4234b41c879985d97a12e44f44  sweep/TPI1.cfg.json
09824cef0b43abc5c5c937e5335d080afc384fd64772dbfa7eceba3d37de44a0  sweep/TPI1.seqs.json
59bb877d7d4412584d901e13abdb2fd7474af172388c1c6be6ac8cdeede2568e  sweep/TSC1.cfg.json
ea10edd2206000932274a30668f96cf25b16000a688bb8fe0e03adae4fbb2e0c  sweep/TSC1.seqs.json
1399074c7af1748b3d69b24fe01c4d937bcd067bed3d2472896be0aa8365bc07  sweep/VEGFA.cfg.json
adf9fb767b5b33be1aa7f11d78c68d1aacd45437539d95dec14272b83c0984ba  sweep/VEGFA.seqs.json
```

</details>

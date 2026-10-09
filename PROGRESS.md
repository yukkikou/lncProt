| step | status | time | key numbers |
|---|---|---|---|
| 02_parse_gencode | RUNNING | 2026-09-30 15:14:04 |  |
| 02_parse_gencode | DONE(42s) | 2026-09-30 15:14:46 |  | genes=63086 transcripts=254070 exons=1668627; groups: pcg=20065, lnc=19258, pseudo=15210, small=7082, other=1471; summary -> /media/share/workdir/1610305236/lncProt/step_summaries/02_parse_gencode.md |
| 01_download | FAILED(code=1) | 2026-09-30 16:30:15 | line=01_download cmd=return |
| 03a_cpg_validation | RUNNING | 2026-09-30 16:49:19 |  |
| 03a_cpg_validation | DONE(0s) | 2026-09-30 16:49:19 |  | lines=30344 mapped=6771 rate=0.223 gap_spanning=48 |
| 03b_chromhmm | FAILED(code=255) | 2026-09-30 16:49:20 | line=03b_chromhmm cmd=N_HIT=$(bt |
| test_chainlift | RUNNING | 2026-09-30 16:49:53 |  |
| test_chainlift | DONE(0s) | 2026-09-30 16:49:53 |  | lines=3 mapped=0 rate=0.000 gap_spanning=0 |
| test_chainlift2 | RUNNING | 2026-09-30 16:53:08 |  |
| test_chainlift2 | DONE(0s) | 2026-09-30 16:53:08 |  | lines=3 mapped=3 rate=1.000 gap_spanning=0 |
| 03a_cpg_validation | RUNNING | 2026-09-30 16:55:13 |  |
| 03a_cpg_validation | DONE(0s) | 2026-09-30 16:55:13 |  | lines=30344 mapped=26574 rate=0.876 gap_spanning=148 |
| 03b_chromhmm | FAILED(code=255) | 2026-09-30 16:55:14 | line=03b_chromhmm cmd=N_HIT=$(bt |
| 04_classify_loci | RUNNING | 2026-09-30 16:55:15 |  |
| 04_classify_loci | FAILED | 2026-09-30 16:55:22 |   File "/media/share/workdir/1610305236/lncProt/scripts/04_classify_loci.py", line 430, in main     tssmap = major_tss(genes, tx, cage)    File "/media/share/workdir/1610305236/lncProt/scripts/04_classify_loci.py", line 128, in major_tss     write_bed(cand, os.path.join(DATA, "tmp_cand_tss.bed"),    |
| 03a_cpg_validation | RUNNING | 2026-09-30 16:58:09 |  |
| 03a_cpg_validation | DONE(0s) | 2026-09-30 16:58:09 |  | lines=30344 mapped=26574 rate=0.876 gap_spanning=148 |
| 03a_lift_E003 | RUNNING | 2026-09-30 16:58:11 |  |
| 03a_lift_E003 | DONE(3s) | 2026-09-30 16:58:14 |  | lines=451435 mapped=439096 rate=0.973 gap_spanning=3666 |
| 03a_lift_E017 | RUNNING | 2026-09-30 16:58:17 |  |
| 03a_lift_E017 | DONE(4s) | 2026-09-30 16:58:21 |  | lines=635733 mapped=619920 rate=0.975 gap_spanning=4386 |
| 03a_lift_E067 | RUNNING | 2026-09-30 16:58:24 |  |
| 03a_lift_E067 | DONE(2s) | 2026-09-30 16:58:26 |  | lines=371167 mapped=360653 rate=0.972 gap_spanning=3424 |
| 03a_lift_E071 | RUNNING | 2026-09-30 16:58:28 |  |
| 03a_lift_E071 | DONE(3s) | 2026-09-30 16:58:32 |  | lines=475489 mapped=463123 rate=0.974 gap_spanning=3616 |
| 03a_lift_E088 | RUNNING | 2026-09-30 16:58:34 |  |
| 03a_lift_E088 | DONE(3s) | 2026-09-30 16:58:37 |  | lines=499248 mapped=486752 rate=0.975 gap_spanning=3979 |
| 03b_chromhmm_qc | RUNNING | 2026-09-30 16:58:39 |  |
| 03b_chromhmm_qc | DONE(4s) | 2026-09-30 16:58:43 |  | E003: TssActive/TssFlnk at PCG TSS=0.00; E017: TssActive/TssFlnk at PCG TSS=0.00; E067: TssActive/TssFlnk at PCG TSS=0.00; E071: TssActive/TssFlnk at PCG TSS=0.00; E088: TssActive/TssFlnk at PCG TSS=0.00; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/03b_chromhmm.md |
| 03b_chromhmm | DONE | 2026-09-30 16:58:43 | cells=5 |
| 04_classify_loci | RUNNING | 2026-09-30 16:58:44 |  |
| 04_classify_loci | FAILED | 2026-09-30 16:59:02 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 77, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/04_classify_loci.py", line 435, in main     meta = classify(genes, tx, ex, tssmap, st)    File "/medi |
| 03a_cpg_validation | RUNNING | 2026-09-30 17:03:52 |  |
| 03a_cpg_validation | DONE(0s) | 2026-09-30 17:03:52 |  | lines=30344 mapped=26574 rate=0.876 gap_spanning=148 |
| 03a_lift_E003 | RUNNING | 2026-09-30 17:03:54 |  |
| 03a_lift_E003 | DONE(3s) | 2026-09-30 17:03:57 |  | lines=451435 mapped=439096 rate=0.973 gap_spanning=3666 |
| 03a_lift_E017 | RUNNING | 2026-09-30 17:04:00 |  |
| 03a_lift_E017 | DONE(4s) | 2026-09-30 17:04:04 |  | lines=635733 mapped=619920 rate=0.975 gap_spanning=4386 |
| 03a_lift_E067 | RUNNING | 2026-09-30 17:04:07 |  |
| 03a_lift_E067 | DONE(2s) | 2026-09-30 17:04:09 |  | lines=371167 mapped=360653 rate=0.972 gap_spanning=3424 |
| 03a_lift_E071 | RUNNING | 2026-09-30 17:04:12 |  |
| 03a_lift_E071 | DONE(3s) | 2026-09-30 17:04:15 |  | lines=475489 mapped=463123 rate=0.974 gap_spanning=3616 |
| 03a_lift_E088 | RUNNING | 2026-09-30 17:04:17 |  |
| 03a_lift_E088 | DONE(3s) | 2026-09-30 17:04:21 |  | lines=499248 mapped=486752 rate=0.975 gap_spanning=3979 |
| 03b_chromhmm_qc | RUNNING | 2026-09-30 17:04:22 |  |
| 03b_chromhmm_qc | DONE(4s) | 2026-09-30 17:04:26 |  | E003: TssActive/TssFlnk at PCG TSS=0.39; E017: TssActive/TssFlnk at PCG TSS=0.36; E067: TssActive/TssFlnk at PCG TSS=0.42; E071: TssActive/TssFlnk at PCG TSS=0.40; E088: TssActive/TssFlnk at PCG TSS=0.32; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/03b_chromhmm.md |
| 03b_chromhmm | DONE | 2026-09-30 17:04:26 | cells=5 |
| 04_classify_loci | RUNNING | 2026-09-30 17:04:27 |  |
| 04_classify_loci | DONE(44s) | 2026-09-30 17:05:11 |  | genes=63086; cage_peaks=209911; cage_refined=27479 (43.6%); exon_overlap_pairs=5521; contained_in_pcg=5523; enhancer_overlap=1850 promoters; lnc classes: intergenic=10035, intronic=3403, antisense=2778, divergent=2692, sense_overlapping=350; PCG cage_refined=91.0%; promoter sets written (feat + fimo); scaffold-missing-seq genes=0; QC GAS5: group=lnc class=divergent refined=1; QC H19: group=lnc class=intergenic refined=1; QC SNHG1: group=lnc class=intergenic refined=1; QC NEAT1: group=lnc class=intergenic refined=1; QC MALAT1: group=lnc class=intergenic refined=1; QC LINC-PINT: group=lnc class=antisense refined=1; QC TSIX: group=lnc class=intergenic refined=0; QC XIST: group=lnc class=intergenic refined=1; QC PEG3: group=pcg class=PCG refined=1; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/04_classify_loci.md |
| 10_consv_chromatin | RUNNING | 2026-09-30 17:07:40 |  |
| 05_seq_features | RUNNING | 2026-09-30 17:07:40 |  |
| 07a_fimo_prepare | DONE | 2026-09-30 17:07:40 |  |
| 10_consv_chromatin | FAILED | 2026-09-30 17:07:41 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 77, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/10_consv_chromatin.py", line 49, in main     sing(SIF_PBW, "python",    File "/media/share/workdir/16 |
| 09_cage_activity | RUNNING | 2026-09-30 17:07:42 |  |
| 09_cage_activity | FAILED | 2026-09-30 17:07:43 |  During handling of the above exception, another exception occurred:   Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 77, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/09_cage_activity.py", line 182 |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:07:50 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:07:50 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:07:50 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:07:50 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:07:51 | line=07a_fimo_scan cmd=return |
| 05_seq_features | FAILED | 2026-09-30 17:07:52 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 77, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/05_seq_features.py", line 75, in main     cl = bt("closest", "-a", os.path.join(RES, "promoters_feat. |
| 06_te_promoters | RUNNING | 2026-09-30 17:07:53 |  |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:08:01 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:08:01 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:08:01 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:08:01 | line=07a_fimo_scan cmd=return |
| 06_te_promoters | FAILED | 2026-09-30 17:08:01 |   File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packages/pandas/io/parsers/readers.py", line 1026, in read_csv     return _read(filepath_or_buffer, kwds)    File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packages/pandas |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:08:01 | line=07a_fimo_scan cmd=return |
| 07b_remap | FAILED(code=137) | 2026-09-30 17:10:06 | line=07b_remap cmd=awk |
| 07a_fimo_prepare | DONE | 2026-09-30 17:13:43 |  |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:13:50 | line=07a_fimo_scan cmd=return |
| 10_consv_chromatin | RUNNING | 2026-09-30 17:23:43 |  |
| 05_seq_features | RUNNING | 2026-09-30 17:23:43 |  |
| 10_consv_chromatin | FAILED | 2026-09-30 17:23:44 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/10_consv_chromatin.py", line 115, in main     st.note("chromhmm by class: " + "; ".join(    File "/me |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:23:52 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:23:52 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:23:52 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:23:52 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:23:52 | line=07a_fimo_scan cmd=return |
| 05_seq_features | FAILED | 2026-09-30 17:23:55 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/05_seq_features.py", line 75, in main     cl = bt("closest", "-a", os.path.join(RES, "promoters_feat. |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:24:02 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:24:02 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:24:02 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:24:02 | line=07a_fimo_scan cmd=return |
| 07a_fimo | FAILED(code=1) | 2026-09-30 17:24:03 | line=07a_fimo_scan cmd=return |
| 07a_fimo_prepare | DONE | 2026-09-30 17:27:00 |  |
| 10_consv_chromatin | RUNNING | 2026-09-30 17:28:06 |  |
| 05_seq_features | RUNNING | 2026-09-30 17:28:06 |  |
| 10_consv_chromatin | DONE(2s) | 2026-09-30 17:28:08 |  | phastCons mean (PCG)=nan; chromhmm cells available=['E003', 'E017', 'E067', 'E071', 'E088']; E003: assigned=0 TssA=0; E017: assigned=0 TssA=0; E067: assigned=0 TssA=0; E071: assigned=0 TssA=0; E088: assigned=0 TssA=0; chromhmm by class: PCG: prom=nan; divergent: prom=nan; antisense: prom=nan; intronic: prom=nan; sense_overlapping: prom=nan; intergenic: prom=nan; pseudogene: prom=nan; transcribed_pseudogene: prom=nan; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/10_consv_chromatin.md |
| 05_seq_features | DONE(14s) | 2026-09-30 17:28:20 |  | promoters scored=54533; class medians: PCG: GC=0.60 o/e=0.65 TATA=0.01 HCP=0.60; divergent: GC=0.64 o/e=0.77 TATA=0.00 HCP=0.78; antisense: GC=0.50 o/e=0.35 TATA=0.00 HCP=0.13; intronic: GC=0.46 o/e=0.28 TATA=0.00 HCP=0.05; sense_overlapping: GC=0.53 o/e=0.44 TATA=0.01 HCP=0.27; intergenic: GC=0.48 o/e=0.31 TATA=0.00 HCP=0.10; pseudogene: GC=0.44 o/e=0.28 TATA=0.00 HCP=0.03; transcribed_pseudogene: GC=0.49 o/e=0.35 TATA=0.00 HCP=0.14; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/05_seq_features.md |
| 06_te_promoters | RUNNING | 2026-09-30 17:28:21 |  |
| 06_te_promoters | FAILED | 2026-09-30 17:28:29 |   File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packages/pandas/io/parsers/readers.py", line 1026, in read_csv     return _read(filepath_or_buffer, kwds)    File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packages/pandas |
| 10_consv_chromatin | RUNNING | 2026-09-30 17:30:08 |  |
| 04b_antisense_lnc_flag | RUNNING | 2026-09-30 17:30:17 |  |
| 04b_antisense_lnc_flag | DONE(7s) | 2026-09-30 17:30:24 |  | antisense_to_lnc flagged=1772; TSIX flag=[1]; master table updated atomically |
| 10_consv_chromatin | DONE(30s) | 2026-09-30 17:30:38 |  | phastCons mean (PCG)=0.242; chromhmm cells available=['E003', 'E017', 'E067', 'E071', 'E088']; E003: assigned=51828 TssA=13561; E017: assigned=51829 TssA=13367; E067: assigned=51686 TssA=16045; E071: assigned=51694 TssA=15079; E088: assigned=51712 TssA=12939; chromhmm by class: PCG: prom=0.55; divergent: prom=0.74; antisense: prom=0.10; intronic: prom=0.05; sense_overlapping: prom=0.27; intergenic: prom=0.06; pseudogene: prom=0.01; transcribed_pseudogene: prom=0.07; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/10_consv_chromatin.md |
| 05_seq_features | RUNNING | 2026-09-30 17:33:34 |  |
| 05_seq_features | DONE(13s) | 2026-09-30 17:33:48 |  | promoters scored=54533; class medians: PCG: GC=0.60 o/e=0.65 TATA=0.06 HCP=0.60; divergent: GC=0.64 o/e=0.77 TATA=0.02 HCP=0.78; antisense: GC=0.50 o/e=0.35 TATA=0.03 HCP=0.13; intronic: GC=0.46 o/e=0.28 TATA=0.02 HCP=0.05; sense_overlapping: GC=0.53 o/e=0.44 TATA=0.03 HCP=0.27; intergenic: GC=0.48 o/e=0.31 TATA=0.03 HCP=0.10; pseudogene: GC=0.44 o/e=0.28 TATA=0.02 HCP=0.03; transcribed_pseudogene: GC=0.49 o/e=0.35 TATA=0.03 HCP=0.14; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/05_seq_features.md |
| 06_te_promoters | RUNNING | 2026-09-30 17:33:48 |  |
| 06_te_promoters | FAILED | 2026-09-30 17:33:56 |   File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packages/pandas/io/parsers/readers.py", line 1026, in read_csv     return _read(filepath_or_buffer, kwds)    File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packages/pandas |
| 06_te_promoters | RUNNING | 2026-09-30 17:34:54 |  |
| 06_te_promoters | DONE(136s) | 2026-09-30 17:37:10 |  | TE overlap rows=54381; TSS inside TE (unique genes)=20481; top TSS-donor TEs: MIRb=778, AluSx1=644, L2=624, AluSx=608, AluSz=582, MIR=549, AluJb=542, L2b_3end=491; TE by class: PCG: cov=0.090 tssTE=0.263 LTR=0.030; divergent: cov=0.048 tssTE=0.069 LTR=0.012; antisense: cov=0.288 tssTE=0.333 LTR=0.102; intronic: cov=0.399 tssTE=0.399 LTR=0.124; sense_overlapping: cov=0.232 tssTE=0.337 LTR=0.117; intergenic: cov=0.406 tssTE=0.433 LTR=0.219; pseudogene: cov=0.364 tssTE=0.559 LTR=0.121; transcribed_pseudogene: cov=0.226 tssTE=0.421 LTR=0.081; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/06_te_promoters.md |
| 07b_remap | DONE | 2026-09-30 17:37:52 | rows=38359161 |
| 11_report_figures | RUNNING | 2026-09-30 17:45:47 |  |
| 11_report_figures | DONE(5s) | 2026-09-30 17:45:52 |  | figures written; REPORT.md written |
| 11_report_figures | RUNNING | 2026-09-30 17:48:13 |  |
| 11_report_figures | DONE(5s) | 2026-09-30 17:48:18 |  | figures written; REPORT.md written |
| 07a_fimo_scan_3 | DONE | 2026-09-30 17:48:33 |  |
| 07a_fimo_scan_4 | DONE | 2026-09-30 17:48:36 |  |
| 07a_fimo_scan_1 | DONE | 2026-09-30 17:48:40 |  |
| 07a_fimo_scan_2 | DONE | 2026-09-30 17:48:41 |  |
| 07a_fimo_scan_0 | DONE | 2026-09-30 17:48:41 |  |
| 07a_fimo_scan_6 | DONE | 2026-09-30 18:08:56 |  |
| 07a_fimo_scan_5 | DONE | 2026-09-30 18:08:59 |  |
| 07a_fimo_scan_8 | DONE | 2026-09-30 18:09:06 |  |
| 07a_fimo_scan_7 | DONE | 2026-09-30 18:09:20 |  |
| 07a_fimo_scan_9 | DONE | 2026-09-30 18:09:35 |  |
| 08_enrichment | RUNNING | 2026-09-30 18:09:36 |  |
| 08_enrichment | FAILED | 2026-09-30 18:09:39 |   File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packages/pandas/io/parsers/readers.py", line 1923, in read     ) = self._engine.read(  # type: ignore[attr-defined]    File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packa |
| 08_enrichment | RUNNING | 2026-09-30 18:11:31 |  |
| 08_enrichment | FAILED | 2026-09-30 18:11:34 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/08_enrichment.py", line 102, in main     hits = load_fimo()    File "/media/share/workdir/1610305236/ |
| 08_enrichment | RUNNING | 2026-09-30 18:15:21 |  |
| 08_enrichment | FAILED | 2026-09-30 18:15:46 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/08_enrichment.py", line 102, in main     hits = load_fimo()    File "/media/share/workdir/1610305236/ |
| 08_enrichment | RUNNING | 2026-09-30 18:18:39 |  |
| 08_enrichment | DONE(401s) | 2026-09-30 18:25:20 |  | fimo hits=18,090,665 motifs=1176; presence matrix: 54533 genes x 1131 TFs; motif_enrichment rows=7917, sig(FDR<0.05)=6385; remap hits=38,359,161 tfs=1242; binding rows=8628 sig=7450; candidate TFs (rows)=2257; dual-evidence=91; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/08_enrichment.md |
| 11_report_figures | RUNNING | 2026-09-30 18:27:06 |  |
| 11_report_figures | DONE(6s) | 2026-09-30 18:27:12 |  | figures written; REPORT.md written |
| 09_cage_activity | RUNNING | 2026-10-01 00:12:51 |  |
| 09_cage_activity | FAILED | 2026-10-01 00:15:55 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/09_cage_activity.py", line 86, in main     act["active"] = np.where(ok, active[idx[ok].astype(int).va |
| 09_cage_activity | RUNNING | 2026-10-01 00:18:07 |  |
| 09_cage_activity | DONE(180s) | 2026-10-01 00:21:07 |  | loading TPM matrix (large, may take a while); matrix: 209913 peaks x 1829 samples; promoters linked to TPM rows: 26447 (48%); divergent pairs used=1656 median rho=0.36; control median rho=0.35; activity by class: PCG: active=0.91; divergent: active=0.62; antisense: active=0.28; intronic: active=0.20; sense_overlapping: active=0.58; intergenic: active=0.33; pseudogene: active=0.09; transcribed_pseudogene: active=0.17; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/09_cage_activity.md |
| 11_report_figures | RUNNING | 2026-10-01 00:25:28 |  |
| 11_report_figures | DONE(7s) | 2026-10-01 00:25:35 |  | figures written; REPORT.md written |
| 12a_pub_stats | RUNNING | 2026-10-01 10:06:50 |  |
| 12a_pub_stats | FAILED | 2026-10-01 10:06:52 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/12a_pub_stats.py", line 126, in main     + df_to_md(sig, ".3f", max_rows=60)    File "/media/share/wo |
| 12a_pub_stats | RUNNING | 2026-10-01 10:07:06 |  |
| 12a_pub_stats | DONE(2s) | 2026-10-01 10:07:09 |  | tests=98; significant(FDR<0.05)=94; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/12a_pub_stats.md; divergent: 13 significant tests of 14; intergenic: 14 significant tests of 14 |
| 12b_te_subfamily_age | RUNNING | 2026-10-01 10:07:57 |  |
| 12b_te_subfamily_age | FAILED | 2026-10-01 10:07:57 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/12b_te_subfamily_age.py", line 54, in main     tt["cls"] = tt["cls"].fillna(    File "/media/share/ho |
| 12b_te_subfamily_age | RUNNING | 2026-10-01 10:08:09 |  |
| 12b_te_subfamily_age | DONE(0s) | 2026-10-01 10:08:10 |  | age mix: divergent: young=0.06 ancient=0.45; antisense: young=0.10 ancient=0.20; intronic: young=0.07 ancient=0.25; sense_overlapping: young=0.07 ancient=0.25; intergenic: young=0.10 ancient=0.21; pseudogene: young=0.09 ancient=0.12; transcribed_pseudogene: young=0.10 ancient=0.17; PCG: young=0.07 ancient=0.32; top intergenic donors: MIRb(144), L2(138), MIR(115), AluSx1(98), AluJb(89), AluSx(84), AluSz(84), L2b_3end(84); summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/12b_te_subfamily_age.md |
| 12d_robustness | RUNNING | 2026-10-01 10:08:40 |  |
| 12d_robustness | DONE(1s) | 2026-10-01 10:08:41 |  | divergent sensitivity: {'500bp': 2118, '1000bp (primary)': 2692, '2000bp': 3201}; active-only rows=12 (n active PCG=18260); GC-controlled conservation: divergent: 0.046->-0.003; antisense: -0.140->-0.033; intronic: -0.283->-0.124; sense_overlapping: -0.056->0.011; intergenic: -0.424->-0.283; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/12d_robustness.md |
| 12e_partner_pcgs | RUNNING | 2026-10-01 10:12:00 |  |
| 12e_partner_pcgs | DONE(5s) | 2026-10-01 10:12:04 |  | unique partner PCGs=2642; partner traits: gc: 0.643 vs 0.593; cpg_oe: 0.796 vs 0.625; phast_mean: 0.271 vs 0.238; active: 0.980 vs 0.899; breadth: 1376.936 vs 1059.846; top GO: protein binding (p=2.1e-58); regulation of primary metabolic process (p=5.9e-45); regulation of biological process (p=3.2e-34); regulation of nucleobase-containing compound metabolic process (p=9.4e-34); biological regulation (p=9.5e-34); summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/12e_partner_pcgs.md |
| 12c_gc_matched_motifs | RUNNING | 2026-10-01 10:13:22 |  |
| 12c_gc_matched_motifs | FAILED | 2026-10-01 10:13:26 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/12c_gc_matched_motifs.py", line 129, in main     raise RuntimeError(r.stderr[:800])  RuntimeError: Er |
| 12c_gc_matched_motifs | RUNNING | 2026-10-01 10:13:45 |  |
| 12c_gc_matched_motifs | FAILED | 2026-10-01 10:13:46 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/12c_gc_matched_motifs.py", line 84, in main     name = b[0].split()[1]  AttributeError: 'list' object |
| 12c_gc_matched_motifs | RUNNING | 2026-10-01 10:14:33 |  |
| 12c_gc_matched_motifs | FAILED | 2026-10-01 10:14:33 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/12c_gc_matched_motifs.py", line 84, in main     name = b[0].strip().split()[1]  AttributeError: 'list |
| 12c_gc_matched_motifs | RUNNING | 2026-10-01 10:16:28 |  |
| 12c_gc_matched_motifs | FAILED | 2026-10-01 10:16:32 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/12c_gc_matched_motifs.py", line 132, in main     raise RuntimeError(r.stderr[:800])  RuntimeError: Er |
| 12c_gc_matched_motifs | RUNNING | 2026-10-01 10:20:13 |  |
| 12c_gc_matched_motifs | FAILED | 2026-10-01 10:20:17 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/12c_gc_matched_motifs.py", line 133, in main     raise RuntimeError(r.stderr[:800])  RuntimeError: Er |
| 13a_toil | DONE | 2026-10-01 10:22:26 | bytes=1323254426 |
| 13b_parse_toil | RUNNING | 2026-10-01 10:22:27 |  |
| 13b_parse_toil | FAILED | 2026-10-01 10:26:40 |   File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packages/pandas/io/parsers/readers.py", line 626, in _read     return parser.read(nrows)    File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packages/pandas/io/parsers/reade |
| 12c_gc_matched_motifs | RUNNING | 2026-10-01 10:34:35 |  |
| 12c_gc_matched_motifs | FAILED | 2026-10-01 10:34:40 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/12c_gc_matched_motifs.py", line 133, in main     raise RuntimeError(r.stderr[:800])  RuntimeError: Er |
| 12c_gc_matched_motifs | RUNNING | 2026-10-01 10:35:11 |  |
| 12c_gc_matched_motifs | DONE(42s) | 2026-10-01 10:35:54 |  | classes=['divergent', 'intergenic'] n_top_TFs=30; motifs found for 4/30 TFs in HOCOMOCO; divergent observed: 1/2692 promoters with hits; divergent gc_matched: 1/2692 promoters with hits; intergenic observed: 1/10035 promoters with hits; intergenic gc_matched: 1/10035 promoters with hits; ratios: divergent: 1.00; intergenic: 1.00; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/12c_gc_matched_motifs.md |
| 12c_gc_matched_motifs | RUNNING | 2026-10-01 10:36:14 |  |
| 12c_gc_matched_motifs | DONE(41s) | 2026-10-01 10:36:56 |  | classes=['divergent', 'intergenic'] n_top_TFs=30; motifs found for 4/30 TFs in HOCOMOCO; divergent observed: 2643/2692 promoters with hits; divergent gc_matched: 2486/2692 promoters with hits; intergenic observed: 8108/10035 promoters with hits; intergenic gc_matched: 6278/10035 promoters with hits; ratios: divergent: 1.06; intergenic: 1.29; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/12c_gc_matched_motifs.md |
| 13b_parse_toil | RUNNING | 2026-10-01 10:37:40 |  |
| 13b_parse_toil | FAILED | 2026-10-01 10:41:54 |   File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/13b_parse_toil.py", line 88, in main     cw = ss.merge(stage, on="file_id", how="outer", suffixes=("", "_stage"))    File "/media/share/h |
| 13b_parse_toil | RUNNING | 2026-10-01 10:44:51 |  |
| 13b_parse_toil | DONE(698s) | 2026-10-01 10:56:29 |  | matrix: 48275 class-genes x 19131 samples; crosswalk rows=22550; with stage=11275; summary -> /media/share/workdir/1610305236/lncProt/results/tissue_disease/step_summaries/13b_sample_counts.md; accounting rows=5; GTEx tissues(n>=20)=0; |
| 13c_gtex_tissue | RUNNING | 2026-10-01 10:56:32 |  |
| 13c_gtex_tissue | FAILED | 2026-10-01 10:56:47 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/13c_gtex_tissue.py", line 75, in main     mean_by_t = np.vstack([d[:, c].mean(axis=1) for t, c in gco |
| 13c_gtex_tissue | RUNNING | 2026-10-01 11:01:10 |  |
| 13c_gtex_tissue | DONE(252s) | 2026-10-01 11:05:23 |  | tissues used (n>=20): 49; tissue rows=49; tau: divergent=nan; antisense=nan; intronic=nan; sense_overlapping=nan; intergenic=nan; PCG=nan; summary -> /media/share/workdir/1610305236/lncProt/results/tissue_disease/step_summaries/13c_gtex_tissue.md |
| 13d_merge | RUNNING | 2026-10-01 11:05:24 |  |
| 13d_merge | FAILED | 2026-10-01 11:05:24 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/13d_tcga_disease.py", line 127, in merge     shards = sorted(f for f in os.listdir(SHARDS) if f.start |
| 13d_TCGA-ESCA | RUNNING | 2026-10-01 11:05:25 |  |
| 13d_TCGA-DLBC | RUNNING | 2026-10-01 11:05:25 |  |
| 13d_TCGA-COAD | RUNNING | 2026-10-01 11:05:25 |  |
| 13d_TCGA-CESC | RUNNING | 2026-10-01 11:05:25 |  |
| 13d_TCGA-CHOL | RUNNING | 2026-10-01 11:05:25 |  |
| 13d_TCGA-BRCA | RUNNING | 2026-10-01 11:05:25 |  |
| 13d_TCGA-ACC | RUNNING | 2026-10-01 11:05:25 |  |
| 13d_TCGA-BLCA | RUNNING | 2026-10-01 11:05:25 |  |
| 13d_TCGA-DLBC | DONE(15s) | 2026-10-01 11:05:39 |  | tumor=47 normal=0; classes=6 |
| 13d_TCGA-CESC | DONE(15s) | 2026-10-01 11:05:40 |  | tumor=306 normal=3; classes=6 |
| 13d_TCGA-ACC | DONE(16s) | 2026-10-01 11:05:40 |  | tumor=77 normal=0; classes=6 |
| 13d_TCGA-GBM | RUNNING | 2026-10-01 11:05:41 |  |
| 13d_TCGA-HNSC | RUNNING | 2026-10-01 11:05:41 |  |
| 13d_TCGA-KICH | RUNNING | 2026-10-01 11:05:42 |  |
| 13d_TCGA-CHOL | DONE(46s) | 2026-10-01 11:06:10 |  | tumor=35 normal=9; classes=6 |
| 13d_TCGA-ESCA | DONE(47s) | 2026-10-01 11:06:11 |  | tumor=182 normal=13; classes=6 |
| 13d_TCGA-KIRC | RUNNING | 2026-10-01 11:06:12 |  |
| 13d_TCGA-COAD | DONE(47s) | 2026-10-01 11:06:12 |  | tumor=289 normal=41; classes=6 |
| 13d_TCGA-BLCA | DONE(48s) | 2026-10-01 11:06:12 |  | tumor=405 normal=19; classes=6 |
| 13d_TCGA-KIRP | RUNNING | 2026-10-01 11:06:13 |  |
| 13d_TCGA-LGG | RUNNING | 2026-10-01 11:06:13 |  |
| 13d_TCGA-LIHC | RUNNING | 2026-10-01 11:06:14 |  |
| 13d_TCGA-BRCA | DONE(50s) | 2026-10-01 11:06:15 |  | tumor=1099 normal=113; classes=6 |
| 13d_TCGA-LUAD | RUNNING | 2026-10-01 11:06:16 |  |
| 13d_TCGA-HNSC | DONE(46s) | 2026-10-01 11:06:27 |  | tumor=520 normal=44; classes=6 |
| 13d_TCGA-GBM | DONE(47s) | 2026-10-01 11:06:28 |  | tumor=153 normal=5; classes=6 |
| 13d_TCGA-LGG | DONE(15s) | 2026-10-01 11:06:28 |  | tumor=509 normal=0; classes=6 |
| 13d_TCGA-LUSC | RUNNING | 2026-10-01 11:06:29 |  |
| 13d_TCGA-MESO | RUNNING | 2026-10-01 11:06:29 |  |
| 13d_TCGA-OV | RUNNING | 2026-10-01 11:06:29 |  |
| 13d_TCGA-KICH | DONE(48s) | 2026-10-01 11:06:29 |  | tumor=66 normal=25; classes=6 |
| 13d_TCGA-PAAD | RUNNING | 2026-10-01 11:06:31 |  |
| 13d_TCGA-MESO | DONE(15s) | 2026-10-01 11:06:44 |  | tumor=87 normal=0; classes=6 |
| 13d_TCGA-OV | DONE(15s) | 2026-10-01 11:06:44 |  | tumor=419 normal=0; classes=6 |
| 13d_TCGA-PCPG | RUNNING | 2026-10-01 11:06:45 |  |
| 13d_TCGA-PRAD | RUNNING | 2026-10-01 11:06:45 |  |
| 13d_TCGA-PAAD | DONE(15s) | 2026-10-01 11:06:46 |  | tumor=179 normal=4; classes=6 |
| 13d_TCGA-READ | RUNNING | 2026-10-01 11:06:47 |  |
| 13d_TCGA-KIRC | DONE(46s) | 2026-10-01 11:06:58 |  | tumor=530 normal=72; classes=6 |
| 13d_TCGA-KIRP | DONE(46s) | 2026-10-01 11:06:59 |  | tumor=288 normal=32; classes=6 |
| 13d_TCGA-SARC | RUNNING | 2026-10-01 11:06:59 |  |
| 13d_TCGA-LIHC | DONE(46s) | 2026-10-01 11:06:59 |  | tumor=369 normal=50; classes=6 |
| 13d_TCGA-PCPG | DONE(15s) | 2026-10-01 11:07:00 |  | tumor=179 normal=3; classes=6 |
| 13d_TCGA-SKCM | RUNNING | 2026-10-01 11:07:00 |  |
| 13d_TCGA-STAD | RUNNING | 2026-10-01 11:07:01 |  |
| 13d_TCGA-TGCT | RUNNING | 2026-10-01 11:07:01 |  |
| 13d_TCGA-LUAD | DONE(48s) | 2026-10-01 11:07:04 |  | tumor=513 normal=59; classes=6 |
| 13d_TCGA-THCA | RUNNING | 2026-10-01 11:07:06 |  |
| 13d_TCGA-SARC | DONE(15s) | 2026-10-01 11:07:14 |  | tumor=259 normal=2; classes=6 |
| 13d_TCGA-SKCM | DONE(15s) | 2026-10-01 11:07:15 |  | tumor=468 normal=1; classes=6 |
| 13d_TCGA-THYM | RUNNING | 2026-10-01 11:07:16 |  |
| 13d_TCGA-TGCT | DONE(15s) | 2026-10-01 11:07:16 |  | tumor=148 normal=0; classes=6 |
| 13d_TCGA-LUSC | DONE(48s) | 2026-10-01 11:07:16 |  | tumor=498 normal=50; classes=6 |
| 13d_TCGA-UCEC | RUNNING | 2026-10-01 11:07:17 |  |
| 13d_TCGA-UCS | RUNNING | 2026-10-01 11:07:17 |  |
| 13d_TCGA-UVM | RUNNING | 2026-10-01 11:07:18 |  |
| 13d_TCGA-THYM | DONE(15s) | 2026-10-01 11:07:30 |  | tumor=119 normal=2; classes=6 |
| 13d_TCGA-UCS | DONE(15s) | 2026-10-01 11:07:32 |  | tumor=57 normal=0; classes=6 |
| 13d_TCGA-PRAD | DONE(47s) | 2026-10-01 11:07:32 |  | tumor=496 normal=52; classes=6 |
| 13d_TCGA-UVM | DONE(15s) | 2026-10-01 11:07:32 |  | tumor=79 normal=0; classes=6 |
| 13d_TCGA-READ | DONE(46s) | 2026-10-01 11:07:33 |  | tumor=92 normal=10; classes=6 |
| 13d_TCGA-STAD | DONE(48s) | 2026-10-01 11:07:48 |  | tumor=410 normal=35; classes=6 |
| 13d_TCGA-THCA | DONE(48s) | 2026-10-01 11:07:54 |  | tumor=512 normal=59; classes=6 |
| 13d_TCGA-UCEC | DONE(46s) | 2026-10-01 11:08:02 |  | tumor=180 normal=23; classes=6 |
| 13d_merge | RUNNING | 2026-10-01 11:11:42 |  |
| 13d_merge | DONE(0s) | 2026-10-01 11:11:42 |  | cancers=32 rows=192; top LTR derepression: TCGA-CHOL: d=0.30; TCGA-KIRP: d=0.12; TCGA-LIHC: d=0.12; TCGA-KICH: d=0.10; TCGA-BLCA: d=0.01; summary -> /media/share/workdir/1610305236/lncProt/results/tissue_disease/step_summaries/13d_tcga_disease.md |
| 13d_stage | RUNNING | 2026-10-01 11:11:43 |  |
| 13d_stage | DONE(15s) | 2026-10-01 11:11:58 |  | no cancer met stage thresholds (n>=10 early & advanced) |
| 13d_stage | RUNNING | 2026-10-01 11:13:28 |  |
| 13d_stage | DONE(25s) | 2026-10-01 11:13:53 |  | stage rows=75 across 15 cancers |
| 13e_integrate_td | RUNNING | 2026-10-01 11:13:54 |  |
| 13e_integrate_td | DONE(1s) | 2026-10-01 11:13:55 |  | summary: divergent: share=0.25; antisense: share=0.25; intronic: share=0.16; sense_overlapping: share=0.03; intergenic: share=0.31; PCG: share=nan; FANTOM tau vs GTEx tau: rho=nan (n=48275); summary -> /media/share/workdir/1610305236/lncProt/results/tissue_disease/step_summaries/13e_integration.md |
| 14_td_figures | RUNNING | 2026-10-01 11:13:56 |  |
| 14_td_figures | DONE(4s) | 2026-10-01 11:13:59 |  | tissue/disease figures done |
| 14_td_figures | RUNNING | 2026-10-01 11:14:37 |  |
| 14_td_figures | DONE(4s) | 2026-10-01 11:14:41 |  | tissue/disease figures done |
| 13d_TCGA-CHOL | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-BLCA | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-CESC | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-ACC | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-BRCA | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-COAD | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-GBM | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-HNSC | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-DLBC | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-ESCA | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-KICH | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-KIRC | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-UVM | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-UCS | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-TGCT | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-UCEC | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-READ | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-THCA | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-SARC | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-THYM | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-SKCM | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-PRAD | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-STAD | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-KIRP | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-LGG | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-PAAD | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-LIHC | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-MESO | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-OV | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-LUAD | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-PCPG | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-LUSC | RUNNING | 2026-10-01 11:15:34 |  |
| 13d_TCGA-UVM | DONE(26s) | 2026-10-01 11:16:00 |  | tumor=79 normal=0; classes=6 |
| 13d_TCGA-UCS | DONE(26s) | 2026-10-01 11:16:00 |  | tumor=57 normal=0; classes=6 |
| 13d_TCGA-CESC | DONE(30s) | 2026-10-01 11:16:04 |  | tumor=306 normal=3; classes=6 |
| 13d_TCGA-ACC | DONE(31s) | 2026-10-01 11:16:05 |  | tumor=77 normal=0; classes=6 |
| 13d_TCGA-DLBC | DONE(32s) | 2026-10-01 11:16:05 |  | tumor=47 normal=0; classes=6 |
| 13d_TCGA-SKCM | DONE(35s) | 2026-10-01 11:16:09 |  | tumor=468 normal=1; classes=6 |
| 13d_TCGA-TGCT | DONE(35s) | 2026-10-01 11:16:09 |  | tumor=148 normal=0; classes=6 |
| 13d_TCGA-MESO | DONE(35s) | 2026-10-01 11:16:09 |  | tumor=87 normal=0; classes=6 |
| 13d_TCGA-PAAD | DONE(35s) | 2026-10-01 11:16:10 |  | tumor=179 normal=4; classes=6 |
| 13d_TCGA-OV | DONE(36s) | 2026-10-01 11:16:10 |  | tumor=419 normal=0; classes=6 |
| 13d_TCGA-PCPG | DONE(36s) | 2026-10-01 11:16:10 |  | tumor=179 normal=3; classes=6 |
| 13d_TCGA-LGG | DONE(36s) | 2026-10-01 11:16:10 |  | tumor=509 normal=0; classes=6 |
| 13d_TCGA-THYM | DONE(37s) | 2026-10-01 11:16:11 |  | tumor=119 normal=2; classes=6 |
| 13d_TCGA-SARC | DONE(37s) | 2026-10-01 11:16:11 |  | tumor=259 normal=2; classes=6 |
| 13d_TCGA-CHOL | DONE(57s) | 2026-10-01 11:16:31 |  | tumor=35 normal=9; classes=6 |
| 13d_TCGA-KICH | DONE(57s) | 2026-10-01 11:16:31 |  | tumor=66 normal=25; classes=6 |
| 13d_TCGA-GBM | DONE(63s) | 2026-10-01 11:16:37 |  | tumor=153 normal=5; classes=6 |
| 13d_TCGA-READ | DONE(65s) | 2026-10-01 11:16:39 |  | tumor=92 normal=10; classes=6 |
| 13d_TCGA-THCA | DONE(66s) | 2026-10-01 11:16:40 |  | tumor=512 normal=59; classes=6 |
| 13d_TCGA-HNSC | DONE(67s) | 2026-10-01 11:16:41 |  | tumor=520 normal=44; classes=6 |
| 13d_TCGA-COAD | DONE(67s) | 2026-10-01 11:16:41 |  | tumor=289 normal=41; classes=6 |
| 13d_TCGA-LIHC | DONE(67s) | 2026-10-01 11:16:41 |  | tumor=369 normal=50; classes=6 |
| 13d_TCGA-BLCA | DONE(68s) | 2026-10-01 11:16:42 |  | tumor=405 normal=19; classes=6 |
| 13d_TCGA-PRAD | DONE(68s) | 2026-10-01 11:16:43 |  | tumor=496 normal=52; classes=6 |
| 13d_TCGA-UCEC | DONE(69s) | 2026-10-01 11:16:43 |  | tumor=180 normal=23; classes=6 |
| 13d_TCGA-KIRP | DONE(69s) | 2026-10-01 11:16:44 |  | tumor=288 normal=32; classes=6 |
| 13d_TCGA-ESCA | DONE(71s) | 2026-10-01 11:16:44 |  | tumor=182 normal=13; classes=6 |
| 13d_TCGA-KIRC | DONE(71s) | 2026-10-01 11:16:45 |  | tumor=530 normal=72; classes=6 |
| 13d_TCGA-LUSC | DONE(73s) | 2026-10-01 11:16:47 |  | tumor=498 normal=50; classes=6 |
| 13d_TCGA-BRCA | DONE(75s) | 2026-10-01 11:16:48 |  | tumor=1099 normal=113; classes=6 |
| 13d_TCGA-LUAD | DONE(75s) | 2026-10-01 11:16:49 |  | tumor=513 normal=59; classes=6 |
| 13d_TCGA-STAD | DONE(80s) | 2026-10-01 11:16:54 |  | tumor=410 normal=35; classes=6 |
| 13d_merge | RUNNING | 2026-10-01 11:25:52 |  |
| 13d_merge | DONE(0s) | 2026-10-01 11:25:52 |  | cancers=32 rows=192; top LTR derepression: TCGA-KICH: d=0.19; TCGA-THCA: d=0.09; TCGA-GBM: d=0.04; TCGA-UCEC: d=0.00; TCGA-BLCA: d=0.00; summary -> /media/share/workdir/1610305236/lncProt/results/tissue_disease/step_summaries/13d_tcga_disease.md |
| 13d_stage | RUNNING | 2026-10-01 11:25:53 |  |
| 13d_stage | DONE(26s) | 2026-10-01 11:26:19 |  | stage rows=75 across 15 cancers |
| 13e_integrate_td | RUNNING | 2026-10-01 11:26:20 |  |
| 13e_integrate_td | DONE(1s) | 2026-10-01 11:26:21 |  | summary: divergent: share=0.25; antisense: share=0.25; intronic: share=0.16; sense_overlapping: share=0.03; intergenic: share=0.31; PCG: share=nan; FANTOM tau vs GTEx tau: rho=nan (n=48275); summary -> /media/share/workdir/1610305236/lncProt/results/tissue_disease/step_summaries/13e_integration.md |
| 14_td_figures | RUNNING | 2026-10-01 11:26:22 |  |
| 14_td_figures | DONE(4s) | 2026-10-01 11:26:26 |  | tissue/disease figures done |
| 13d_stage | RUNNING | 2026-10-01 11:26:42 |  |
| 13d_stage | DONE(26s) | 2026-10-01 11:27:07 |  | stage rows=75 across 15 cancers |
| 11_report_figures | RUNNING | 2026-10-01 11:27:38 |  |
| 11_report_figures | DONE(7s) | 2026-10-01 11:27:45 |  | figures written; REPORT.md written |
| 15_final_gaps | RUNNING | 2026-10-02 09:42:48 |  |
| 15_final_gaps | DONE(133s) | 2026-10-02 09:45:02 |  | (a) FANTOM~GTEx tau rho=nan p=nan n=23746; (b) testis LTR share 0.075 vs others 0.047; OR=1.64 p=1.6e-22; (c) pan-cancer LTR delta=0.005 95%CI[-0.025,0.036] P(<=0)=0.383; (d) DUX4: LTR-frac bound=0.371 bg=0.218 p=6.6e-09; PHOX2B: LTR-frac bound=0.105 bg=0.218 p=9.4e-12; BATF3: LTR-frac bound=0.151 bg=0.218 p=4.2e-02; CTCF: LTR-frac bound=0.141 bg=0.218 p=3.0e-111; MAZ: LTR-frac bound=0.090 bg=0.218 p=2.9e-137; (e) normal-activity vs dysregulation rho=0.70 p=0.19 (class level, n=5); summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/15_final_gaps.md |
| 16_submission_figures | RUNNING | 2026-10-02 09:47:56 |  |
| 16_submission_figures | FAILED | 2026-10-02 09:47:59 |   File "/media/share/workdir/1610305236/lncProt/scripts/16_submission_figures.py", line 431, in main     fig1(); fig2(); fig3(); fig3b(); fig4(); fig4b(); fig5(); fig6()    File "/media/share/workdir/1610305236/lncProt/scripts/16_submission_figures.py", line 200, in fig3b     old = [age.loc[c, "frac |
| 16_submission_figures | RUNNING | 2026-10-02 09:48:09 |  |
| 16_submission_figures | DONE(11s) | 2026-10-02 09:48:20 |  | submission figures -> /media/share/workdir/1610305236/lncProt/figures/sub |
| 19_ccle_cgga | RUNNING | 2026-10-02 10:41:38 |  |
| 19_ccle_cgga | FAILED | 2026-10-02 10:41:55 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/19_ccle_cgga.py", line 196, in main     msg1, ccle_res = ccle(sym2cls, sym2ltr)    File "/media/share |
| 19_ccle_cgga | RUNNING | 2026-10-02 10:43:28 |  |
| 19_ccle_cgga | FAILED | 2026-10-02 10:43:44 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/19_ccle_cgga.py", line 196, in main     msg1, ccle_res = ccle(sym2cls, sym2ltr)    File "/media/share |
| 19_ccle_cgga | RUNNING | 2026-10-02 10:45:35 |  |
| 19_ccle_cgga | FAILED | 2026-10-02 10:45:51 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/19_ccle_cgga.py", line 199, in main     msg1, ccle_res = ccle(sym2cls, sym2ltr)    File "/media/share |
| 19_ccle_cgga | RUNNING | 2026-10-02 10:47:56 |  |
| 19_ccle_cgga | DONE(74s) | 2026-10-02 10:49:11 |  | CCLE: 24 lineages, 1719 lines; DUX4~LTR-lnc rel-activation rho=-0.26 (p=5.6e-29, hi-lo=-0.05, MWU p=1.1e-23, n=1719); skip: CGGA matrices absent; summary -> /media/share/workdir/1610305236/lncProt/results/ccle_cgga/step_summaries/19_ccle_cgga.md |
| 19_ccle_cgga | RUNNING | 2026-10-02 10:50:31 |  |
| 19_ccle_cgga | DONE(78s) | 2026-10-02 10:51:49 |  | CCLE: 24 lineages, 1719 lines; DUX4~LTR-lnc rel-activation rho=-0.26 (p=5.6e-29, hi-lo=-0.05, MWU p=1.1e-23, n=1719); skip: CGGA matrices absent; summary -> /media/share/workdir/1610305236/lncProt/results/ccle_cgga/step_summaries/19_ccle_cgga.md |
| 19_ccle_cgga | RUNNING | 2026-10-02 10:52:53 |  |
| 19_ccle_cgga | DONE(78s) | 2026-10-02 10:54:11 |  | CCLE: 24 lineages, 1719 lines; DUX4~LTR-lnc rel-activation rho=-0.26 (p=5.6e-29, hi-lo=-0.05, MWU p=1.1e-23, n=1719); CGGA: 1018 samples; Grade WHO II: LTR-Δ=-0.35; Grade WHO III: LTR-Δ=-0.27; Grade WHO IV: LTR-Δ=-0.25; IDH Mutant: LTR-Δ=-0.35; IDH Wildtype: LTR-Δ=-0.27; summary -> /media/share/workdir/1610305236/lncProt/results/ccle_cgga/step_summaries/19_ccle_cgga.md |
| 20_methylation | RUNNING | 2026-10-02 10:58:30 |  |
| 20_methylation | DONE(6s) | 2026-10-02 10:58:36 |  | methylation: 765,084 probes x 29 samples; paired samples: 28; global beta vs LTR-relative activation: rho=-0.34 p=8.1e-02; global beta vs LTR-lnc absolute mean: rho=-0.03 p=8.9e-01; summary -> /media/share/workdir/1610305236/lncProt/results/ccle_cgga/step_summaries/20_methylation.md |
| 23_scrna | RUNNING | 2026-10-02 11:14:25 |  |
| 23_scrna | FAILED | 2026-10-02 11:14:27 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/23_scrna.py", line 111, in main     sb = base_sym.get(r.pc_gid.split(".")[0])    File "/media/share/h |
| 23_scrna | RUNNING | 2026-10-02 11:16:51 |  |
| 23_scrna | DONE(2s) | 2026-10-02 11:16:53 |  | matrix: 36601 genes x 765 cells; cell types: {'Neoplastic': 417, 'Oligodendrocyte': 143, 'Myeloid': 101, 'Astrocyte': 66, 'LowQC': 22, 'T_cell': 15, 'B_cell': 1}; celltype rows=28; divergent pairs present in matrix: 0; summary -> /media/share/workdir/1610305236/lncProt/results/scrna/step_summaries/23_scrna.md |
| 23_scrna | RUNNING | 2026-10-02 11:18:47 |  |
| 23_scrna | DONE(2s) | 2026-10-02 11:18:49 |  | matrix: 36601 genes x 765 cells; cell types: {'Neoplastic': 417, 'Oligodendrocyte': 143, 'Myeloid': 101, 'Astrocyte': 66, 'LowQC': 22, 'T_cell': 15, 'B_cell': 1}; celltype rows=28; divergent pairs present in matrix: 0; summary -> /media/share/workdir/1610305236/lncProt/results/scrna/step_summaries/23_scrna.md |
| 23_scrna | RUNNING | 2026-10-02 11:20:47 |  |
| 23_scrna | DONE(14s) | 2026-10-02 11:21:01 |  | matrix: 36601 genes x 765 cells; cell types: {'Neoplastic': 417, 'Oligodendrocyte': 143, 'Myeloid': 101, 'Astrocyte': 66, 'LowQC': 22, 'T_cell': 15, 'B_cell': 1}; celltype rows=28; divergent pairs present in matrix: 586; pairs tested=287; median obs/exp ratio=1.31 vs random-pair control=0.00; FDR<0.01 via p<1e-2 in 46%; summary -> /media/share/workdir/1610305236/lncProt/results/scrna/step_summaries/23_scrna.md |
| 23_scrna | RUNNING | 2026-10-02 11:21:43 |  |
| 23_scrna | RUNNING | 2026-10-02 11:38:56 |  |
| 23_scrna | DONE(54s) | 2026-10-02 11:39:50 |  | matrix: 36601 genes x 765 cells; cell types: {'Neoplastic': 417, 'Oligodendrocyte': 143, 'Myeloid': 101, 'Astrocyte': 66, 'LowQC': 22, 'T_cell': 15, 'B_cell': 1}; celltype rows=28; divergent pairs present in matrix: 586; pairs tested=287; median obs/exp ratio=1.31 vs random-pair control=1.24; FDR<0.01 via p<1e-2 in 46%; summary -> /media/share/workdir/1610305236/lncProt/results/scrna/step_summaries/23_scrna.md |
| 22_gbm_own | RUNNING | 2026-10-02 13:40:53 |  |
| 22_gbm_own | DONE(46s) | 2026-10-02 13:41:39 |  | salmon quants: 127; gene x sample matrix: (0, 127); per-sample table written; summary -> /media/share/workdir/1610305236/lncProt/results/gbm_own/step_summaries/22_gbm_own.md |
| 22_gbm_own | RUNNING | 2026-10-02 14:15:00 |  |
| 22_gbm_own | DONE(47s) | 2026-10-02 14:15:47 |  | salmon quants: 135; gene x sample matrix: (0, 135); per-sample table written; summary -> /media/share/workdir/1610305236/lncProt/results/gbm_own/step_summaries/22_gbm_own.md |
| 22_gbm_own | RUNNING | 2026-10-02 14:16:21 |  |
| 22_gbm_own | DONE(86s) | 2026-10-02 14:17:46 |  | salmon quants: 135; gene x sample matrix: (62648, 135); per-sample table written; summary -> /media/share/workdir/1610305236/lncProt/results/gbm_own/step_summaries/22_gbm_own.md |
| 22_gbm_own | RUNNING | 2026-10-02 14:18:23 |  |
| 22_gbm_own | DONE(79s) | 2026-10-02 14:19:42 |  | salmon quants: 135; gene x sample matrix: (62648, 135); per-sample table written; summary -> /media/share/workdir/1610305236/lncProt/results/gbm_own/step_summaries/22_gbm_own.md |
| 11_report_figures | RUNNING | 2026-10-02 15:30:40 |  |
| 11_report_figures | DONE(7s) | 2026-10-02 15:30:47 |  | figures written; REPORT.md written |
| 24_ml_extensions | RUNNING | 2026-10-02 15:49:11 |  |
| 24_ml_extensions | FAILED | 2026-10-02 15:50:17 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/24_ml_extensions.py", line 267, in main     C = task_c(st)    File "/media/share/workdir/1610305236/l |
| 24_ml_extensions | RUNNING | 2026-10-02 15:51:15 |  |
| 24_ml_extensions | DONE(62s) | 2026-10-02 15:52:17 |  | A: macro-F1=0.336 acc=0.610 (5-fold CV, n=19258); top features: prom_like(+0.060), phast_mean(+0.026), gc(+0.006), te_frac_LTR(+0.005), te_frac_LINE(+0.005); B: sequence-only 6-mer AUROC=0.701 (n_pos=1595, n=7975); top words: CAGTCG, CGTCAT, GCCCGA, CTTGAT, GACATC, ATTTGG; C: TCGA-LGG tumor: -0.071 (n=509); TCGA-GBM tumor: -0.067 (n=153); summary -> /media/share/workdir/1610305236/lncProt/results/ml/step_summaries/24_ml.md |
| 22_gbm_own | RUNNING | 2026-10-02 18:42:59 |  |
| 22_gbm_own | DONE(87s) | 2026-10-02 18:44:26 |  | salmon quants: 135; gene x sample matrix: (62648, 135); clinical crosswalk: 131/135 matched (sample_Name key 116, Sample key 15); abandoned: ['526393', '7801', '7802', '7812']; Grade 2: LTR-delta=-0.148 (n=19); Grade 3: LTR-delta=-0.143 (n=26); Grade 4: LTR-delta=-0.135 (n=81); grade IV vs II MWU p=3.28e-02; trend Spearman rho=0.22 p=1.1e-02 (n=126); GBM vs non-tumor: -0.139 (n=123) vs -0.123 (n=4); MWU p=3.88e-01; matched-tumor class hierarchy: PCG 69%, divergent 29%, antisense 14%, intronic 18%, sense_overlapping 29%, intergenic 13%; per-sample table written; summary -> /media/share/workdir/1610305236/lncProt/results/gbm_own/step_summaries/22_gbm_own.md |
| 16_submission_figures | RUNNING | 2026-10-02 18:53:12 |  |
| 16_submission_figures | DONE(11s) | 2026-10-02 18:53:23 |  | submission figures -> /media/share/workdir/1610305236/lncProt/figures/sub |
| 25_sample_accounting | RUNNING | 2026-10-02 19:08:42 |  |
| 25_sample_accounting | DONE(1s) | 2026-10-02 19:08:43 |  | summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/25_sample_accounting.md; rows=24 -> results/sample_accounting.tsv |
| 25_sample_accounting | RUNNING | 2026-10-02 19:14:27 |  |
| 25_sample_accounting | DONE(0s) | 2026-10-02 19:14:28 |  | summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/25_sample_accounting.md; rows=24 -> results/sample_accounting.tsv |
| 26_cns_cohorts | RUNNING | 2026-10-02 23:21:09 |  |
| 26_cns_cohorts | DONE(150s) | 2026-10-02 23:23:39 |  | MB quants: 206; MB samples=175 hierarchy: PCG 68%, divergent 30%, antisense 15%, intronic 19%, sense_overlapping 29%, intergenic 13%; MB LTR-delta median -0.144 (IQR -0.164..-0.125); SFT quants: 107; SFT samples=107 hierarchy: PCG 68%, divergent 39%, antisense 26%, intronic 28%, sense_overlapping 34%, intergenic 23%; SFT LTR-delta median -0.132; SFT grade 1: LTR-delta -0.139 (n=6); SFT grade 2: LTR-delta -0.138 (n=21); SFT grade 3: LTR-delta -0.120 (n=10); SFT III vs I MWU p=5.62e-01; trend rho=0.09 p=5.8e-01 (n=37); summary -> /media/share/workdir/1610305236/lncProt/results/cns_cohorts/step_summaries/26_cns_cohorts.md |
| 26_cns_cohorts | RUNNING | 2026-10-02 23:24:50 |  |
| 26_cns_cohorts | DONE(174s) | 2026-10-02 23:27:44 |  | MB quants: 206; MB samples=175 hierarchy: PCG 68%, divergent 30%, antisense 15%, intronic 19%, sense_overlapping 29%, intergenic 13%; MB LTR-delta median -0.144 (IQR -0.164..-0.125); SFT quants: 107; SFT samples=107 hierarchy: PCG 68%, divergent 39%, antisense 26%, intronic 28%, sense_overlapping 34%, intergenic 23%; SFT LTR-delta median -0.132; SFT grade 1: LTR-delta -0.139 (n=6); SFT grade 2: LTR-delta -0.138 (n=21); SFT grade 3: LTR-delta -0.120 (n=10); SFT III vs I MWU p=5.62e-01; trend rho=0.09 p=5.8e-01 (n=37); summary -> /media/share/workdir/1610305236/lncProt/results/cns_cohorts/step_summaries/26_cns_cohorts.md |
| 25_sample_accounting | RUNNING | 2026-10-03 17:34:12 |  |
| 25_sample_accounting | DONE(1s) | 2026-10-03 17:34:13 |  | summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/25_sample_accounting.md; rows=26 -> results/sample_accounting.tsv |
| 27_phylop | RUNNING | 2026-10-08 15:15:04 |  |
| 27_phylop | DONE(84s) | 2026-10-08 15:16:29 |  | phylop_raw rows=54533; nan mean=175; A: rows=21; significant(FDR<0.05)=21; B phyloP GC-controlled: divergent: 0.032->-0.001; antisense: -0.117->-0.027; intronic: -0.258->-0.121; sense_overlapping: -0.048->0.006; intergenic: -0.394->-0.278; C LTR contrasts: lnc/frac_accel_2=0.0317 v 0.0485 (p=1.48e-100); intergenic/frac_accel_2=0.0311 v 0.0459 (p=2.41e-59); D cross-track spearman(all)=0.806; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/27_phylop.md |
| 27_phylop | RUNNING | 2026-10-08 15:20:07 |  |
| 27_phylop | DONE(21s) | 2026-10-08 15:20:28 |  | phylop_raw rows=54533; nan mean=106; A: rows=21; significant(FDR<0.05)=21; B phyloP GC-controlled: divergent: 0.032->-0.001; antisense: -0.117->-0.027; intronic: -0.258->-0.121; sense_overlapping: -0.048->0.006; intergenic: -0.394->-0.278; C LTR contrasts: lnc/frac_accel_2=0.0319 v 0.0485 (p=3.64e-98); intergenic/frac_accel_2=0.0314 v 0.0461 (p=5.19e-58); D cross-track spearman(all)=0.805; summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/27_phylop.md |
| 16_submission_figures | RUNNING | 2026-10-08 15:28:06 |  |
| 16_submission_figures | DONE(15s) | 2026-10-08 15:28:21 |  | submission figures -> /media/share/workdir/1610305236/lncProt/figures/sub |
| 11_report_figures | RUNNING | 2026-10-08 15:28:25 |  |
| 11_report_figures | DONE(9s) | 2026-10-08 15:28:34 |  | figures written; REPORT.md written |
| 28_evo_A | RUNNING | 2026-10-08 16:02:37 |  |
| 28_evo_A | DONE(43s) | 2026-10-08 16:03:20 |  | A ablation: F0_base19: F1=0.3357 (d=+0.0000); F1_minus_phast: F1=0.3139 (d=-0.0217); F2_plus_phylop: F1=0.3669 (d=+0.0312); F3_phylop_replaces_phast: F1=0.3642 (d=+0.0285); F0 gate passed (0.33567447793175975); F2 top features: prom_like(+0.059), phast_mean(+0.027), phylop_max(+0.023), phylop_mean(+0.019), gc(+0.017), cpg_oe(+0.011); summary -> /media/share/workdir/1610305236/lncProt/results/evo_ext/step_summaries/28_a.md |
| 28_evo_B | RUNNING | 2026-10-08 16:03:25 |  |
| 28_evo_B | DONE(3s) | 2026-10-08 16:03:27 |  | B gtex lnc n=4705 q25=0.834 q75=1.043; B cage lnc n=6630 q25=19 q75=360; B key tests: gtex_lnc_overall phylop_mean: p=9.04e-05; gtex_lnc_overall tss_in_ltr: p=2.59e-03; cage_lnc_overall phylop_mean: p=1.15e-44; cage_lnc_overall tss_in_ltr: p=9.48e-26; summary -> /media/share/workdir/1610305236/lncProt/results/evo_ext/step_summaries/28_b.md |
| 28_evo_C | RUNNING | 2026-10-08 16:03:31 |  |
| 28_evo_C | DONE(773s) | 2026-10-08 16:16:25 |  | C cancers with >=20 tumor & >=5 normal: 18 (BLCA,BRCA,CHOL,COAD,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LIHC,LUAD,LUSC,PRAD,READ,STAD,THCA,UCEC); TCGA-BLCA: tumor=405 normal=19 up(lnc)=265; TCGA-BRCA: tumor=1099 normal=113 up(lnc)=174; TCGA-CHOL: tumor=35 normal=9 up(lnc)=730; TCGA-COAD: tumor=289 normal=41 up(lnc)=320; TCGA-ESCA: tumor=182 normal=13 up(lnc)=347; TCGA-GBM: tumor=153 normal=5 up(lnc)=238; TCGA-HNSC: tumor=520 normal=44 up(lnc)=297; TCGA-KICH: tumor=66 normal=25 up(lnc)=184; TCGA-KIRC: tumor=530 normal=72 up(lnc)=365; TCGA-KIRP: tumor=288 normal=32 up(lnc)=390; TCGA-LIHC: tumor=369 normal=50 up(lnc)=390; TCGA-LUAD: tumor=513 normal=59 up(lnc)=251; TCGA-LUSC: tumor=498 normal=50 up(lnc)=360; TCGA-PRAD: tumor=496 normal=52 up(lnc)=265; TCGA-READ: tumor=92 normal=10 up(lnc)=232; TCGA-STAD: tumor=410 normal=35 up(lnc)=433; TCGA-THCA: tumor=512 normal=59 up(lnc)=122; TCGA-UCEC: tumor=180 normal=23 up(lnc)=223; C tiers (lnc n=13909): tier0=11864; tier1=885; tier2=804; tier3=356; spot-check 18x6 vs 13d shards PASSED; summary -> /media/share/workdir/1610305236/lncProt/results/evo_ext/step_summaries/28_c.md |
| 28_evo_R | RUNNING | 2026-10-08 16:16:30 |  |
| 28_evo_R | DONE(0s) | 2026-10-08 16:16:30 |  | REPORT.md: appended 3 evolution-extension rows |
| 25_sample_accounting | RUNNING | 2026-10-08 16:40:02 |  |
| 25_sample_accounting | DONE(1s) | 2026-10-08 16:40:02 |  | summary -> /media/share/workdir/1610305236/lncProt/results/step_summaries/25_sample_accounting.md; rows=22 -> results/sample_accounting.tsv |
| 16_submission_figures | RUNNING | 2026-10-08 16:40:05 |  |
| 16_submission_figures | DONE(14s) | 2026-10-08 16:40:18 |  | submission figures -> /media/share/workdir/1610305236/lncProt/figures/sub |
| 16_submission_figures | RUNNING | 2026-10-08 16:44:23 |  |
| 16_submission_figures | DONE(12s) | 2026-10-08 16:44:35 |  | submission figures -> /media/share/workdir/1610305236/lncProt/figures/sub |
| 29_ext_B | RUNNING | 2026-10-08 17:09:26 |  |
| 29_ext_B | DONE(0s) | 2026-10-08 17:09:26 |  | B: age medians phylop ancient=0.028; primate_old=-0.006; primate_young=-0.036; LTR subset rho(age, phylop)=-0.173; summary -> /media/share/workdir/1610305236/lncProt/results/promoter_ext/step_summaries/29_b.md |
| 29_ext_C | RUNNING | 2026-10-08 17:09:28 |  |
| 29_ext_C | DONE(1s) | 2026-10-08 17:09:29 |  | C eta^2 gain from divergent split: phylop_mean: +0.0821 (+677%); phast_mean: +0.1099 (+1021%); tau: +0.1823 (+557%); cpg_oe: +0.2278 (+823%); te_cov: +0.1252 (+1071%); summary -> /media/share/workdir/1610305236/lncProt/results/promoter_ext/step_summaries/29_c.md |
| 29_ext_A | RUNNING | 2026-10-08 17:09:48 |  |
| 29_ext_A | DONE(5s) | 2026-10-08 17:09:53 |  | A: top_width x cpg_oe rho=0.636 p=0.0e+00; PCG median width=45.0 n_peaks=3.0; divergent 49.0 / intergenic 0.0; summary -> /media/share/workdir/1610305236/lncProt/results/promoter_ext/step_summaries/29_a.md |
| 29_ext_A | RUNNING | 2026-10-08 17:10:52 |  |
| 29_ext_A | DONE(5s) | 2026-10-08 17:10:57 |  | A: top_width x cpg_oe rho=0.636 p=0.0e+00; PCG median width=48.0 n_peaks=3.0; divergent 49.0 / intergenic 0.0; summary -> /media/share/workdir/1610305236/lncProt/results/promoter_ext/step_summaries/29_a.md |
| 30_perturb | RUNNING | 2026-10-08 17:19:28 |  |
| 30_perturb | FAILED | 2026-10-08 17:19:28 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/30_perturb.py", line 135, in main     r = prep_response(c.reset_index(), t, k, "x")    File "/media/s |
| 30_perturb | RUNNING | 2026-10-08 17:19:55 |  |
| 30_perturb | DONE(3s) | 2026-10-08 17:19:58 |  | perturb LTR-vs-other: GSE56787_FSHD_vs_ctrl: med 0.003339438487616762 v 0.0032853939904249757, p=5.82e-01; GSE337536_decitabine: med 0.07313559813971784 v 0.0673924779010853, p=1.88e-03; GSE305426_DNMTi_vs_DMSO: med 0.00566456354210873 v 0.0018769795231018166, p=1.07e-25; summary -> /media/share/workdir/1610305236/lncProt/results/perturb/step_summaries/30_perturb.md |
| 30_perturb | RUNNING | 2026-10-08 17:23:52 |  |
| 30_perturb | DONE(2s) | 2026-10-08 17:23:54 |  | perturb LTR-vs-other: GSE56787_FSHD_vs_ctrl: med 0.003339438487616762 v 0.0032853939904249757, p=5.82e-01; GSE337536_decitabine: med 0.07313559813971784 v 0.0673924779010853, p=1.88e-03; GSE305426_DNMTi_vs_DMSO: med 0.00566456354210873 v 0.0018769795231018166, p=1.07e-25; summary -> /media/share/workdir/1610305236/lncProt/results/perturb/step_summaries/30_perturb.md |
| 31_tf_A | RUNNING | 2026-10-09 11:20:38 |  |
| 31_tf_A | DONE(50s) | 2026-10-09 11:21:28 |  | A: pairs=950 (top: antisense:OLIG1+BHLHE23; antisense:OLIG2+BHLHE22; antisense:BHLHA15+BHLHE22); cofactor divergent EP300=0.92 intergenic EP300=0.32; summary -> /media/share/workdir/1610305236/lncProt/results/tf_contexts/step_summaries/31_a.md |
| 31_tf_T | RUNNING | 2026-10-09 11:21:43 |  |
| 31_tf_T | FAILED | 2026-10-09 11:22:21 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/31_tf_contexts.py", line 138, in mode_t     tf_expr = pd.DataFrame({t: d[i, gtex.index].mean(axis=0)  |
| 31_tf_T | RUNNING | 2026-10-09 11:40:11 |  |
| 31_tf_T | FAILED | 2026-10-09 11:40:38 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/31_tf_contexts.py", line 144, in mode_t     need_rows = np.array(sorted(set(tf_idx.values())  TypeErr |
| 31_tf_T | RUNNING | 2026-10-09 11:55:41 |  |
| 31_tf_T | FAILED | 2026-10-09 11:56:08 | Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/31_tf_contexts.py", line 144, in mode_t     need_rows = np.array(sorted(set(tf_idx.values())  TypeErr |
| 31_tf_T | RUNNING | 2026-10-09 13:32:04 |  |
| 31_tf_T | DONE(35s) | 2026-10-09 13:32:39 |  | TF panel mapped: 441/884; T: gtex couplings sig=990; top: HINFP antisense rho=0.93; ZNF211 antisense rho=0.91; ZNF692 antisense rho=0.89; HINFP intronic rho=0.89; POU2F1 sense_overlapping rho=0.88; summary -> /media/share/workdir/1610305236/lncProt/results/tf_contexts/step_summaries/31_t.md |
| 30_perturb | RUNNING | 2026-10-09 13:38:24 |  |
| 30_perturb | FAILED | 2026-10-09 13:38:32 |   File "/media/share/workdir/1610305236/lncProt/scripts/30_perturb.py", line 131, in gse85461_fastq     f = d.merge(lab, on="gene_id")    File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packages/pandas/core/frame.py", line 10859, in merge     return merge(    F |
| 30_perturb | RUNNING | 2026-10-09 13:39:05 |  |
| 30_perturb | FAILED | 2026-10-09 13:39:12 |   File "/media/share/workdir/1610305236/lncProt/scripts/30_perturb.py", line 131, in gse85461_fastq     f = d.merge(lab, on="gene_id")    File "/media/share/home/1610305236/.local/share/mamba/envs/Jupyter/lib/python3.10/site-packages/pandas/core/frame.py", line 10859, in merge     return merge(    F |
| 30_perturb | RUNNING | 2026-10-09 13:39:59 |  |
| 30_perturb | FAILED | 2026-10-09 13:40:07 |  The above exception was the direct cause of the following exception:   Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/30_perturb.py", line 223, in  |
| 30_perturb | RUNNING | 2026-10-09 13:40:33 |  |
| 30_perturb | FAILED | 2026-10-09 13:40:41 |  The above exception was the direct cause of the following exception:   Traceback (most recent call last):    File "/media/share/workdir/1610305236/lncProt/scripts/pylib/progress.py", line 93, in step     yield st    File "/media/share/workdir/1610305236/lncProt/scripts/30_perturb.py", line 223, in  |
| 30_perturb | RUNNING | 2026-10-09 13:41:07 |  |
| 30_perturb | DONE(8s) | 2026-10-09 13:41:15 |  | perturb LTR-vs-other: GSE56787_FSHD_vs_ctrl: med 0.003339438487616762 v 0.0032853939904249757, p=5.82e-01; GSE337536_decitabine: med 0.07313559813971784 v 0.0673924779010853, p=1.88e-03; GSE305426_DNMTi_vs_DMSO: med 0.00566456354210873 v 0.0018769795231018166, p=1.07e-25; GSE85461_iDUX4_fastq: med 0.0 v 0.0, p=4.26e-10; GSE85461_lentiDUX4_14h: med 0.0 v 0.0, p=3.90e-08; summary -> /media/share/workdir/1610305236/lncProt/results/perturb/step_summaries/30_perturb.md |
| 30_perturb | RUNNING | 2026-10-09 13:41:37 |  |
| 30_perturb | DONE(7s) | 2026-10-09 13:41:44 |  | perturb LTR-vs-other: GSE56787_FSHD_vs_ctrl: med 0.003339438487616762 v 0.0032853939904249757, p=5.82e-01; GSE337536_decitabine: med 0.07313559813971784 v 0.0673924779010853, p=1.88e-03; GSE305426_DNMTi_vs_DMSO: med 0.00566456354210873 v 0.0018769795231018166, p=1.07e-25; GSE85461_iDUX4_fastq: med 0.0 v 0.0, p=4.26e-10; GSE85461_lentiDUX4_14h: med 0.0 v 0.0, p=3.90e-08; summary -> /media/share/workdir/1610305236/lncProt/results/perturb/step_summaries/30_perturb.md |
| 30_perturb | RUNNING | 2026-10-09 13:42:10 |  |
| 30_perturb | DONE(7s) | 2026-10-09 13:42:17 |  | perturb LTR-vs-other: GSE56787_FSHD_vs_ctrl: med 0.003339438487616762 v 0.0032853939904249757, p=5.82e-01; GSE337536_decitabine: med 0.07313559813971784 v 0.0673924779010853, p=1.88e-03; GSE305426_DNMTi_vs_DMSO: med 0.00566456354210873 v 0.0018769795231018166, p=1.07e-25; GSE85461_iDUX4_fastq: med 0.0 v 0.0, p=4.26e-10; GSE85461_lentiDUX4_14h: med 0.0 v 0.0, p=3.90e-08; summary -> /media/share/workdir/1610305236/lncProt/results/perturb/step_summaries/30_perturb.md |

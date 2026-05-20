Plan: Germline-Informed Reverse Translation Package

 Context

 abstar requires nucleotide input, but researchers often have only amino acid antibody sequences (e.g., from mass spec or translated proteomics). The
 goal is a companion package that reverse-translates AA sequences to NT sequences in the most germline-faithful way possible — so that downstream abstar
 output (V/J assignment, mutation counts, CDR/FWR regions) is biologically meaningful, not just an artifact of arbitrary codon choice.

 The user described one valid pipeline: naive reverse-translate → run abstar → align NT to germline → revert all silent mutations. This plan proposes a
 single-pass direct algorithm that achieves the same provably optimal result without requiring an abstar call in the loop — critical for scaling to
 millions of sequences.

 ---
 Core Mathematical Insight

 The optimization problem ("find the minimum-Hamming-distance NT sequence consistent with this AA sequence") decomposes into independent per-position
 sub-problems because:
 - Codons don't overlap
 - Hamming distance is additive across positions
 - Global minimum = sum of per-position minima

 For each codon position i aligned to a germline:
 optimal_codon = argmin_{c : translate(c) == target_aa[i]}  Hamming(c, germline_codon[i])
 Ties broken by human codon frequency. This is pre-computable as a 1280-entry lookup table (20 AAs × 64 possible germline codons). O(1) per position at
 runtime.

 ---
 Algorithm

 1. Batch MMseqs2 protein-protein search of all input AA sequences against
    V germline AA database  →  V assignment per sequence
 2. Build post-V AA FASTA (aa_seq[v_qend:] per sequence)
 3. Batch MMseqs2 protein-protein search against J germline AA database
    →  J assignment per sequence
 4. Parallel per-sequence reconstruction (ProcessPoolExecutor):
    a. 5' overhang (before V alignment): human codon frequency
    b. V region:  optimal_codon(target_aa, germline_codon) per position
    c. CDR3 (between V end and J start): human codon frequency
    d. J region:  optimal_codon(target_aa, germline_codon) per position
    e. 3' overhang (after J alignment): human codon frequency
 5. Return list[abutils.Sequence] (NT) → feed directly into abstar.run()

 ---
 Package Structure (standalone, pip install ab_revtrans)

 abverse/
 ├── pyproject.toml
 ├── abverse/
 │   ├── __init__.py          # exports: reverse_translate, build_germline_aa_db
 │   ├── _codons.py           # lookup table, hamming helper, fallback codon
 │   ├── _germline_db.py      # build + cache AA germline databases
 │   ├── _search.py           # MMseqs2 protein-protein search wrappers + result parsing
 │   ├── _reconstruct.py      # per-sequence NT reconstruction (pure, picklable)
 │   ├── _pipeline.py         # orchestration: input normalization → output
 │   └── tests/
 │       ├── test_codons.py
 │       ├── test_germline_db.py
 │       ├── test_reconstruct.py
 │       └── test_pipeline.py

 Companion abstar/pp.py addition (graceful import):
 try:
     from abverse import reverse_translate
 except ImportError:
     pass

 ---
 Key Modules

 _codons.py

 - CODON_TABLE: 64→20+stops (from abutils/utils/codons.py)
 - SYNONYMOUS_CODONS: inverted, built at import
 - HUMAN_PREFERRED_CODON: dict[aa → most_frequent_human_codon] — hardcoded from DNAchisel's h_sapiens table (so DNAchisel is build-time only, not
 runtime)
 - build_optimal_codon_lookup() → dict[(target_aa, germline_codon), optimal_codon] — built once at import, 1280 entries
 - optimal_codon(target_aa, germline_codon) → str — O(1) lookup
 - fallback_codon(target_aa) → str — most frequent human codon (used for CDR3/overhangs)

 _germline_db.py

 Translates abstar's ungapped NT germlines → AA FASTA; builds MMseqs2 protein databases; caches under ~/.abverse/germline_dbs/.

 Cache layout:
 ~/.abverse/germline_dbs/bcr/human/
   v_aa/  sequences.fasta | mmseqs_db | frame_map.json | checksum.sha256
   j_aa/  sequences.fasta | mmseqs_db | frame_map.json | checksum.sha256

 Cache invalidation: SHA-256 of source NT FASTA. Rebuilt if abstar updates germlines.

 frame_map.json: V genes always frame 1 (all 342 human V genes verified). J genes: detect frame by finding conserved WG.G (IGH) or FG.G (IGK/IGL) motif;
 fallback = longest stop-codon-free translation.

 Source NT files:
 - abstar/germline_dbs/bcr/human/ungapped/v.fasta (342 genes)
 - abstar/germline_dbs/bcr/human/ungapped/j.fasta (23 genes)
 - Same for IGK, IGL chains

 Key functions:
 - build_germline_aa_db(species, receptor, force_rebuild) → dict — returns v/j DB paths + NT seq dicts + frame maps
 - get_germline_aa_db_path(species, receptor, segment) → str
 - _build_aa_fasta(nt_fasta, output_fasta, frame_map) — translates + detects frames
 - _make_mmseqs_protein_db(fasta, db_prefix) — runs mmseqs createdb

 _search.py

 MMseqs2 protein-protein search wrappers (mirrors abstar/assigners/mmseqs.py structure).

 V search parameters:
 abutils.tl.mmseqs_search(search_type=1, max_seqs=25, sensitivity=7.5,
     format_output="query,target,nident,qstart,qend,tstart,tend,qseq,tseq",
     additional_cli_args="--min-aln-len 10 --alignment-mode 3")

 J search parameters (looser — J genes are short):
 abutils.tl.mmseqs_search(search_type=1, max_evalue=1000.0,
     additional_cli_args="--min-aln-len 5 --alignment-mode 3")

 Result parsing via Polars: sort by nident desc, keep best hit per query (unique(subset=["query"], keep="first")).

 Key functions:
 - search_v_germline(query_fasta, v_db_path, output_path, threads) → pl.DataFrame
 - build_j_query_fasta(v_results, input_aa_seqs, output_path) — writes aa_seq[v_qend:] per sequence
 - search_j_germline(j_query_fasta, j_db_path, output_path, threads) → pl.DataFrame
 - merge_vj_results(v_df, j_df) → pl.DataFrame

 Coordinate note: MMseqs2 easy-search reports qstart/qend as 0-indexed inclusive. v_qend must be treated as inclusive in parsing; j_qstart is relative to
  the post-V query subsequence and must be offset by v_qend + 1 to get absolute AA position.

 _reconstruct.py

 Pure functions (no I/O), picklable for ProcessPoolExecutor.

 Main function:
 def reconstruct_sequence(
     seq_id, aa_seq,
     v_call, v_qstart, v_qend, v_tstart,
     j_call, j_qstart, j_qend, j_tstart,    # j coords relative to post-V region
     v_nt_seqs, j_nt_seqs, v_frame_map, j_frame_map,
 ) -> abutils.Sequence

 Internal coordinate math:
 # V region
 v_germ_offset = (v_frame - 1) + (v_tstart + i) * 3
 germ_codon = v_nt[v_germ_offset : v_germ_offset + 3]

 # J region (j_qstart is relative to post-V query)
 j_abs_start = (v_qend + 1) + j_qstart
 j_germ_offset = (j_frame - 1) + (j_tstart + i) * 3
 germ_codon = j_nt[j_germ_offset : j_germ_offset + 3]

 Validation: assert abutils.tl.translate(nt_out) == aa_in before returning.

 Returns abutils.Sequence(nt_out, id=seq_id) with annotations: {v_call, j_call, reconstruction_method}.

 Edge cases:

 ┌──────────────────────────────────┬────────────────────────────────────────────────────────────┐
 │            Situation             │                          Handling                          │
 ├──────────────────────────────────┼────────────────────────────────────────────────────────────┤
 │ No V assignment                  │ fallback_codon for all positions; method='codon_frequency' │
 ├──────────────────────────────────┼────────────────────────────────────────────────────────────┤
 │ No J assignment                  │ V lookup for V region; fallback for rest                   │
 ├──────────────────────────────────┼────────────────────────────────────────────────────────────┤
 │ 5'/3' overhangs                  │ fallback_codon                                             │
 ├──────────────────────────────────┼────────────────────────────────────────────────────────────┤
 │ Germline codon truncated at edge │ fallback_codon                                             │
 ├──────────────────────────────────┼────────────────────────────────────────────────────────────┤
 │ Non-standard AA (X, B, Z)        │ 'NNN' (don't crash)                                        │
 ├──────────────────────────────────┼────────────────────────────────────────────────────────────┤
 │ Stop codon in input AA           │ ValueError with position and seq_id                        │
 ├──────────────────────────────────┼────────────────────────────────────────────────────────────┤
 │ V/J alignment overlap            │ V takes priority; J starts after V end                     │
 └──────────────────────────────────┴────────────────────────────────────────────────────────────┘

 _pipeline.py

 def reverse_translate(
     sequences,              # FASTA path | list[str] | list[abutils.Sequence]
     species="human",
     receptor="bcr",
     n_processes=None,       # defaults to cpu_count()
     threads=None,           # MMseqs2 threads
     chunksize=500,
     force_rebuild_db=False,
     output_fasta=None,
     verbose=False,
 ) -> list[abutils.Sequence]

 Steps: normalize input → build/load germline AA DBs → write AA query FASTA → MMseqs2 V search → build J query FASTA → MMseqs2 J search → merge V/J →
 dispatch _reconstruct_batch via ProcessPoolExecutor(mp_context=spawn) → collect + sort by input order → optional FASTA write → return.

 All temp files cleaned up in try/finally.

 ---
 Dependencies

 abutils >= 0.5.1   # Sequence, mmseqs_search, translate, alignment, io
 polars >= 0.20     # batch DataFrame parsing of MMseqs output
 DNAchisel is build-time only (for generating the embedded HUMAN_PREFERRED_CODON table). Not a runtime dependency.

 ---
 Integration

 import abverse
 import abstar

 aa_seqs = abutils.io.read_fasta("antibodies_aa.fasta")
 nt_seqs = abverse.reverse_translate(aa_seqs)   # list[Sequence] with NT
 results  = abstar.run(nt_seqs)                      # normal abstar annotation

 ---
 Implementation Sequencing

 1. _codons.py — implement + test all 1280 lookup entries and known examples (no external deps)
 2. _germline_db.py — translate + frame-detect, verify all V/J genes produce clean AA, build MMseqs DB, test cache invalidation
 3. _search.py — V + J protein search, result parsing, coordinate validation on known antibodies
 4. _reconstruct.py — per-sequence reconstruction, all edge cases, validate translate(out) == in for all tests
 5. _pipeline.py — orchestration, input normalization, parallel dispatch
 6. Integration test — 100 known human BCR AA sequences: reverse_translate → abstar.run() → verify V/J calls match expected germline assignments
 7. Scaling benchmark — 10k sequences, measure throughput, ensure no iterative abstar calls in critical path

 ---
 Anticipated Challenges

 - J-gene frame detection: The motif heuristic (WG.G / FG.G) covers canonical alleles; unusual alleles may need the longest-stop-free fallback. Validate
 against all 23 human J genes.
 - MMseqs coordinate off-by-one: Confirm qend inclusive/exclusive by printing aa_seq[v_qend] for a known V gene and verifying it matches the expected
 last V residue. Adjust +1 in reconstruct_sequence() accordingly.
 - Short CDR3 / overlap: When CDR3 is 0–1 AA, j_abs_start ≤ v_abs_end. Guard required; V coordinates are authoritative.

 ---
 Critical Source Files

 - abstar/assigners/mmseqs.py — mirror the search/parse pattern
 - abstar/annotation/germline.py — get_germline_database_path() to locate source NT FASTs
 - abstar/germline_dbs/bcr/human/ungapped/{v,j}.fasta — source NT germlines to translate
 - /opt/conda/lib/python3.12/site-packages/abutils/utils/codons.py — codon table source
 - /opt/conda/lib/python3.12/site-packages/abutils/tools/search.py — mmseqs_search signature
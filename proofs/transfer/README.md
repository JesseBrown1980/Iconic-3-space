# Transfer 1: finite subgroups of SU(2) reproduce continuous-SU(2) lattice gauge theory through their design degree

**Status: proven, with corrections.** An independent second reader checked the proof line by line ([`second-seat/CHECK.md`](second-seat/CHECK.md)) and gave the verdict PROVEN_WITH_CORRECTIONS. The corrections sit beside the theorem in [`CORRECTIONS.hbp`](CORRECTIONS.hbp). The full statement and proof are in [`TRANSFER-THEOREM.md`](TRANSFER-THEOREM.md).

## The theorem

**Setting.**
- A finite lattice whose plaquettes are closed 4-step words.
- Link variables in K, which is either SU(2) with Haar measure or a finite subgroup G with the uniform measure.
- Wilson weight exp(β Σ_p tr(U_p)/2), and Z_K(β) as a power series in β.

**Hypotheses.**
- **H1:** every plaquette word runs through four distinct links.
- **H2′:** inside the star of each link l, every plaquette through l has a private link, that is, a link l′ ≠ l that occurs exactly once in the plaquette and lies in no other plaquette through l.

Distinct plaquettes sharing at most one link implies H2′, given H1. Open hypercubic boxes satisfy H1 and H2′, and so do periodic ones with every side ≥ 2.

**Claims.** If G is a t-design on SU(2), that is, its moments M_n equal the Haar moments for n ≤ t:
- **(A)** Under H1: log Z_G and log Z_SU(2) agree through βᵗ.
- **(B)** Under H1 and H2′: at β^(t+1) they differ by exactly N_p · ι_G / (2^(t+1) (t+1)!). Here ι_G = M^G_(t+1) − M^Haar_(t+1) is a positive integer: the number of independent G-invariants of degree t + 1. Only single-plaquette powers contribute.
- **(C)** An observable of link-degree d agrees through β^(t−d).
- **(C′)** Under H1 and H2′: each plaquette expectation first differs at βᵗ, by ι_G / (2^(t+1) t!).
- **(D)** The same holds between two finite subgroups, with t the smaller of their two design degrees.

**On the ladder.** Q8 and Dic₃ reach t = 3, 2T t = 5, 2O t = 7 and 2I t = 11. Dic₃ is the SU(2) rung over the finite repository's S₃.

## What was checked

- **Exact enumeration.** The prover's Rust 1.81 enumerator, gated with `clippy -D warnings`, ran with exact rational bookkeeping on:
  - the cube surface;
  - the 2×2 and 3×3 tori;
  - the 2×1 strip and the 1×1 torus;
  - a two-plaquette "pillow";
  - an open 1×1×2 stack, with up to 5,159,780,352 configurations.

  The Haar side comes from the character formula, validated against enumeration. All 106 checks on lattices that satisfy the hypotheses hold. Three lattices that break a hypothesis show why each hypothesis is needed.
- **The cube.** On the cube, Dic₃ and 2T first differ at β⁴ by 6/384 = 1/64. GPT-6's review derived the same value analytically.
- **The second reader.** It recomputed all 45 frozen pair rows from its own moments, with 0 mismatches. It also ran its own single-plaquette enumeration.

## The frozen test still to run

[`../../predictions/transfer/`](../../predictions/transfer/) holds the prediction file (`8ada45c8…`). It was written in the same run as the enumeration, so the rows checked there were not frozen first ([`FREEZE.hbp`](../../predictions/transfer/FREEZE.hbp) records this). The rows still unmeasured become a public frozen test from this commit on:
- every row on the 4×4 torus;
- 14 of the 15 rows on the 3×3 torus;
- every 2I lattice row.

## Scope

The transfer carries exact strong-coupling Taylor coefficients at β = 0 on every fixed finite lattice satisfying the hypotheses, with a computable first deviation. Bounds at a fixed β, orders above t + 1, and the weak-coupling regime where the continuum limit lives are OPEN.

## Files

- **`check/`:** the prover's checker and enumerator, with every job's input and output.
- **`second-seat/`:** the second reader's script, report and receipt.

Both scripts use Python's exact `fractions` module, so both folders are pinned as evidence in [`../../tools/FLOAT-EVIDENCE.txt`](../../tools/FLOAT-EVIDENCE.txt). An integer-pair port is the follow-up.

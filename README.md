# QEAR: Quantum Exposure Assessment and Ranking

Analysis code for the article *Quantum Exposure Assessment and Ranking of Lightweight
Authentication Protocols* (submitted to Quantum Information Processing).

The code derives, in closed form, the quantum resources needed for a Grover key-recovery
attack on five ultralightweight RFID authentication protocols (LMAP, EMAP, M2AP, SASI,
URAP) and an Ascon-Hash256-based scheme. It produces every table and figure in the article.

## Repository

https://github.com/khayyam2302/qear-quantum-exposure

## Requirements

Python 3.9 or later.

```bash
git clone https://github.com/khayyam2302/qear-quantum-exposure.git
cd qear-quantum-exposure
pip install -r requirements.txt
```

## Reproduce the results

```bash
python3 qear_model.py   # prints all results, writes results/numbers.json
python3 make_figs.py    # writes Fig3 to Fig8 to figures/ (PDF and 600 dpi PNG), numbered as in the article
```

Figures 1 and 2 are TikZ diagrams in `figures/tikz/`. Build them with `pdflatex`.

The sensitivity analysis uses the fixed seed 20260916, so every run gives identical numbers.

## Files

| File | Contents |
|---|---|
| `qear_model.py` | Cost rules, protocol models, attack cost (Theorem 1), QES, sensitivity analysis, NIST margins, minimum word size, surface-code model |
| `make_figs.py` | Figures 3 to 8 |
| `figures/tikz/` | LaTeX sources of Figures 1 and 2 |
| `results/numbers.json` | All numbers reported in the article |

## Where each result comes from

| Article item | Source in `qear_model.py` |
|---|---|
| Table 5 (cost rules) | `op_cost`, `ascon_round` |
| Tables 12, 13, 15, 16 (resources) | `ulw_circuit`, `ascon_circuit` |
| Tables 17, 18, 21 (QES, leave-one-out, sensitivity) | `qes`, main block |
| Tables 19, 20 (attack cost, rankings) | `attack` |
| Physical resources, NIST margins, minimum word size | `surface`, `nist`, `machines`, main block |

## Changing the assumptions

* Word size: argument `b` of `ulw_circuit` (default 32).
* Ascon permutation calls: `ASCON_CALLS` (default 8).
* Rotation costing: `rot="ROT"` (with rotate-back), `"ROT_NOUNDO"`, or `"ROT_STATIC"`.

## Citation

If you use this code, please cite the article (details to be added on publication) and
this repository through its Zenodo DOI.

## License

MIT, see `LICENSE`.

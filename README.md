# genome-atlas

[![CI](https://github.com/ahmedanees-m/genome-atlas/workflows/CI/badge.svg)](https://github.com/ahmedanees-m/genome-atlas/actions)
[![codecov](https://codecov.io/gh/ahmedanees-m/genome-atlas/branch/main/graph/badge.svg)](https://codecov.io/gh/ahmedanees-m/genome-atlas)
[![PyPI](https://img.shields.io/pypi/v/genome-atlas.svg)](https://pypi.org/project/genome-atlas/)
[![Docs](https://img.shields.io/badge/docs-readthedocs-blue)](https://genome-atlas.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A unified knowledge graph and embedding space for programmable genome-writing enzymes — CRISPR, CAST, bridge recombinases, Fanzor, and more — with cross-system selection decision support.

Part of the [PEN-STACK](https://github.com/ahmedanees-m/pen-stack) infrastructure for non-destructive genome engineering.

## Install

```bash
pip install genome-atlas
```

## Quickstart

```python
from genome_atlas import Atlas

atlas = Atlas.load()
recs = atlas.select_editor(
    cell_type="HEK293T",
    edit_type="insertion",
    cargo_size_bp=1500,
    delivery="AAV",
    top_k=5,
)
for r in recs:
    print(r.system, r.pen_score, r.aav_fit)
```

## Citation

See [CITATION.cff](CITATION.cff) or cite the preprint:
> Ahmed, A. (2026). *GENOME-ATLAS: A unified knowledge graph for programmable genome-writing enzymes.* bioRxiv.

## License

MIT — see [LICENSE](LICENSE).

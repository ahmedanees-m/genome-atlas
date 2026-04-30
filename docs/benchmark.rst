Benchmark Results
=================

Primary Benchmark — Protein → Domain Link Prediction (v0.6.0)
--------------------------------------------------------------

All results use an 80/10/10 train/val/test split with type-consistent negative
sampling. Bootstrap 95% CIs: 1,000 resamples, seed 42.

.. list-table::
   :header-rows: 1
   :widths: 30 20 30 20

   * - Model
     - AUROC
     - 95% CI
     - AUPRC
   * - Node2Vec (inductive)
     - **0.9868**
     - [0.9806, 0.9921]
     - 0.9639
   * - GraphSAGE (2-layer, ESM-2 init)
     - 0.9707
     - [0.9627, 0.9780]
     - 0.9721
   * - GAT (residual connections)
     - 0.9685
     - [0.9597, 0.9770]
     - 0.9663
   * - Node2Vec (transductive †)
     - 0.9965
     - —
     - —
   * - Quantum kernel ‡
     - 0.8847
     - —
     - —
   * - Classical RBF-SVM ‡
     - 0.8761
     - —
     - —

† Transductive Node2Vec uses full-graph walks including test edges; reported as
a topology upper-bound only (supplementary material).

‡ Quantum kernel and RBF-SVM evaluated on a subgraph sample; not directly
comparable to inductive GNN results.

Selection Validation
--------------------

Top-3 accuracy on 10 published therapeutic scenarios: **70%** (7/10).
See `VALIDATION.md <https://github.com/ahmedanees-m/genome-atlas/blob/main/VALIDATION.md>`_
for the full per-scenario breakdown and miss analysis.

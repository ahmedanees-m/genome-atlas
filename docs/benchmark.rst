Benchmark Results
=================

Primary Benchmark (Protein → Domain Link Prediction)
----------------------------------------------------

.. list-table::
   :header-rows: 1

   * - Model
     - AUROC
     - AUPRC
   * - GAT (1 head, residual)
     - 0.9705 [0.9446–0.9964]
     - 0.9421
   * - GraphSAGE (2-layer)
     - 0.9664 [0.9405–0.9923]
     - 0.9184
   * - Classical RBF
     - 0.9331
     - —
   * - Quantum Kernel
     - 0.8731
     - —
   * - Node2Vec
     - 0.8202 [0.8052–0.8342]
     - 0.8905

Topology Consistency Checks
---------------------------

Structure → Protein edges yield AUROC ≈ 0.995–0.997, reflecting near-deterministic
topology (mean structure degree 0.72).

Selection Validation
--------------------

Top-3 accuracy on 10 published therapeutic scenarios: **70%** (7/10).

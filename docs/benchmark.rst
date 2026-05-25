Benchmark Results
=================

Primary Benchmark - Protein -> Domain Link Prediction (v0.7.1 / v0.7.2)
------------------------------------------------------------------------

All results use an 80/10/10 train/val/test split with type-consistent negative
sampling. Bootstrap 95% CIs: 1,000 resamples, seed 42. Test set: n~1,908
(positive + matched negatives). Graph: 28 foundational systems (ISCro4 added
in v0.7.1; canonical name adopted in v0.7.2), 13,401 nodes, 11,817 edges.
Rebuilt 2026-05-23.

.. list-table::
   :header-rows: 1
   :widths: 30 20 30 20

   * - Model
     - AUROC
     - 95% CI
     - AUPRC
   * - GraphSAGE (2-layer, ESM-2 init)
     - **0.9714**
     - [0.9625, 0.9797]
     - **0.9451**
   * - GAT (residual connections, 4-head)
     - 0.9690
     - [0.9590, 0.9778]
     - 0.9331
   * - Node2Vec (inductive †)
     - 0.9890
     - [0.9825, 0.9940]
     - 0.9675
   * - Node2Vec (transductive ‡)
     - 0.9965
     - [0.9924, 0.9997]
     - -
   * - Quantum kernel §
     - 0.8731
     - -
     - 0.8429
   * - Classical RBF-SVM §
     - 0.9331
     - -
     - 0.8984

GraphSAGE and GAT are statistically tied (AUROC delta = 0.0024; 95% CIs
substantially overlap). Both models remain within v0.7.0 CIs after adding
ISCro4 (28th foundational system; preprint label IS622 deprecated in v0.7.2).
v0.7.0 CIs: GraphSAGE [0.9629, 0.9800], GAT [0.9593, 0.9775].

† Node2Vec inductive: random walks on train-split only (test/val HAS_DOMAIN
edges withheld). Topology without ESM-2 features. Excluded from primary Table 1
because input modality (topology-only) differs from inductive GNNs (ESM-2
+ topology). Reported as a supplementary topology-only baseline.

‡ Transductive Node2Vec uses full-graph walks including test edges; reported as
a topology upper-bound only (supplementary material). Incomparable to inductive
GNN results.

§ Quantum kernel and RBF-SVM evaluated on a subgraph sample; not directly
comparable to inductive GNN results. Carried forward from v0.6.0.

Secondary Task - Structure -> Protein Link Prediction
-----------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 20 30

   * - Model
     - AUROC
     - 95% CI
   * - GraphSAGE
     - **0.9971**
     - [0.9913, 1.0000]
   * - GAT
     - **0.9971**
     - [0.9913, 1.0000]
   * - Node2Vec
     - 0.9739
     - [0.9538, 0.9891]

Selection Validation
--------------------

Top-3 accuracy on 10 published therapeutic scenarios: **70%** (7/10).
See `VALIDATION.md <https://github.com/ahmedanees-m/genome-atlas/blob/main/VALIDATION.md>`_
for the full per-scenario breakdown and miss analysis.

Authoritative CI file: ``reproduction/bootstrap_ci_v7.json`` as raw data.

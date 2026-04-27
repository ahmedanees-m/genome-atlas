GENOME-ATLAS Documentation
==========================

GENOME-ATLAS is a programmable knowledge graph for genome-writing systems,
combining structural biology, protein language models, graph neural networks,
and quantum kernel methods to support therapeutic editor selection.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api
   selection
   benchmark

Quick Start
-----------

.. code-block:: python

   from genome_atlas.api import Atlas

   atlas = Atlas(
       graph_path="atlas.gpickle",
       embeddings_path="embeddings.parquet",
       targets_path="targets.parquet",
   )

   # Query a system
   print(atlas.query_system("System_SpCas9"))

   # Select an editor for a therapeutic scenario
   recs = atlas.select_editor(
       cell_type="HEK293T",
       edit_type="deletion",
       cargo_size_bp=0,
       delivery="AAV",
       prefer_dsb_free=True,
       top_k=5,
   )
   for r in recs:
       print(r.system, r.pen_score)

CLI Usage
---------

.. code-block:: bash

   genome-atlas --help
   genome-atlas query-system System_SpCas9
   genome-atlas select --cell HEK293T --edit deletion --top-k 5

Modules
-------

* :mod:`genome_atlas.api` — Public query API
* :mod:`genome_atlas.selection` — Therapeutic editor selection engine
* :mod:`genome_atlas.models.graphsage` — Heterogeneous GNN models
* :mod:`genome_atlas.graph.build` — PyG graph builder

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

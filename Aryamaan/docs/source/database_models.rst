.. _database-schema:

Database Schema
===============

.. toctree::
   :maxdepth: 2

Tables
------

.. _model_input-table:

model_input Table
~~~~~~~~~~~~~~~~~

The `model_input` table stores information about input models.

+-----------------------+------------------+----------------------------------------------------------+
| Field                 | Type             | Description                                              |
+=======================+==================+==========================================================+
| id                    | UUID (Primary)   | Unique identifier, md5 hash of                           |
|                       |                  | "{scaffold_SMILE}_{log_d_min}_{log_d_max}_{chem_param}". |
+-----------------------+------------------+----------------------------------------------------------+
| scaffold_SMILE        | Text (Not Null)  | The scaffold SMILE (Simplified Molecular                 |
|                       |                  | Input Line Entry) of the model input.                    |
+-----------------------+------------------+----------------------------------------------------------+
| log_p_min             | Float (Not Null) | The minimum logP (partition coefficient)                 |
|                       |                  | value for the model input.                               |
+-----------------------+------------------+----------------------------------------------------------+
| log_p_max             | Float (Not Null) | The maximum logP (partition coefficient)                 |
|                       |                  | value for the model input.                               |
+-----------------------+------------------+----------------------------------------------------------+

.. note:: The `id` field is a primary key and must be unique. The `scaffold_SMILE` field is also unique and cannot be null. The `log_p_min` and `log_p_max` fields are required.

.. _model_inference-table:

model_inference Table
~~~~~~~~~~~~~~~~~~~~~

The `model_inference` table stores information about model inference.

+------------------------+------------------+---------------------------------+
| Field                  | Type             | Description                     |
+========================+==================+=================================+
| id                     | UUID (Primary)   | Unique identifier.              |
+------------------------+------------------+---------------------------------+
| generated_SMILE_array  | Text[]           | An array of generated SMILEs.   |
|                        |                  | Limited to 500 elements.        |
+------------------------+------------------+---------------------------------+

.. note:: The `id` field is a primary key. The `generated_SMILE_array` field stores an array of generated SMILEs.

.. _generated_molecule-table:

generated_molecule Table
~~~~~~~~~~~~~~~~~~~~~~~~

The `generated_molecule` table stores information about generated molecules.

+-----------------+------------------+----------------------------------------+
| Field           | Type             | Description                            |
+=================+==================+========================================+
| SMILE           | Text (Primary)   | Unique identifier for the molecule.    |
+-----------------+------------------+----------------------------------------+
| predicted_log_p | Float (Not Null) | Predicted logP                         |
|                 |                  | value for the molecule.                |
+-----------------+------------------+----------------------------------------+

.. note:: The `SMILE` field is a primary key and must be unique. The `predicted_log_p`, `metadata`, `metadata2`, and `metadata3` fields are required.


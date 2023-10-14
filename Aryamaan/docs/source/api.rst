*************
API Endpoints
*************
.. _api-endpoints:

===========
Model Input
===========
.. _model-input:

Insert Model Input
~~~~~~~~~~~~~~~~~~~
..  http:post:: /insert_model_input

**Description:** Insert data into the database.

**Request:**

.. code-block:: json

    {
        "scaffold_smile": "C1CCC1",
        "log_p_min": 0.0,
        "log_p_max": 1.0
    }

**Response:**

- 200 OK: Data inserted successfully.
- 400 Bad Request: Malformed data.

List Model Inputs
~~~~~~~~~~~~~~~~~
..  http:get:: /list_model_inputs

**Description:** List model inputs.

**Response:**

- List of model inputs as JSON.

.. code-block:: json

    [
        {
            "id": "c18000c08234a682bf435d747b0e0bf6",
            "scaffold_smile": "C1CCC1",
            "log_p_min": 0.0,
            "log_p_max": 1.0
        },
        {
            "id": "45ecc4dd568420521f285e291939ecf1",
            "scaffold_smile": "C1CCCC1",
            "log_p_min": 0.5,
            "log_p_max": 1.5
        }
    ]

================
Model Inferences
================

Insert Model Inference
~~~~~~~~~~~~~~~~~~~~~~
..  http:post:: /insert_model_inference

**Description:** Insert model inference data into the database.

**Request:**

.. code-block:: json

    {
        "generated_smile_array": ["CCO", "C1CCCC1N", "O=C(C)Oc1ccccc1C(=O)O"],
        "id": "31a400c00b28c398f2490a8f05b6e168"
    }


**Response:**

- 200 OK: Data inserted successfully.

- 400 Bad Request: Malformed data.

List Model Inferences
~~~~~~~~~~~~~~~~~~~~~
..  http:get:: /list_model_inferences

**Description:** List model inferences.

**Response:**

- List of model inferences as JSON.

.. code-block:: json

    [
        {
            "id": "f510e45d2ae35afc90198f04430386b7",
            "generated_smile_array": ["C1CCC1", "CC(C)(C)C1CC1"]
        },
        {
            "id": "c48b89fc935a194a8b115712484fc331",
            "generated_smile_array": ["C1CCCC1", "CCOC(=O)C1CCCC1"]
        }
    ]

=======================
Generated Molecules
=======================

Insert Generated Molecule
~~~~~~~~~~~~~~~~~~~~~~~~~
..  http:post:: /insert_generated_molecule

**Description:** Insert generated molecule data into the database.

**Request:**

.. code-block:: json

    {
        "smile": "C1CCC1",
        "predicted_log_p": 0.8
    }

**Response:**

- 200 OK: Data inserted successfully.

- 400 Bad Request: Malformed data.


List Generated Molecules
~~~~~~~~~~~~~~~~~~~~~~~~~
..  http:get:: /list_generated_molecules

**Description:** List generated molecules.

**Response:**

- List of generated molecules as JSON.

.. code-block:: json

    [
        {
            "smile": "C1CCC1",
            "predicted_log_p": 0.8
        },
        {
            "smile": "CC(C)(C)C1CC1",
            "predicted_log_p": 1.2
        }
    ]

===================
Model Inference Run
===================

Run Model Inference
~~~~~~~~~~~~~~~~~~~
..  http:post:: /run_model_inference

**Description:** Run model inference with different model types.

**Request:**

.. code-block:: json

    {
        "model_type": "scaffold_constrained",
        "payload": {
            "scaffold_smile": "CC(C)(C(=O)O)c1ccc(cc1)C(O)CCCN2CCC(CC2)C(O)(*)c4ccccc4",
            "log_p_min": 5.1,
            "log_p_max": 5.9
        }
    }
**Response:**

- 200 OK: Successful response.

- 400 Bad Request: Malformed data.

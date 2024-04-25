# molecular-generative-modeling

This repository contains the code for the Molecular Generative Modeling Data Mine project. The project is a collaboration between the [Data Mine](https://datamine.purdue.edu/) and Merck & Co. The goal of the project is to develop a full-stack platform for molecular generative modeling. The platform will include a web-based user interface for data input and model training, a database for storing molecular data, and a machine learning model for generating new molecules.

## Project Structure
```
.
├── README.md - This file.
├── backend
│   ├── containers
│   │   ├── multiobj-rationale - Docker container for the multi-objective rationale aka Graph Translation model.
│   │   ├── scaffold_constrained - Docker container for the scaffold constrained aka SAMOA model.
│   │   ├── sqs_processor - Docker container for the SQS processor.
│   │   ├── utils - Utils package for filtering molecules.
│   │   └── vae_gan - Docker container for the VAE-GAN model.
│   └── lambdas
│       ├── enqueue-sqs.py - Lambda function for enqueuing messages to the SQS queue.
│       └── poll-molecule-db.py - Lambda function for polling the molecule database.
├── dev_models - Models still under development for the project, lacking documentation.
├── front-end - Front-end code for the project.
└── prop_pred
    └── p_chem_CEVR - Model to predict Log D.
```

### Backend
The backend of the project is divided into two main components: the containers and the lambdas. The containers are Docker containers that run the machine learning models. The lambdas are AWS Lambda functions that handle the communication between the containers and the database.

Each Container can be built locally with the `test_locally.sh` script. Additionally, the `push_to_aws.sh` script can be used to build containers and upload them to ECR and deploy them using Sagemaker Severless Inference Endpoints.

All containers are made with the [sagemaker-inference-toolkit](https://github.com/aws/sagemaker-inference-toolkit). Additionally, see [here](https://docs.aws.amazon.com/sagemaker/latest/dg/serverless-endpoints.html) for more details on how to deploy a container as a serverless endpoint.

Note that the deployed backend uses lambda functions as an intermediary between the front-end and the containers. The lambda functions are responsible for sending messages to the SQS queue and polling the molecule database. Following the documentation in the sagemaker-inference-toolkit repo, one can test the endpoints locally.

When deployed, dynamoDB is used to store the molecules. The `poll-molecule-db.py` lambda function is responsible for polling the database and sending the molecules to the appropriate container for processing.

See backend_diagram.pdf for a visual representation of the backend architecture.

### Front-end
The front-end of the project is a web-based user interface for data input and model training. The front-end is built using React and is hosted on AWS Amplify. The front-end code can be found in the `front-end` directory.

Use `npm install` to install the necessary dependencies and `npm run dev` to run the front-end locally.

A static export can be generated with `npm run build` and the output can be found in the `out` directory.

## Getting Started
First, clone the repository to your local machine:
```
git clone https://github.com/jag-prabhakaran/molecular-generative-modeling
``` 
Next, with [Docker](https://www.docker.com/) installed, run the following command to build the Docker image for each container:
```bash
./test_locally.sh
```

To deploy the containers to AWS, run the following command:
```bash
./push_to_aws.sh
```

To run the front-end locally, navigate to the `front-end` directory and run the following commands:
```bash
npm install
npm run dev
```
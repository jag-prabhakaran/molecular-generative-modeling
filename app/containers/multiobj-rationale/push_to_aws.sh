#!/bin/bash
algorithm_name=sagemaker-multiobj-rationale
endpoint_name=multiobj-rationale

echo "Building container and pushing to AWS..."

account=$(aws sts get-caller-identity --query Account --output text)

# Get the region defined in the current configuration (default to us-west-2 if none defined)
region=$(aws configure get region)
region=${region:-us-west-2}

fullname="${account}.dkr.ecr.${region}.amazonaws.com/${algorithm_name}:latest"

# If the repository doesn't exist in ECR, create it.
aws ecr describe-repositories --repository-names "${algorithm_name}" > /dev/null 2>&1

if [ $? -ne 0 ]
then
    aws ecr create-repository --repository-name "${algorithm_name}" > /dev/null
fi

# Get the login command from ECR and execute it directly
aws ecr get-login-password --region ${region}|docker login --username AWS --password-stdin ${fullname}

# Build the Docker image locally with the image name and then push it to ECR
# with the full name.

docker build -t ${algorithm_name} --platform linux/amd64 -f Dockerfile ..
docker tag ${algorithm_name} ${fullname}

docker push ${fullname}

aws sagemaker list-endpoints --output text | awk -v search="$endpoint_name" '$4 == search {exit 1}';
if [ $? -eq 1 ]; then
    echo "Endpoint already exists. Deleting it..."
    aws sagemaker delete-endpoint --endpoint-name $endpoint_name
    aws sagemaker wait endpoint-deleted --endpoint-name $endpoint_name
    echo "Creating new endpoint..."
    aws sagemaker create-endpoint --endpoint-name $endpoint_name --endpoint-config-name $endpoint_name
    echo "Waiting for endpoint to be in service..."
    aws sagemaker wait endpoint-in-service --endpoint-name $endpoint_name
else
    echo "Endpoint does not exist. Creating it..."
    aws sagemaker create-endpoint --endpoint-name $endpoint_name --endpoint-config-name $endpoint_name
    echo "Waiting for endpoint to be in service..."
    aws sagemaker wait endpoint-in-service --endpoint-name $endpoint_name
fi

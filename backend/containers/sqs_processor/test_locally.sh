#!/bin/bash
container_name="sqs-manager"
docker build -t $container_name --platform linux/amd64 -f Dockerfile . && docker run --platform linux/amd64 $container_name

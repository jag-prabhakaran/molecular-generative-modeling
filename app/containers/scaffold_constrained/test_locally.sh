#!/bin/bash
algorithm_name="scaffold-constrained"
docker build -t $algorithm_name --platform linux/amd64 -f Dockerfile .. && docker run --platform linux/amd64 $algorithm_name

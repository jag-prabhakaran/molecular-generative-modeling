FROM mambaorg/micromamba

# Create the environment:
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml "/tmp/"
RUN micromamba install -y -n base -f env.yml && \
    micromamba clean --all --yes

USER root

# Install java for MMS
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -y install default-jre-headless && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV JAVA_HOME /usr/lib/jvm/java-17-openjdk-amd64/
RUN export JAVA_HOME

# Copy required files:
COPY --chown=$MAMBA_USER:$MAMBA_USER ckpt "/tmp/model_files/ckpt"
COPY --chown=$MAMBA_USER:$MAMBA_USER ckpt "/opt/ml/model/ckpt"
COPY --chown=$MAMBA_USER:$MAMBA_USER fuseprop "/tmp/model_files/fuseprop"
COPY --chown=$MAMBA_USER:$MAMBA_USER decode.py "/tmp/model_files/decode.py"
COPY --chown=$MAMBA_USER:$MAMBA_USER generate.sh "/tmp/model_files/generate.sh"
COPY --chown=$MAMBA_USER:$MAMBA_USER sagemaker_model_handler.py "/tmp/model_files/sagemaker_model_handler.py"
COPY --chown=$MAMBA_USER:$MAMBA_USER rationale.txt "/tmp/model_files/rationale.txt"
COPY --chown=$MAMBA_USER:$MAMBA_USER entrypoint.py "/tmp/model_files/entrypoint.py"
#COPY --chown=$MAMBA_USER:$MAMBA_USER scaffold-constrained.mar "/tmp/model_files/multiobj-rationale.mar"

# Run the application:
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/tmp/model_files/entrypoint.py"]
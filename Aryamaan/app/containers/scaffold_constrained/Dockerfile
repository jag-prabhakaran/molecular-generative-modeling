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
COPY --chown=$MAMBA_USER:$MAMBA_USER data "/tmp/data/"
COPY --chown=$MAMBA_USER:$MAMBA_USER model_files "/tmp/model_files/"
COPY --chown=$MAMBA_USER:$MAMBA_USER data "/opt/ml/model/data/"
COPY --chown=$MAMBA_USER:$MAMBA_USER model_files "/opt/ml/model/model_files/"
COPY --chown=$MAMBA_USER:$MAMBA_USER entrypoint.py "/tmp/"
COPY --chown=$MAMBA_USER:$MAMBA_USER sagemaker_model_handler.py "/tmp/"
COPY scaffold-constrained.mar "/opt/ml/model/"

# Run the application:
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/tmp/entrypoint.py"]
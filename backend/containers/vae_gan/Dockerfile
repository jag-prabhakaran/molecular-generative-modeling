FROM mambaorg/micromamba

# Create the environment:
COPY --chown=$MAMBA_USER:$MAMBA_USER ./vae_gan/env.yml "/tmp/"
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
COPY --chown=$MAMBA_USER:$MAMBA_USER ./vae_gan/generator.py "/tmp/"
COPY --chown=$MAMBA_USER:$MAMBA_USER ./vae_gan/model_weights.h5 "/opt/ml/model/"
COPY --chown=$MAMBA_USER:$MAMBA_USER ./vae_gan/app.py "/tmp/"
COPY --chown=$MAMBA_USER:$MAMBA_USER ./vae_gan/sagemaker_model_handler.py "/tmp/"
COPY --chown=$MAMBA_USER:$MAMBA_USER ./vae_gan/entrypoint.py "/tmp/"
COPY --chown=$MAMBA_USER:$MAMBA_USER ./vae_gan/vae_gan.mar "/opt/ml/model/"
COPY --chown=$MAMBA_USER:$MAMBA_USER ./utils "/tmp/utils/"

ARG

# Run the application:
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/tmp/entrypoint.py"]
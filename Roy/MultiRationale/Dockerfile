FROM mambaorg/micromamba

# Create the environment:
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml "/tmp/"
RUN micromamba install -y -n base -f env.yml && \
    micromamba clean --all --yes

# Copy required files:
COPY --chown=$MAMBA_USER:$MAMBA_USER multiobj-rationale "/tmp/model_files/"

# Install Files
apt-get install vim


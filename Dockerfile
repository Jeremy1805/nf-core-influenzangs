<<<<<<< HEAD
FROM nfcore/base
LABEL authors="Jeremy Guntoro, Aengus Stewart" \
      description="Docker image containing all requirements for nf-core/influenzangs pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-influenzangs-1.0dev/bin:$PATH
=======
FROM nfcore/base:1.9
LABEL authors="Jeremy Guntoro, Aengus Stewart" \
      description="Docker image containing all software requirements for the nf-core/influenzangs pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-influenzangs-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-influenzangs-1.0dev > nf-core-influenzangs-1.0dev.yml
>>>>>>> TEMPLATE

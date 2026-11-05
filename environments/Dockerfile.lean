# Dockerfile.lean -- minimal QIIME 2 environment for the classification step only.
#
# Reproduces the environment measured in Section 4.6 of the manuscript: 5 plugins
# against the 27 in the full distribution, which removes the plugin-import
# component of the fixed per-invocation cost.
#
#   docker build --platform=linux/amd64 \
#       -f environments/Dockerfile.lean -t qiime2-lean:2026.7 .
#   docker run --rm -v "$PWD":/data qiime2-lean:2026.7 \
#       qiime feature-classifier classify-sklearn \
#         --i-classifier /data/classifier.qza \
#         --i-reads /data/rep-seqs.qza \
#         --o-classification /data/taxonomy.qza
#
# This image CANNOT denoise. q2-dada2 is deliberately absent. Use the official
# QIIME 2 distribution for denoising and this image for classification only.
#
# Build for linux/amd64: QIIME 2 2026.7 publishes linux-64 builds, so pass
# --platform=linux/amd64 on the build command (shown above) rather than
# pinning it here. The scikit-learn pin must match the version the classifier
# artifact was serialised under or it is refused at load time.
FROM condaforge/miniforge3:25.9.1-0

LABEL org.opencontainers.image.title="qiime2-lean" \
      org.opencontainers.image.version="2026.7" \
      org.opencontainers.image.description="Minimal QIIME 2 environment for classify-sklearn"

COPY environments/qiime2-lean.yml /tmp/qiime2-lean.yml

RUN conda env create -y -f /tmp/qiime2-lean.yml \
 && conda clean -afy \
 && rm /tmp/qiime2-lean.yml

ENV PATH=/opt/conda/envs/qiime2-lean-2026.7/bin:$PATH
ENV MPLBACKEND=Agg

# Single-threaded BLAS. Section 4.7 shows classification is slower with more
# workers at this scale, so the default is one.
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

WORKDIR /data
ENTRYPOINT []
CMD ["qiime", "info"]

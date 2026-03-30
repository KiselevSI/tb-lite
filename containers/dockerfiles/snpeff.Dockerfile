# snpeff -- SnpEff variant annotation from local installation
# Build context: project root (tb-lite/)

FROM eclipse-temurin:21-jre-jammy

LABEL maintainer="TB-Lite pipeline"

COPY assets/h37rv.fa /ref/h37rv.fa
COPY snpEff_latest_core/snpEff /opt/snpEff

RUN apt-get update && apt-get install -y --no-install-recommends \
        bash ca-certificates tabix bcftools \
    && rm -rf /var/lib/apt/lists/* \
    && printf '#!/bin/bash\nexec java -jar /opt/snpEff/snpEff.jar -c /opt/snpEff/snpEff.config "$@"\n' > /usr/local/bin/snpEff \
    && chmod +x /usr/local/bin/snpEff

ENV SNPEFF_HOME=/opt/snpEff
ENV PATH=/usr/local/bin:$PATH

CMD ["snpEff"]

FROM r-base:3.6.0

RUN mkdir /app

COPY dataframes.RData /app
COPY workflow.R /app

WORKDIR /app

CMD ["Rscript", "workflow.R"]

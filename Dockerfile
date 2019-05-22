FROM r-base:3.6.0

RUN mkdir /app

COPY useful_files/dataframes.RData /app
COPY useful_files/workflow.R /app

WORKDIR /app

CMD ["Rscript", "workflow.R"]

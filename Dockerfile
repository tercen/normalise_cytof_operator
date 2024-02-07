FROM tercen/runtime-flowsuite:3.15-4

COPY . /operator
WORKDIR /operator

RUN R -e "install.packages('tercen', repos = c(TERCEN='https://cran.tercen.com/api/v1/rlib/tercen', CRAN='https://cran.tercen.com/api/v1/rlib/CRAN'))"

ENV TERCEN_SERVICE_URI https://tercen.com

ENTRYPOINT ["R", "--no-save", "--no-restore", "--no-environ", "--slave", "-f", "main.R", "--args"]
CMD ["--taskId", "someid", "--serviceUri", "https://tercen.com", "--token", "sometoken"]
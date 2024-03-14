FROM golang
COPY DNAAnalysis/* /opt/DNAAnalysis/
COPY *.go /opt/
COPY go.mod /opt/
WORKDIR /opt/
RUN go get
ENTRYPOINT ["go", "run", "."]
FROM python:3.10

LABEL author="Rachel Duffin" \
      description="seglh-polyedge:v1.1.0" \
      maintainer="rachel.duffin2@nhs.net"

RUN mkdir -p /polyedge /data /outputs
ADD ./ /polyedge
RUN apt-get update -y && apt-get install wkhtmltopdf -y
RUN pip3 install -r /polyedge/package-requirements.txt
RUN ln /polyedge/polyedge.py /usr/local/bin/polyedge.py
WORKDIR /outputs
ENTRYPOINT [ "python3","/polyedge/polyedge.py" ]

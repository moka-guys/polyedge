FROM python:3.10

LABEL author="Rachel Duffin" \
      description="seglh-polyedge v1.1.1.0" \
      maintainer="rachel.duffin2@nhs.net"

RUN mkdir -p /polyedge
ADD ./ /polyedge
RUN pip3 install -r /polyedge/package-requirements.txt
# RUN ln /polyedge/polyedge.py /usr/local/bin/polyedge.py
ENTRYPOINT [ "python3","/polyedge/polyedge.py" ]

FROM python:3.10
RUN mkdir /code /data
ADD ./* /code/
RUN pip3 install -r /code/package-requirements.txt
ENTRYPOINT [ "python","/code/polyedge.py" ]


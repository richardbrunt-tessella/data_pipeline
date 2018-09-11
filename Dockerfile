FROM python:2.7
# install fresh these requirements.
# do this before copying the code to minimize image layer rebuild for dev
COPY ./requirements.txt /usr/src/app/
RUN pip install --no-cache-dir -r /usr/src/app/requirements.txt

#put the application in the right place
WORKDIR /usr/src/app
COPY . /usr/src/app

#To ensure that docker stop will signal any long running ENTRYPOINT 
#executable correctly, you need to remember to start it with exec
ENTRYPOINT [ "exec", "python", "-m", "mrtarget.CommandLine" ]

CMD [ "--dry-run", "master" ]

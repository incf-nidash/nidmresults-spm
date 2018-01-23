FROM mtmiller/octave-snapshot

ADD code /code

RUN apt-get update && apt-get install -y git

RUN cd /code && git clone https://github.com/spm/spm12.git && cd spm12/src && make PLATFORM=octave && make install PLATFORM=octave && cd ../../..

RUN rm /code/spm12/@file_array/private/file2mat.m

RUN octave --no-window-system --eval "addpath('/code/spm12'); path; savepath;"
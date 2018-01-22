FROM mtmiller/octave-snapshot

RUN apt-get update && apt-get install -y git

RUN git clone https://github.com/spm/spm12.git && cd spm12/src && make PLATFORM=octave && make install PLATFORM=octave
#!/bin/bash

python -c "import alf; alf.initialize(engine='pycharmm')"
python -c "import alf; alf.runflat(1,50,13000,39000,engine='pycharmm')"
python -c "import alf; alf.runflat(51,60,125000,375000,engine='pycharmm')"

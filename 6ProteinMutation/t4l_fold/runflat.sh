#!/bin/bash

python -c "import alf; alf.initialize(engine='bladelib')"
python -c "import alf; alf.runflat(1,100,13000,39000,engine='bladelib')"
python -c "import alf; alf.runflat(101,110,125000,375000,engine='bladelib')"

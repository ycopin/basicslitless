.. $Id: script.txt,v 1 2020/07/26 11:34:46 saljibouri $

.. _script:

runSimulation (script)
==============================

**Usage:** `runSimulation.py [options]`

Simulates the spectrogram for a given star based on the given parameters. The simulation is a basic geometric model of slitless spectroscopy, using circle for the order 0 and rectangles to represent the other orders. The simulator also searches for the best angle to use for the dispersion grating used in the spectrograph and will show a graph of the total contamination for each dispersion angle.

**Options:**
  --version             show program's version number and exit
  -h, --help            show this help message and exit

positional arguments:

  config                Configuration file path

  star                  Star identifier

optional arguments:
  -s SEEING, --seeing SEEING
                        Seeing in arcsec, default seeing is 1 arcsec
  -o ORDERS, --orders ORDERS
                        Orders to take into account as a comma separated list inside parenthesis
                        (default : (-1,0,1))

.. figure:: Sim.*
   :align: center
   :width: 80%

   Graphical output of :command:`runSimulation.py ../config/auxtel.ini HD107696`.

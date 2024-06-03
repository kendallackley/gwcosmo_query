First tried this on the LIGO cluster, but I can't set up a database there, so I'll set it up locally


curl -s -O http://elysium.elte.hu/~dalyag/GLADE+.txt

python version: /cvmfs/software.igwn.org/conda/envs/igwn-py39/bin/python
directory: /home/kendall.ackley/projects/gals/query_methods
GLADE+ = home_dir + ../GLADE+.txt

git clone git@git.ligo.org:lscsoft/gwcosmo.git
./gwcosmo/scripts_galaxy_catalogs/GLADE+/select_columnts.sh < edited to point to where GLADE+ lives
./gwcosmo/scripts_galaxy_catalogs/GLADE+/create_glade+.py


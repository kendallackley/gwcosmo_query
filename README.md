

# Repository Name

This repository contains astro tools for querying GW cosmology data.

## Prerequisites

Before you can use this repository, you need to have the following installed:

- [Poetry](https://python-poetry.org/): A dependency management and packaging tool for Python.
- [PostgreSQL 15](https://www.postgresql.org/): A powerful, open-source object-relational database system.

# Install PostgreSQL and create the database
Installation tutorial: https://www.postgresql.org/docs/current/tutorial-install.html

On MacOS:
brew install postgresql@15
brew services start postgresql@15
createdb gwcosmo_db

# Add gist extensions to database
psql -d gwcosmo_db -c "CREATE EXTENSION IF NOT EXISTS btree_gist;"

# Catalogs
- [GLADE+ Catalog](http://elysium.elte.hu/~dalyag/GLADE+.txt)

See RunScripts.ipynb for a tutorial on how to run the full script.

Notes: 
- This example to set up and make use of a PostgreSQL database are meant for local installation. The databases cannot be created on the LIGO Data Grid.
- If using ARM64 architecture (MacOS M1/M2 chips), the [Q3C extension](https://github.com/segasai/q3c) cannot be applied to the database. Work-arounds are planned until future versions of PostgreSQL are developed to include more complete support for the ARM64 architecture.
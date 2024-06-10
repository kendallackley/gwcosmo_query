import pandas as pd
import sqlalchemy as sa  # type: ignore
from ..classes import Base, GalaxySim
from ..helpers import glade_filter


def read_in_catalog(filename="catalogs/GLADE+_reduced.txt"):
    columns = [
        "glade",
        "hyperleda",
        "wise",
        "quasar_cluster_flag",
        "ra",
        "dec",
        "m_B",
        "Bflag",
        "m_K",
        "m_W1",
        "W1_flag",
        "redshift_cmb",
        "redshift_corrected_pec_vel",
        "error_peculiar_velocity",
        "redshift_lum_distance_flag",
    ]

    glade_chunk = pd.read_csv(filename, chunksize=5000, names=columns, sep=r"\s+")
    return glade_chunk


def add_galaxies_sim(df, engine, table):

    rows = df.rename(
        columns={
            "glade": "id",
            "coord": "hpx",
            "redshift_cmb": "redshift",
            "zerror": "redshift_error",
            "m_B": "B_mag",
            "m_K": "K_mag",
            "m_W1": "W1_mag",
        }
    ).to_dict("records")

    with sa.orm.Session(engine) as session:
        session.bulk_insert_mappings(table, rows)
        session.commit()
    return


def create_fake_catalogs() -> None:

    engine = sa.create_engine("postgresql://localhost/gwcosmo_db")
    Base.metadata.create_all(engine)
    glade_chunk = read_in_catalog()

    for chunk in glade_chunk:
        chunk_filtered = glade_filter(chunk, radians=True)
        add_galaxies_sim(chunk_filtered, engine, GalaxySim)

    with sa.orm.Session(engine) as session:
        session.execute(sa.text("ANALYZE"))

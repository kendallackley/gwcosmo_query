import numpy as np
import healpy as hp
from astropy.table import Table
from ..classes import Skymap, SkymapTile
from ..helpers import nest2uniq

import sqlalchemy as sa  # type: ignore
from sqlalchemy.exc import IntegrityError  # type: ignore


def make_skymap_tiles(nside, skymap_probability=None, gaussian=False):
    ipix = np.arange(hp.nside2npix(nside), dtype=np.uint64)
    uniq = nest2uniq(nside, ipix)
    skymap_probability_density = skymap_probability / np.sum(
        skymap_probability * hp.nside2pixarea(nside)
    )

    moc_data = Table(
        [uniq, skymap_probability_density], names=("UNIQ", "PROBDENSITY")
    )
    if gaussian:
        # Reduce size of skymap (after Gaussian smoothing) for tile generation
        moc_data = moc_data[moc_data["PROBDENSITY"] > 1e-10]

    # Bug that pandas doesn't keep the type of the column through the apply
    # So doing it by list comprehension instead
    # TODO: Fix this
    tiles = [
        SkymapTile(
            hpx=int(row["UNIQ"]), probdensity=np.float64(row["PROBDENSITY"])
        )
        for row in moc_data
    ]

    return tiles


def create_uniform_skymap(nside=1024, db_id=None):

    # Make fake Uniform skymap
    skymap = hp.ma(np.random.random(hp.nside2npix(nside))).data
    skymap /= np.sum(skymap)
    tiles = make_skymap_tiles(nside, skymap)

    engine = sa.create_engine("postgresql://localhost/gwcosmo_db")
    with sa.orm.Session(engine) as session:
        if db_id is None:
            # Get the maximum id in the Skymap table
            max_id = session.query(sa.func.max(Skymap.id)).scalar()
            if max_id is None:
                max_id = 0
            # Take the next number
            db_id = max_id + 1
        session.add(Skymap(id=db_id, tiles=tiles))
        session.commit()

        session.execute(sa.text("ANALYZE"))


def generate_skymap(nside, ra, dec, degrees=True):

    if degrees:
        ra = np.radians(ra)
        dec = np.radians(dec)
    theta = np.pi / 2 - dec
    phi = ra
    # Keep in RING ordering
    center = hp.ang2pix(nside, theta, phi)

    npix = hp.nside2npix(nside)
    skymap = np.zeros(npix)
    skymap[center] = 1
    skymap = hp.sphtfunc.smoothing(skymap, sigma=np.pi / 64)
    skymap /= np.sum(skymap)
    skymap = hp.pixelfunc.reorder(skymap, r2n=True)
    return skymap


def create_gaussian_skymap(ra, dec, nside=1024, db_id=None, degrees=True):

    # Make fake Gaussian skymap
    skymap = generate_skymap(nside, ra, dec, degrees=degrees)
    tiles = make_skymap_tiles(nside, skymap, gaussian=True)

    engine = sa.create_engine("postgresql://localhost/gwcosmo_db")
    with sa.orm.Session(engine) as session:
        if db_id is None:
            # Get the maximum id in the Skymap table
            max_id = session.query(sa.func.max(Skymap.id)).scalar()
            if max_id is None:
                max_id = 0
            db_id = max_id + 1
        while True:
            try:
                print("Inserting skymap with id:", db_id)
                session.add(Skymap(id=db_id, tiles=tiles))
                session.commit()
                break
            except IntegrityError as e:
                session.rollback()  # Rollback the session to a clean state
                # If it fails, get the maximum id in the Skymap table and add 1 to it
                max_id = session.query(sa.func.max(Skymap.id)).scalar()
                db_id = max_id + 1
                print(
                    "Skymap with id already exists in db. Trying again with db_id:",
                    db_id,
                )

        session.execute(sa.text("ANALYZE"))

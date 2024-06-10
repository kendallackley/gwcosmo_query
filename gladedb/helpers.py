import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy_healpix import HEALPix
import time
import sqlalchemy as sa  # type: ignore
from .classes import Galaxy, SkymapTile


def nest2uniq(nside, ipix):
    """
    https://gitlab.com/burstcube/mhealpy/-/blob/master/mhealpy/pixelfunc/moc.py
    """
    ipix = np.array(ipix)
    return 4 * nside * nside + ipix


def get_healpix_index_from_coord(ra, dec, nside, lonlat=False):
    """
    Arguments:
    ra: float, right ascension
    dec: float, declination
    nside: int, healpix nside
    lonlat: bool, if True, coordinates assumed to be in degrees,
            otherwise, radians. Default is False

    Returns:
    value: int, healpix index
    """
    hp = HEALPix(nside=nside, order="nested", frame="icrs")
    coord = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    if lonlat:
        return hp.lonlat_to_healpix(coord.ra, coord.dec)
    return hp.skycoord_to_healpix(coord)


def get_healpix_range_from_coord(ra, dec, nside, max_level=29):
    """
    Arguments:
    ra: float, right ascension in degrees
    dec: float, declination in degrees
    nside: int, healpix nside
    max_level: int, maximum level of the healpix grid. Default is 29

    Returns:
    value: tuple, (start, end) healpix range at the maximum level
    """
    ipix = get_healpix_index_from_coord(ra, dec, nside=nside, lonlat=True)
    level = np.log2(nside).astype(np.int64)
    shift = 2 * (max_level - level)
    value = (ipix << shift, (ipix + 1) << shift)
    return value


def get_healpix_max_index_from_coord(ra, dec, nside, max_level=29):
    """
    Arguments:
    ra: float, right ascension in degrees
    dec: float, declination in degrees
    nside: int, healpix nside
    max_level: int, maximum level of the healpix grid. Default is 29

    Returns:
    value: int, healpix index at the maximum level
    """
    ipix = get_healpix_index_from_coord(ra, dec, nside=nside, lonlat=True)
    level = np.log2(nside).astype(np.int64)
    shift = 2 * (max_level - level)
    value = ipix << shift
    return value


def find_galaxies_in_range(ra, dec, engine, GalaxyClass, nside=1024):
    """
    Find galaxies in a healpix range given by a coordinate

    Arguments:
    ra: float, right ascension in degrees
    dec: float, declination in degrees
    engine: sqlalchemy engine
    GalaxyClass: sqlalchemy class. Can be Galaxy or GalaxySim
    nside: int, healpix nside. Default is 1024

    Returns:
    galaxy_ids: list, list of galaxy ids
    """

    hpx_range = get_healpix_range_from_coord(ra, dec, nside=nside)

    with sa.orm.Session(engine) as session:
        start_time = time.time()
        session.execute(sa.text("ANALYZE"))
        query = session.query(GalaxyClass.id).filter(
            sa.func.int8range(int(hpx_range[0]), int(hpx_range[1]), "[)").op("@>")(
                GalaxyClass.hpx
            ),
        )

        galaxy_ids = query.all()
        end_time = time.time()

    print(
        f"Query run time: {np.round(end_time - start_time, 3)} seconds, "
        f"Number of galaxies: {len(galaxy_ids)}, nside: {nside}"
    )
    return galaxy_ids


def find_galaxies_in_skymap(
    skymap_id, engine, contour_level=None, GalaxyClass=Galaxy, SkymapTileClass=SkymapTile
):
    """
    Find the galaxies in a skymap given by a coordinate

    Arguments:
    skymap_id: int, skymap id
    engine: sqlalchemy engine
    contour_level: float, probability contour level. Default is None. Can take values between 0 and 1.
    GalaxyClass: sqlalchemy class. Can be Galaxy or GalaxySim
    SkymapTileClass: sqlalchemy class. Can be SkymapTile or SkymapTileSim

    Returns:
    galaxy_ids: list, list of galaxy ids
    """

    # hpx_index = get_healpix_index_from_coord(ra, dec, nside=nside, lonlat=True)
    with sa.orm.Session(engine) as session:
        start_time = time.time()
        session.execute(sa.text("ANALYZE"))
        min_probdensity = 1e-10
        if contour_level is not None:
            # Syntax adapted from healpix-alchemy example
            # https://github.com/skyportal/healpix-alchemy/
            cumulative_prob = sa.func.sum(
                    SkymapTileClass.probdensity * SkymapTileClass.hpx.area
                ).over(
                    order_by=SkymapTile.probdensity.desc()
                ).label(
                    'cumulative_prob'
                )
            subquery = sa.select(
                    SkymapTileClass.probdensity,
                    cumulative_prob
                ).filter(
                    SkymapTileClass.id == skymap_id
                ).subquery()
            min_probdensity = sa.select(
                    sa.func.min(subquery.columns.probdensity)
                ).filter(
                    subquery.columns.cumulative_prob <= contour_level
                ).scalar_subquery()
        query = sa.select(GalaxyClass.id, SkymapTileClass.hpx).filter(
            # SkymapTileClass.hpx.contains(hpx_index),
            SkymapTileClass.hpx.contains(GalaxyClass.hpx),
            SkymapTileClass.id == skymap_id,
            SkymapTileClass.probdensity >= min_probdensity,
        )

        result = session.execute(query).fetchall()
        end_time = time.time()

    galaxy_ids, skymaptile_hpxs = zip(*result)
    unique_galaxy_ids = set(galaxy_ids)
    unique_skymaptile_hpxs = set(skymaptile_hpxs)
    print(
        f"Query run time: {np.round(end_time - start_time, 3)} seconds, "
        f"Number of galaxies: {len(unique_galaxy_ids)} over {len(unique_skymaptile_hpxs)} healpix indexes"
    )
    return galaxy_ids


def glade_filter(catalog_chunk, radians=False):
    """
    Filter the GLADE catalog
    Based on scripts from https://git.ligo.org/lscsoft/gwcosmo/-/tree/master/scripts_galaxy_catalogs/GLADE+

    Arguments:
    catalog_chunk: pandas dataframe, GLADE catalog
    radians: bool, if True, the radian values of ra and dec should be used as the new coordinates,
            otherwise, use ra and dec in degrees. Default is False.

    Returns:
    catalog_chunk: pandas dataframe, filtered GLADE catalog
    """
    # Remove galaxies with no reported redshift
    catalog_chunk = catalog_chunk[catalog_chunk["redshift_cmb"].notnull()]

    # Remove QSOs and galaxy clusters.
    catalog_chunk = catalog_chunk[catalog_chunk["quasar_cluster_flag"] != "C"]
    catalog_chunk = catalog_chunk[catalog_chunk["quasar_cluster_flag"] != "Q"]

    # Make sure redshifts are positive
    catalog_chunk = catalog_chunk[catalog_chunk["redshift_cmb"] > 0]

    # Remove galaxies without peculiar velocity corrections
    # unless the galaxy is above redshift 0.05
    mask = (catalog_chunk["redshift_corrected_pec_vel"] == 0) & (
        catalog_chunk["redshift_cmb"] <= 0.05
    )
    catalog_chunk = catalog_chunk[~mask]

    # Remove galaxies without measured redshift and distance (0) and
    # with redshift from distance measure (2)
    catalog_chunk = catalog_chunk[
        (catalog_chunk["redshift_lum_distance_flag"] != 0)
        | (catalog_chunk["redshift_lum_distance_flag"] != 2)
    ]

    catalog_chunk["zerror"] = np.zeros(len(catalog_chunk))
    catalog_chunk["error_peculiar_velocity"] = catalog_chunk[
        "error_peculiar_velocity"
    ].fillna(0)

    # Redshift error for galaxies with spectroscopic errors
    # This is an absolute error
    zerror_sp = 1.5e-4
    zs_sp = catalog_chunk["redshift_lum_distance_flag"] == 3
    catalog_chunk.loc[zs_sp, "zerror"] = np.sqrt(
        zerror_sp**2 + catalog_chunk.loc[zs_sp, "error_peculiar_velocity"] ** 2
    )

    # Z error for galaxies with WISE,
    # This is a relative error of 4% (From Maciej discussion)
    zerror_wise = 0.04
    zs_wise = catalog_chunk["wise"].notnull()
    z_wise = catalog_chunk.loc[zs_wise, "redshift_cmb"]
    catalog_chunk.loc[zs_wise, "zerror"] = np.sqrt(
        (zerror_wise * (1 + z_wise)) ** 2
        + catalog_chunk.loc[zs_wise, "error_peculiar_velocity"] ** 2
    )

    # Z error for galaxies with photoz
    # This is an absolute error
    zerror_ph = 1.5e-2
    zs_ph = catalog_chunk["redshift_lum_distance_flag"] == 1
    catalog_chunk.loc[zs_ph, "zerror"] = np.sqrt(
        zerror_ph**2 + catalog_chunk.loc[zs_ph, "error_peculiar_velocity"] ** 2
    )

    # Z error for galaxies from Hyperleda
    # This is a relative error.
    # Since these are small redshift the relative error is multiplied only to z
    zerror_hyperleda = 0.36
    zs_hyperleda = catalog_chunk["hyperleda"].notnull()
    z_hyperleda = catalog_chunk.loc[zs_hyperleda, "redshift_cmb"]
    catalog_chunk.loc[zs_hyperleda, "zerror"] = np.sqrt(
        (zerror_hyperleda * (1 + z_hyperleda)) ** 2
        + catalog_chunk.loc[zs_hyperleda, "error_peculiar_velocity"] ** 2
    )

    catalog_chunk["m_B"] = catalog_chunk["m_B"].fillna(0).replace(0, np.nan)
    catalog_chunk["m_K"] = catalog_chunk["m_K"].fillna(0).replace(0, np.nan)
    catalog_chunk["m_W1"] = catalog_chunk["m_W1"].fillna(0).replace(0, np.nan)

    if radians:
        catalog_chunk["ra_rad"] = np.radians(catalog_chunk["ra"].values)
        catalog_chunk["dec_rad"] = np.radians(catalog_chunk["dec"].values)
        catalog_chunk["coord"] = catalog_chunk.apply(
            lambda row: SkyCoord(
                ra=row["ra_rad"] * u.degree, dec=row["dec_rad"] * u.degree
            ),
            axis=1,
        )
        return catalog_chunk[
            [
                "glade",
                "coord",
                "redshift_cmb",
                "zerror",
                "m_B",
                "m_K",
                "m_W1",
                "ra",
                "dec",
                "ra_rad",
                "dec_rad",
            ]
        ]

    catalog_chunk["coord"] = catalog_chunk.apply(
        lambda row: SkyCoord(ra=row["ra"] * u.degree, dec=row["dec"] * u.degree), axis=1
    )
    return catalog_chunk[
        ["glade", "coord", "redshift_cmb", "zerror", "m_B", "m_K", "m_W1", "ra", "dec"]
    ]

from sqlalchemy import orm  # type: ignore
import sqlalchemy as sa  # type: ignore
import healpix_alchemy as ha  # type: ignore
from sqlalchemy import DDL, event  # type: ignore
from sqlalchemy.ext.declarative import declared_attr, as_declarative  # type: ignore


def gist_index(index_type):
    """Class decorator to add a GiST or SP-GiST index to a table."""

    def decorator(cls):
        event.listen(
            cls.__table__,
            "after_create",
            DDL(
                f"CREATE INDEX ix_{cls.__tablename__}_hpx ON {cls.__tablename__} USING {index_type} (hpx);"
            ),
        )
        event.listen(
            cls.__table__,
            "after_create",
            DDL(
                f"CREATE INDEX ix_{cls.__tablename__}_id_hpx ON {cls.__tablename__} USING gist (id gist_int4_ops, hpx);"
            ),
        )
        return cls

    return decorator


@as_declarative()
class Base:

    @declared_attr
    def __tablename__(cls):
        return cls.__name__.lower()


@gist_index("gist")
class Galaxy(Base):
    """
    Each row of the Galaxy table represents a point in a catalog
    """

    id = sa.Column(sa.Integer, primary_key=True)
    hpx = sa.Column(ha.Point, index=True, nullable=False)
    redshift = sa.Column(sa.Float, index=True, nullable=True)
    redshift_error = sa.Column(sa.Float, index=True, nullable=True)
    B_mag = sa.Column(sa.Float, index=True, nullable=True)
    K_mag = sa.Column(sa.Float, index=True, nullable=True)
    W1_mag = sa.Column(sa.Float, index=True, nullable=True)
    ra = sa.Column(sa.Float, nullable=False)
    dec = sa.Column(sa.Float, nullable=False)


class Skymap(Base):
    """
    Each row of the Skymap table represents a LIGO/Virgo
    HEALPix localization map.
    """

    id = sa.Column(sa.Integer, primary_key=True)
    tiles = orm.relationship(lambda: SkymapTile)


@gist_index("spgist")
class SkymapTile(Base):
    """
    Each row of the SkymapTile table represents a multi-resolution
    HEALPix tile within a LIGO/Virgo localization map. There is a
    one-to-many mapping between Skymap and SkymapTile.
    """

    id = sa.Column(sa.ForeignKey(Skymap.id), primary_key=True)
    hpx = sa.Column(ha.Tile, primary_key=True)
    probdensity = sa.Column(sa.Float, nullable=False)


@gist_index("gist")
class GalaxySim(Base):
    """
    Each row of the GalaxySim table represents a point in a catalog
    This is for testing purposes only
    """

    id = sa.Column(sa.Integer, primary_key=True)
    hpx = sa.Column(ha.Point, index=True, nullable=False)
    redshift = sa.Column(sa.Float, index=True, nullable=True)
    redshift_error = sa.Column(sa.Float, index=True, nullable=True)
    B_mag = sa.Column(sa.Float, index=True, nullable=True)
    K_mag = sa.Column(sa.Float, index=True, nullable=True)
    W1_mag = sa.Column(sa.Float, index=True, nullable=True)
    ra = sa.Column(sa.Float, nullable=False)
    dec = sa.Column(sa.Float, nullable=False)


class SkymapSim(Base):
    """
    Each row of the Skymap table represents a LIGO/Virgo
    HEALPix localization map.
    """

    id = sa.Column(sa.Integer, primary_key=True)
    tiles = orm.relationship(lambda: SkymapTileSim)


@gist_index("spgist")
class SkymapTileSim(Base):
    """
    Each row of the SkymapTile table represents a multi-resolution
    HEALPix tile within a LIGO/Virgo localization map. There is a
    one-to-many mapping between Skymap and SkymapTile.
    """

    id = sa.Column(sa.ForeignKey(SkymapSim.id), primary_key=True)
    hpx = sa.Column(ha.Tile, primary_key=True)
    probdensity = sa.Column(sa.Float, nullable=False)

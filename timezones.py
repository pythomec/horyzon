"""
Various routines related to time zones
"""

import warnings
import logging
import datetime

def lonlat2tz(lon, lat):
    """Get time zone for a particular point on Earth.

    TBD. Currently it prints warning and returns UTC.

    :param lon: longitude
    :param lat: latitude
    :return:
    """

    #TODO: implement, e.g. using https://github.com/MrMinimal64/timezonefinder/tree/master/timezonefinder

    logger = logging.getLogger(__name__)
    warning = 'lonlat2tz not implemented, default time for all locations is UTC'
    logger.warning(warning)
    warnings.warn(warning)

    return datetime.timezone.utc


def fill_default_tz(time_or_date, lon=None, lat=None):
    """Fill default time zone into date/datetime object

    Date is converted to datetime, time: 00:00

    If lol/lat coordinates are provided, function lonlat2tz is used to look up default time zone for that
    location, otherwise UTC is used

    :param time_or_date: datetime or date
    :param lon: longitude
    :param lat: latitude
    :return:
    """

    # converte date to datetime
    if isinstance(time_or_date, datetime.date):
        time_or_date =  datetime.datetime(time_or_date.year, time_or_date.month, time_or_date.day)

    # set default datetime if missing
    if time_or_date.tzinfo is None:
        time_or_date = time_or_date.replace(tzinfo = get_default_tz(lon=lon, lat=lat))

    return time_or_date


def get_default_tz(lon=None, lat=None, tzinfo=None):
    """Return default time zone for given position or current tzinfo if not None

    Returns UTC if lon or lat is not provided.

    :param lon:
    :param lat:
    :param tzinfo: current tzinfo
    :return:
    """

    if tzinfo is not None:
        return tzinfo

    if lon is not None and lat is not None:
        tzinfo = lonlat2tz(lon=lon, lat=lat)
    else:
        tzinfo = datetime.timezone.utc

    return tzinfo

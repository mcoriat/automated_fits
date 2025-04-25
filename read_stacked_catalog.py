from astropy.io import fits

def read_stacked_catalog(catalog_file, srcid_ref):
    """
    Read a stacked catalog and return a dict { srcid_ref: [obsid1, obsid2, ...] }.
    Stops reading once you've passed all entries for that srcid.
    """
    with fits.open(catalog_file) as hdul:
        data = hdul[1].data

    mapping = {}
    found = False

    for row in data:
        sid   = row['SRCID']
        obsid = row['OBS_ID']
        # only care about this one srcid
        if sid == srcid_ref:
            mapping.setdefault(sid, []).append(obsid)
            found = True
        elif found:
            # once we've seen all of srcid_ref's rows, quit
            break

    return mapping

